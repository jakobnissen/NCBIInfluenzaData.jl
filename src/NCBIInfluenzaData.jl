"""
    NCBIInfluenzaData

Package for parsing and filtering NCBI's bulk influenza data. This can be used to
extract and work with tens of thousands of sequences in bulk.

It has been confirmed to work with the NCBI data of 2020-10-13
"""
module NCBIInfluenzaData

using BioSequences
using CodecZlib
using ErrorTypes
using FASTX 
using Serialization
using Transducers
using UnicodePlots
using Influenza
import Downloads

import Influenza: ReferenceProtein

function try_split!(
    v::Vector{SubString{String}},
    s::Union{String, SubString{String}},
    sep::UInt8
)::Option{UInt}
    n = UInt(0)
    start = 1
    isempty(v) && return none
    @inbounds for i in 1:ncodeunits(s)
        if codeunit(s, i) == sep
            n += 1
            n >= length(v) && return none
            substr = SubString(s, start, i-1)
            @inbounds v[n] = substr
            start = i + 1
        end
    end
    @inbounds v[n+1] = SubString(s, start, ncodeunits(s))
    some(n + 1)
end

function try_split!(
    v::Vector{SubString{String}},
    s::Union{String, SubString{String}},
    sep::Char
)::Option{UInt}
    try_split!(v, s, UInt8(sep))
end

# This is pretty primitive, but fuck it
function maybe_gzip(f, path::AbstractString)
    if endswith(path, ".gz")
        io = GzipDecompressorStream(open(path))
        try
            f(io)
        finally
            close(io)
        end
    else
        open(f, path)
    end
end

"""
    download_influenza_data(dst::AbstractString, force::Bool=false)

Download data from the NCBI FTP server to an existing directory `dst`.
If `force`, re-download, even if the files already exist in the directory.
"""
function download_influenza_data(dst_dir::AbstractString; force::Bool=false)
    isdir(dst_dir) || error("Path is not a directory: \"$dst_dir\"")
    ftp_address = "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/"
    for filename in [
        "genomeset.dat.gz",
        "influenza.dat.gz",
        "influenza.fna.gz",
        "influenza_aa.dat.gz"
        ]
        dstpath = joinpath(dst_dir, filename)
        if !isfile(dstpath) || force
            Downloads.download(joinpath(ftp_address, filename), dstpath)
        end
    end
end

baremodule Hosts
import Base: @enum
@enum Host::UInt8 human swine avian other
export Host
end # baremodule
using .Hosts

function Base.parse(::Type{Host}, s::Union{String, SubString})
    lcase = lowercase(s)
    lcase == "human" && return Hosts.human
    lcase == "swine" && return Hosts.swine
    lcase == "avian" && return Hosts.avian
    return Hosts.other
end

"""
    IncompleteSegmentData

Placeholder struct used to incrementally build `SegmentData`. You can construct
`SegmentData` using `SegmentData(::IncompleteSegmentData)`, after which the created
struct takes ownership of all fields of the incomplete data.
"""
mutable struct IncompleteSegmentData
    id::String
    host::Host
    segment::Segment
    serotype::Option{SeroType}
    year::Int16
    isolate::String
    proteins::Vector{ReferenceProtein}
    seq::Option{LongDNASeq}
end

"""
    SegmentData

Data structure representing an influenza segment from NCBI, and relevant information
that can be extracted from the files provided by NCBI.

Some fields in some records of input files may be missing. If so, the record
will be skipped, and not parsed to a `SegmentData`.

Notes:
* The SeroType is `none(SeroType)` if the segment is not Influenza A.
* The isolate is the isolate name, e.g. "A/Denmark/1/2021".
"""
struct SegmentData
    id::String
    host::Host
    segment::Segment
    serotype::Option{SeroType}
    year::Int16
    isolate::String
    proteins::Vector{ReferenceProtein}
    seq::LongDNASeq
end

# What we need to check here is the presence of proteins and sequence, since these
# are added separately to the incomplete records. The rest of the fields are checked
# when instantiating the IncompleteSegmentData. The biological plausibility of the
# data is not addressed at instantiation, only whether it's parsable.
function SegmentData(x::IncompleteSegmentData)
    isempty(x.proteins) && error("Cannot construct SegmentData from empty proteins")
    SegmentData(x.id, x.host, x.segment, x.serotype, x.year, x.isolate, x.proteins, unwrap(x.seq))
end

function Influenza.Reference(x::SegmentData)
    Influenza.Reference(x.id, x.segment, copy(x.seq), copy(x.proteins))
end

include("parse_genomeset.jl")
include("parse_fasta.jl")
include("parse_orfs.jl")
include("filter.jl")
include("clustering.jl")

"""
    parse_ncbi_records(genomeset, fasta, influenza_aa_dat, influenza_dat)

Parse NCBI records from the paths of four files: The cleaned "genomeset.dat"
output of `clean_genomeset`, "influenza.fna.gz", "influenza_aa.dat.gz" and
"influenza.dat.gz".
Return a `Dict{String, SegmentData}` with keys being the accession number.

Note that minimal filtering is done on the records. You may want to filter for
`isok_all(::Record)` to keep only those that pass all checks.
"""
function parse_ncbi_records(
    genomeset::String,
    fasta::String,
    influenza_aa_dat::String,
    influenza_dat::String
)
    dict = parse_cleaned_genomeset(genomeset)
    add_sequences!(dict, fasta)
    add_orfs!(dict, influenza_aa_dat, influenza_dat)

    # Now filter away any segments without seqs or orfs
    filter!(dict) do (k, v)
        !is_error(v.seq) && !isempty(v.proteins)
    end

    # This is a fairly slow operation, but necessary
    strip_false_termini!(dict)

    return Dict(k => SegmentData(v) for (k, v) in dict)
end

"""
    run_all(dir::AbstractString, split_dedup=true)

Convenience function: Downloads influenza data to `dir` if not already present,
then cleans the genomeset, then parses the records, filters them, deduplicate them,
and serialize both all records and the deduplicated ones, as well as writes FASTA
files to the directory.

If `split_dedup`, split the deduplicated files to one file per segment.

The executable `cd_hit_est` must be in the Julia PATH.
"""
function run_all(dir::AbstractString, split_dedup::Bool=true)
    Reference = Influenza.Reference
    println("Downloading...")
    download_influenza_data(dir)
    
    println("Cleaning genomeset.dat.gz...")
    cleaned_path = "$dir/genomeset.clean.dat.gz"
    if !isfile(cleaned_path)
        open(GzipDecompressorStream, "$dir/genomeset.dat.gz") do inio
            open(GzipCompressorStream, cleaned_path, "w") do outio
                clean_genomeset(outio, inio)
            end
        end
    end

    println("Parsing records...")
    data = parse_ncbi_records(
        "$dir/genomeset.clean.dat",
        "$dir/influenza.fna.gz",
        "$dir/influenza_aa.dat.gz",
        "$dir/influenza.dat.gz",
    )
    
    println("Filtering records...")
    filter!(data) do (k, v)
        isok_all(v)
    end

    println("Deduplicating...")
    (dedup, _) = cd_hit_deduplicate_all(data)
    
    println("Serializing")
    files = Dict("all" => Reference.(values(data)))
    if split_dedup
        for v in values(dedup)
            push!(get!(Vector{Reference}, files, string(v.segment)), Reference(v))
        end
    else
        files["dedup"] = Reference.(values(dedup))
    end
    for (basename, refs) in files
        store_references("$dir/$basename", refs, true)
    end
    println("Done!")
end

end # module
