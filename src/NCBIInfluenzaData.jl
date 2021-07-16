#=
To do:

* Add clustering + serialization
* Add extra records - maybe using Influenza.jl functionality?

=#

"""
    NCBIInfluenzaData

Package for parsing and filtering NCBI's bulk influenza data. This can be used to
extract and work with tens of thousands of sequences in bulk.
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

function download_influenza_data(dst_dir::AbstractString; force::Bool=false)
    isdir(dst_dir) && !force && return nothing
    isdir(dst_dir) || mkdir(dst_dir)
    ftp_address = "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/"
    for filename in [
        "genomeset.dat.gz",
        "influenza.dat.gz",
        "influenza.fna.gz",
        "influenza_aa.dat.gz"
        ]
        Downloads.download(joinpath(ftp_address, filename), joinpath(dst_dir, filename))
        println("Downloaded $filename")
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
    ProteinORF

A data structure that represents an Influenza A/B protein and the positions of the ORFs
that encode the protein for a given nucleotide sequence.
"""
struct ProteinORF
    variant::Protein
    orfs::Vector{UnitRange{UInt16}}
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
    proteins::Vector{ProteinORF}
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
    proteins::Vector{ProteinORF}
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

include("parse_genomeset.jl")
include("parse_fasta.jl")
include("parse_orfs.jl")
include("filter.jl")
include("extra_records.jl")

"""
    parse_ncbi_records(genomeset, fasta, influenza_aa_dat, influenza_dat)

Parse NCBI records from the paths of four files: The cleaned "genomeset.dat"
output of `clean_genomeset`, "influenza.fna.gz", "influenza_aa.dat.gz" and
"influenza.dat.gz".
Return a `Dict{String, SegmentData}` with keys being the accession number.
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

export clean_genomeset,
    parse_ncbi_records,
    isok_minimum_proteins,
    isok_orf_length,
    isok_seq_length,
    isok_translatable,
    isok_all,
    SegmentData

#=
br = "/Users/jakobnissen/Documents/ssi/projects/flupipe/build_references/"
NCBIInfluenzaData.parse_ncbi_records(
    br * "processed/genomeset.filt.dat",
    br * "download/influenza.fna.gz",
    br * "download/influenza_aa.dat.gz",
    br * "download/influenza.dat.gz"
)
=#

end # module