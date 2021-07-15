#=
TODO:
* It should not be possible to have a SegmentData without a sequence, since
that is probably not useful. Instead, have a dummy type IncompleteSegmentData or
whatever, and then create SegmentData from all input files, discarding any incomplete
ones. Biologically implausible ones may be instantiated, and can be filtered afterwards.
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
    SegmentData

Data structure representing an influenza segment from NCBI, and relevant information
that can be extracted from the files provided by NCBI.

Some fields in some records of input files may be missing. If so, the record
will be skipped, and not parsed to a `SegmentData`.

Notes:
* The SeroType is `none(SeroType)` if the segment is not Influenza A.
* The isolate is the isolate name, e.g. "A/Denmark/1/2021".
"""
mutable struct SegmentData
    id::String
    host::Host
    segment::Segment
    serotype::Option{SeroType}
    year::Int16
    isolate::String
    proteins::Vector{ProteinORF}
    seq::Option{LongDNASeq}
end

include("parse_genomeset.jl")
include("parse_fasta.jl")
include("parse_orfs.jl")
include("filter.jl")

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
        !is_error(v.seq) && !isempty(v.protieins)
    end
end

export clean_genomeset,
    parse_cleaned_genomeset,
    isok_minimum_proteins,
    isok_orf_length,
    isok_seq_length,
    isok_translatable,
    isok_all,
    SegmentData,
    add_sequences!,
    add_orfs!

#=
sd = open(NCBIInfluenzaData.parse_cleaned_genomeset, "/Users/jakobnissen/Documents/ssi/projects/flupipe/build_references/processed/genomeset.filt.dat")
NCBIInfluenzaData.add_sequences!(sd, "/Users/jakobnissen/Documents/ssi/projects/flupipe/build_references/download/influenza.fna.gz")
NCBIInfluenzaData.add_orfs!(
    sd,
    "/Users/jakobnissen/Documents/ssi/projects/flupipe/build_references/download/influenza_aa.dat.gz",
    "/Users/jakobnissen/Documents/ssi/projects/flupipe/build_references/download/influenza.dat.gz"
)
=#

end # module