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

using Base: RefValue

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

include("parse_genomeset.jl")
include("parse_fasta.jl")
include("parse_orfs.jl")
include("filter.jl")

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