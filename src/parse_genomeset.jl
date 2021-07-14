# Functionality to download, clean and parse the text
# files from NCBI into a Julia native datastructure

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

struct ProteinORF
    variant::Protein
    orfs::Vector{UnitRange{UInt16}}
end

mutable struct SegmentData
    id::String
    host::Host
    segment::Segment
    serotype::Option{SeroType}
    clade::Option{String}
    year::Int16
    name::String
    proteins::Vector{ProteinORF}
    seq::Option{LongDNASeq}
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

"""
    clean_genomeset(outpath::AbstractString, inpath::AbstractString)

Reads in `inpath`, a "genomeset.dat.gz" file, and write a cleaned copy to `outpath`.
The cleaning process renames mistyped or invalid data. 

Because it is implemented by manually correcting each observed mistake,
this function may fail, or incompletely clean the data in future versions of 
the NCBI Influenza data. In that case, you can assume that parsing incompletely
cleaned data will fail, and not parse invalid data.
"""
function clean_genomeset(outpath::AbstractString, inpath::AbstractString)
    file = open(inpath)
    outfile = open(outpath, "w")
    decompressed = GzipDecompressorStream(file)

    for line in eachline(decompressed) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        
        # Check correct number of fields
        if length(fields) != 11
            println("Error: Should have 11 fields \"$line\"")
            error()
        end

        gi, host, segment, subtype, country, year, len, name, age, gender, group = fields

        # Filter subtype
        subtype_upper = uppercase(subtype)
        subtype = if startswith(subtype_upper, "MIXED")
            ""
        elseif subtype in ["H1", "H11N9/N2"]
            ""
        elseif subtype in ["H1N2v", "H3N2v"]
            subtype[1:end-1]
        elseif subtype == "H3N6,H3"
            "H3N6"
        elseif subtype == "H6N1,H6"
            "H6N1"
        else
            subtype_upper
        end

        # Filter name
        nn = parse_genomeset_name(name)
        name = if is_error(nn)
            ""
        else
            unwrap(nn)
        end

        # Filter year
        year = if year in ["NON", "Unknown", "unknown"]
            ""
        else
            year
        end 

        println(outfile, join([gi, host, segment, subtype, country, year, len, name, age, gender, group], '\t'))
    end

    close(decompressed)
    close(outfile)
end

# Parses a name in genomeset.dat during cleaning, returning none if it's malformed.
# If it's unexpectedly malformed, throw an error"
function parse_genomeset_name(s::Union{String, SubString})::Option{String}
    isempty(s) && return none

    occursin(r"^Influenza [AB] [vV]irus", s) || error(s)
    ncodeunits(s) == 17 && return none
    rest = SubString(s, 18 + (codeunit(s, 18) == UInt(' ')):ncodeunits(s))

    isempty(rest) && return none
    rest2 = if codeunit(rest, 1) == UInt8('(') && codeunit(rest, ncodeunits(rest)) == UInt8(')')
        SubString(rest, 2:ncodeunits(rest)-1)
    else
        rest
    end

    m = match(r"\([Hh]\d+[nN]\d+\)$", rest2)
    rest3 = if m === nothing
        rest2
    else
        SubString(rest2, 1:ncodeunits(rest2) - ncodeunits(m.match))
    end
    
    return some(String(rest3))
end

###
function try_parse_from_integer(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = tryparse(UInt8, s)
    y === nothing && return none
    (iszero(y) | (y > 0x08)) && return none
    some(reinterpret(Segment, y - 0x01))
end

function try_parse(::Type{SeroType}, s::AbstractString)::Option{SeroType}
    st = tryparse(SeroType, s)
    st === nothing ? none : some(st)
end

function try_parse(::Type{SegmentData}, line::Union{String, SubString{String}})::Option{SegmentData}
    fields = split(line, '\t')
    serotype = if isempty(strip(fields[4]))
        none(SeroType)
    else
        some(@? try_parse(SeroType, fields[4]))
    end
    year = @? try_parse_year(fields[6])
    host = parse(Host, fields[2])
    gi = String(fields[1])
    name = String(fields[8])
    segment = @? try_parse_from_integer(Segment, fields[3])
    some(SegmentData(gi, host, segment, serotype, none(String), year, name, ProteinORF[], none(LongDNASeq)))
end

function try_parse_year(s::Union{String, SubString})::Option{Int16}
    isempty(s) && return none
    pos_slash_found = findfirst(isequal('/'), s)
    last_byte = pos_slash_found === nothing ? ncodeunits(s) : pos_slash_found - 1
    return some(parse(Int16, SubString(s, 1:last_byte)))
end

function parse_cleaned_genomeset(io::IO)::Dict{String, SegmentData}
    result = Dict{String, SegmentData}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        data = @unwrap_or try_parse(SegmentData, line) continue
        gi = data.id
        @assert !haskey(result, gi) "GB identifier $(gi) not unique"
        result[gi] = data
    end
    result
end

"""
    parse_cleaned_genomeset(::IO) -> Dict{String, SegmentData}

Parse the file generated by `clean_genomeset`, and return a `Dict` where the
keys are the GI accession number, and the value is a `SegmentData` object.

This function attempts to validate the input data, but SKIPS any row with critical
fields missing or invalid. Erroring on the first invalid record would make the
dataset cleaning too big a task. Therefore, when updating the dataset, you should check
how many individual records are returned, or if an alarming number are skipped.
"""
parse_cleaned_genomeset(s::AbstractString) = open(parse_cleaned_genomeset, s)