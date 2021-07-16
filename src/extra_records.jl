# This file contains code to add extra records using annotation.tbl and
# a FASTA file.
"""
    NCBIInfluenzaData.Extra

This module contains functionality to create new `SegmentData` objects from existing
sequences not from NCBI, such that influenza sequences from other sources can be
parsed to the same data structure as those from NCBI.

It requires two files: A FASTA file where the header is of the following format:
"IDENTIFIER|HOST|SEGMENT|SEROTYPE|YEAR|ISOLATE".
* Identifier can be any string, not including `|` or trailing/leading whitespace
* Host can be "Human", "Swine", "Avian", or "Other", or in lowercase.
* Segment must be parsable into InfluenzaCore.Segment
* Serotype can be the empty string or parsable into InfluenzaCore.SeroType
* Year must be 4 characters in 0-9
* Isolate must be the standardized e.g. "A/Denmark/1/2021"

The second file is an annotation.tbl file produced by The NCBI Influenza Virus
Sequence Annotation Tool. The headers in the annotation.tbl file must match the
headers in the FASTA file.
"""


function try_parse_fasta_header(
    ::Type{IncompleteSegmentData},
    fields::Vector{SubString{String}},
    line::Union{String, SubString{String}}
)::Option{IncompleteSegmentData}
    # Check there are 6 fields
    (@? try_split!(fields, line, UInt8('|'))) == 6 || return none
    id, s_host, s_segment, s_serotype, s_year, isolate = fields
    # Check ID is not empty
    isempty(id) && return none
    # Check host
    host = parse(Host, s_host)
    host in (Hosts.human, Hosts.avian, Hosts.swine) || return none
    # Check segment
    segment::Segment = let
        s = tryparse(Segment, s_segment)
        s === nothing && return none
        s
    end
    # Check serotype
    serotype::Option{SeroType} = let
        if isempty(s_serotype)
            none(SeroType)
        else
            s = tryparse(SeroType, s_serotype)
            s === nothing && return none
            some(s)
        end
    end
    # Check year
    year = let
        s = tryparse(Int16, s_year)
        s === nothing && return none
        s
    end
    # Check isolate format
    validate_isolate_format(isolate) || return none
    return some(IncompleteSegmentData(
        id, host, segment, serotype, year, isolate, ProteinORF[], none(LongDNASeq)
    ))
end

function try_parse_fasta_header(::Type{IncompleteSegmentData}, line::Union{String, SubString{String}})
    try_parse_fasta_header(IncompleteSegmentData, Vector{SubString{String}}(undef, 6), line)
end

#=
struct Unbounded end
const FeatureRange = Tuple{Union{Unbounded, UInt}, Union{Unbounded, UInt}}
struct Feature
    key::String
    ranges::Vector{FeatureRange}
    qualifiers::Dict{String, String}
end

# This format is absolutely horrendous. The parser here is not a generic parser for
# genbank formats.
function try_parse_genbank_feature_table(lines::Vector{SubString{String}})::Option{Vector{Feature}}
    # First split it into individual features
    # The features consists of N "headers" where the first field is not a tab
    # followed by M lines where the field IS a tab.
    fieldbuffer = Vector{SubString{String}}(undef, 5)
    result = Vector{Feature}()
    feature_groups = Vector{Vector{SubString{String}}}()
    filter!(!isempty, lines)
    start = 1
    last_started_tab = false
    for (i, line) in enumerate(lines)
        this_starts_tab = codeunit(line, 1) === UInt8('\t')
        if last_started_tab && !this_starts_tab
            push!(feature_groups, lines[start:i-1])
            start = i
        end
        last_started_tab = this_starts_tab
    end
    push!(feature_groups, lines[start:end])

    # Then parse each of the features
    for feature in feature_groups
        push!(result, @? try_parse_genbank_feature(fieldbuffer, feature))
    end

    return some(result)
end

# Example:
# <3522    3572    CDS
# 3706    >4197
#                         product        Yip2p
#                         prot_desc      similar to human polyposis locus protein 1 (YPD)
function try_parse_genbank_feature(
    v::Vector{SubString{String}},
    lines::Vector{SubString{String}}
)::Option{Feature}
    length(lines) < 2 && return none

    ranges = FeatureRange[]
    is_first_line = true
    line_num = 1
    key = nothing
    n_fields = (@? try_split!(v, lines[line_num], UInt8('\t')))

    # While the first field is not empty, we are still parsing the feature range
    while !isempty(@inbounds v[1])
        line_num += 1
        line_num > length(lines) && return none
        if is_first_line
            n_fields < 3 && return none
            key = @inbounds v[3]
            isempty(key) && return none
            is_first_line = false
        else
            n_fields > 2 && !isempty(v[3]) && return none 
        end
        fst = @? try_parse_feature_number(@inbounds(v[1]), true)
        lst = @? try_parse_feature_number(@inbounds(v[2]), false)
        push!(ranges, (fst, lst))
        n_fields = (@? try_split!(v, lines[line_num], UInt8('\t')))
    end

    # Parse the qualifier lines
    qualifiers = Dict{String, String}()
    while line_num ≤ length(lines)
        n_fields = (@? try_split!(v, lines[line_num], UInt8('\t')))
        n_fields == 5 || return none
        all(isempty(v[i]) for i in 1:3) || return none
        qualkey, qualval = (@inbounds v[4]), (@inbounds v[5])
        (isempty(qualkey) || isempty(qualval)) && return none
        qualifiers[qualkey] = qualval
        line_num += 1
    end

    return some(Feature(key, ranges, qualifiers))
end

function try_parse_feature_number(
    s::Union{SubString{String}, String},
    isfirst::Bool
)::Option{Union{Unbounded, UInt}}
    isempty(s) && return none
    firstbyte = codeunit(s, 1)
    if (isfirst  & (firstbyte === UInt8('<'))) ||
       (!isfirst & (firstbyte === UInt8('>')))
        return some(Unbounded())
    end
    n = tryparse(UInt, s)
    n === nothing && return none
    return some(n::UInt)
end
=#

#=
function fill_orfs_extra_records!(records::Dict{String, IncompleteSegmentData}, annotpath::String)
    chunks = maybe_gzip(annotpath) do file
        partition_chunks(file)
    end
    filled_ids = Set{String}()
    for chunk in chunks
        parse_annot_segment_data!(records, filled_ids, chunk)
    end

    # If any were missing, error here
    if length(filled_ids) < length(records)
        missing_ids = setdiff(keys(records), filled_ids)
        error("Missing header from annot.tbl: \"$(first(missing_ids))\"")
    end
    return records
end

# Please have a look at the files in annotation.tbl to understand this
# it partitions each individual record which is separated by equal signs.
function partition_chunks(io::IO)::Vector{Vector{String}}
    lines = collect(eachline(io))
    chunks = Vector{Vector{String}}()
    chunk = String[]
    for line in lines
        if startswith(line, "=========")
            push!(chunks, copy(chunk))
            empty!(chunk)
        else
            push!(chunk, line)
        end
    end
    push!(chunks, chunk)
end

function parse_annot_segment_data!(
    records::Dict{String, SegmentData},
    filled_ids::Set{String},
    chunk::Vector{String}
)
    @assert startswith(first(chunk), ">Feature ")
    header = first(chunk)[10:end]
    id = first(split(header, '|'))
    @assert haskey(records, id)
    push!(filled_ids, id)
    data = records[id]
    orfs = UnitRange{UInt16}[]
    reading_cds = false
    for line in @view chunk[2:end]
        isempty(strip(line)) && continue
        fields = split(line)
        if isnumeric(line[1]) && length(fields) > 2 && fields[3] == "CDS"
            reading_cds = true
        end
        if isnumeric(line[1]) && reading_cds
            push!(orfs, parse(UInt16, fields[1]) : parse(UInt16, fields[2]))
        end
        if reading_cds && fields[1] == "gene"
            variant = expect(parse(Protein, fields[2]), "Error on chunk $header")
            reading_cds = false
            push!(data.proteins, ProteinORF(variant, copy(orfs)))
            empty!(orfs)
        end
    end
    records
end
=#

#=
function add_extra_records!(
    segment_data::Dict{String, SegmentData},
    annotpath::String,
    fastapath::String,
    deduplicate::Bool
)
    # id => data from the FASTA file
    records = get_extra_records(fastapath, deduplicate)

    # Add info from annotated ORFS
    fill_orfs_extra_records!(records, annotpath)

    # Merge with the original segment_data dict
    @assert isdisjoint(keys(segment_data), keys(records))
    merge!(segment_data, records)
end

function get_extra_records(fastapath::String, deduplicate::Bool)
    records = Dict{String, SegmentData}()
    maybe_gzip(fastapath) do io
        reader = FASTA.Reader(io)
        record = FASTA.Record()
        while !eof(reader)
            read!(reader, record)
            data = parse_data_from_header(FASTA.header(record)::String, deduplicate)
            data.seq[] = some(FASTA.sequence(LongDNASeq, record))
            @assert !haskey(records, data.id) "Duplicate key $(data.id)"
            records[data.id] = data
        end
    end
    records
end

function parse_data_from_header(header::String, deduplicate::Bool)::SegmentData
    # Looks like this "EPI1843298|avian|NS|H1N1|2021|A/Denmark/1/2021|CLADE"
    # some of the fields may be empty
    fields = split(header, '|')
    @assert length(fields) == 7 header
    id, s_species, s_segment, s_subtype, s_year, name, clade = fields
    @assert !isempty(id)
    @assert !isempty(name)
    species = parse(Species, s_species)
    if species == other
        println(s_species)
    end
    @assert species != other
    return SegmentData(
        id,
        species,
        unwrap(parse(Segment, s_segment)),
        unwrap(parse(SubType, s_subtype)),
        deduplicate,
        isempty(clade) ? none(String) : some(String(clade)),
        parse(Int16, s_year),
        name,
        ProteinORF[],
        RefValue(none(LongDNASeq))
    )
end


=#