# Functionality to download, clean and parse the text
# files from NCBI into a Julia native datastructure

"""
    clean_genomeset(outio::IO, inio::IO)

Reads in data from `inio`, representing the decompressed data of "genomeset.dat.gz"
and write a cleaned copy to `outio`. The cleaning process renames mistyped or invalid data. 

Because it is implemented by manually correcting each observed mistake,
this function may fail, or incompletely clean the data in future versions of 
the NCBI Influenza data. In that case, you can assume that parsing incompletely
cleaned data will fail, and not parse invalid data.
"""
function clean_genomeset(outio::IO, inio::IO)
    fields = Vector{SubString{String}}(undef, 11)

    for line in eachline(inio) |> Map(strip) ⨟ Filter(!isempty)
        # Check the correct number of fields
        if unwrap_or(try_split!(fields, line, '\t'), 0) < 11
            error("Error: Should have 11 fields \"$line\"")
        end

        gi, host, segment, subtype, country, year, len, isolate, age, gender, group = fields

        # Filter isolate
        maybe_alphaname = parse_genomeset_isolate(isolate)
        (isalpha, isolate) = if is_error(maybe_alphaname)
            nothing, ""
        else
            unwrap(maybe_alphaname)
        end

        # Filter subtype - set it to "" if not an influenza A
        subtype_upper = uppercase(subtype)
        subtype = if isalpha === nothing
            ""
        elseif !isalpha
            ""
        elseif startswith(subtype_upper, "MIXED")
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

        # Filter year
        year = if year in ["NON", "Unknown", "unknown"]
            ""
        else
            year
        end 

        println(outio, join([gi, host, segment, subtype, country, year, len, isolate, age, gender, group], '\t'))
    end
end

# Parses a isolate in genomeset.dat during cleaning, returning none if it's malformed.
# If it's unexpectedly malformed, throw an error"
# Returns Option(is_alphavirus, name)
function parse_genomeset_isolate(s::Union{String, SubString})::Option{Tuple{Bool, String}}
    # They look like this: "Influenza A virus (A/swine/Minnesota/66960/2006(H3N2))"
    isempty(s) && return none

    # Strip "Influenza A Virus" beginning off
    occursin(r"^Influenza [AB] [vV]irus", s) || error(s)
    ncodeunits(s) == 17 && return none
    rest = SubString(s, 18 + (codeunit(s, 18) == UInt(' ')):ncodeunits(s))

    # Strip parenthesis off
    isempty(rest) && return none
    rest2 = if codeunit(rest, 1) == UInt8('(') && codeunit(rest, ncodeunits(rest)) == UInt8(')')
        SubString(rest, 2:ncodeunits(rest)-1)
    else
        rest
    end

    # Some of them ends with "(mixed)" - these are useless and not real isolates.
    if endswith(rest2, "(mixed)")
        return none
    end

    # Now it should look like e.g. "A/mallard/Interior Alaska/8BM3586/2008(H3N8)"
    # strip off the serotype designation, if present.
    m = match(r"\([Hh]\d+[nN]\d+\)$", rest2)
    rest3 = if m === nothing
        rest2
    else
        SubString(rest2, 1:ncodeunits(rest2) - ncodeunits(m.match))
    end

    # Now it should look like "A/northern shoveler/Mississippi/09OS168/2009"
    # verify it actually does so!
    validate_isolate_format(rest3) || return none
    isalpha = first(rest3) == 'A'
    return some((isalpha, String(rest3)))
end

function validate_isolate_format(s::AbstractString)
    m = match(r"^([AB])/[^/]+/[^/]+(?:/[^/]+)?/\d{4}$", s)
    return m !== nothing
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

function try_parse(
    ::Type{IncompleteSegmentData},
    fields::Vector{SubString{String}},
    line::Union{String, SubString{String}}
)::Option{IncompleteSegmentData}
    nfields = try_split!(fields, line, '\t')
    unwrap_or(nfields, 0) == 11 || error("Expected 11 fields in \"$line\"")
    serotype = if isempty(strip(fields[4]))
        none(SeroType)
    else
        some(@? try_parse(SeroType, fields[4]))
    end
    year = @? try_parse_year(fields[6])
    host = parse(Host, fields[2])
    gi = String(fields[1])
    isolate = String(fields[8])
    segment = @? try_parse_from_integer(Segment, fields[3])
    some(IncompleteSegmentData(gi, host, segment, serotype, year, isolate, ReferenceProtein[], none(LongDNASeq)))
end

function try_parse(
    ::Type{IncompleteSegmentData},
    line::Union{String, SubString{String}}
)::Option{IncompleteSegmentData}
    try_parse(IncompleteSegmentData, Vector{SubString{String}}(undef, 11), line)
end

function try_parse_year(s::Union{String, SubString})::Option{Int16}
    isempty(s) && return none
    pos_slash_found = findfirst(isequal('/'), s)
    last_byte = pos_slash_found === nothing ? ncodeunits(s) : pos_slash_found - 1
    return some(parse(Int16, SubString(s, 1:last_byte)))
end

function parse_cleaned_genomeset(io::IO)::Dict{String, IncompleteSegmentData}
    result = Dict{String, IncompleteSegmentData}()
    buffer = Vector{SubString{String}}(undef, 11)
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        data = @unwrap_or try_parse(IncompleteSegmentData, buffer, line) continue
        gi = data.id
        @assert !haskey(result, gi) "GB identifier $(gi) not unique"
        result[gi] = data
    end
    result
end

"""
    parse_cleaned_genomeset(::IO) -> Dict{String, IncompleteSegmentData}

Parse the file generated by `clean_genomeset`, and return a `Dict` where the
keys are the GI accession number, and the value is a `IncompleteSegmentData` object.

This function attempts to validate the input data, but SKIPS any row with critical
fields missing or invalid. Erroring on the first invalid record would make the
dataset cleaning too big a task. Therefore, when updating the dataset, you should check
how many individual records are returned, or if an alarming number are skipped.
"""
parse_cleaned_genomeset(s::AbstractString) = open(parse_cleaned_genomeset, s)