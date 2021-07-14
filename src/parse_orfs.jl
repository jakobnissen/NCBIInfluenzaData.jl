function add_orfs!(
    segment_data::Dict{String, SegmentData},
    inf_aa_dat_path::String,
    inf_dat_path::String,
)
    accession_protein_map = maybe_gzip(parse_inf_aa, inf_aa_dat_path)
    maybe_gzip(inf_dat_path) do io
        add_orfs!(segment_data, accession_protein_map, io)
    end
end

function try_parse(::Type{Protein}, s::AbstractString)::Option{Protein}
    p = tryparse(Protein, s)
    p === nothing ? none : some(p)
end

# Parse the inf_aa file and get a map from accession => protein
function parse_inf_aa(io::IO)::Dict{String, Protein}
    result = Dict{String, Protein}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        accession = first(fields)
        variant = @unwrap_or try_parse(Protein, fields[3]) continue
        @assert !haskey(result, accession) "Duplicate key $accession"
        result[accession] = variant
    end
    result
end

function add_orfs!(segment_data::Dict{String, SegmentData}, accession_protein_map::Dict{String, Protein}, io::IO)
    n_updates = 0
    proteinbuffer = ProteinORF[]
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        gb_accession = first(fields)
        haskey(segment_data, gb_accession) || continue
        @assert isodd(length(fields))
        data = segment_data[gb_accession]
        n_updates += 1
        is_bad = false
        empty!(proteinbuffer)
        for (protein_accession, orf_field) in zip(@view(fields[2:2:end]), @view(fields[3:2:end]))

            # For some reason, some of these protein accessions are not actually present
            # in the influenza_aa.dat file.
            if !haskey(accession_protein_map, protein_accession)
                is_bad = true
                break
            end
            variant = accession_protein_map[protein_accession]
            orfs = parse_orf_field(orf_field)
            if is_error(orfs)
                is_bad = true
                break
            end
            push!(proteinbuffer, ProteinORF(variant, unwrap(orfs)))
        end
        if !is_bad
            append!(data.proteins, proteinbuffer)
        end
    end
    return n_updates
end
        
function parse_orf_field(s::Union{String, SubString{String}})::Option{Vector{UnitRange{UInt16}}}
    # If it looks like gb|AB266090:<411->632, the ORF is not present
    # in the reference, and we skip it
    if occursin('>', s) || occursin('<', s)
        return none
    end

    # Multiple ORFs in a gene makes them enclosed in brackets
    # e.g. (gb|AB212651:26-52, 741-1007)
    s1 = strip(s, ['(', ')'])
    p = findfirst(isequal(':'), s1)
    p === nothing && return none
    s2 = @view s1[p+1:ncodeunits(s1)]
    orf_strings = split(s2, ", ")

    orfstrings = split(s2, ", ")
    result = map(parse_range, orfstrings)

    # Return none if it's not divisible by three and thus cant be ORF
    iszero(sum(length, result) % 3) || return none
    some(result)
end

function parse_range(s::AbstractString)::UnitRange{UInt16}
    p2 = findfirst(isequal('-'), s)
    # If it's not a range, it must be a single number
    range = if p2 === nothing
        n = parse(UInt16, s)
        n:n
    else
        start = parse(UInt16, @view s[1:p2-1])
        stop = parse(UInt16, s[p2+1:ncodeunits(s)])
        start:stop
    end
    iszero(range) && error("ORF range \"$s\" cannot contain the position zero.")
end