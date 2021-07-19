"""
    add_orfs!(::Dict{String, IncompleteSegmentData}, inf_aa_dat_path::String, inf_dat_path::String)

Parse the files `influenza_aa.dat` and `influenza.dat`, extract information
about ORFs, and update the `IncompleteSegmentData` `proteins` field for all segment data
present in the input dict. Returns the number of updated `IncompleteSegmentData` instances.
"""
function add_orfs!(
    segment_data::Dict{String, IncompleteSegmentData},
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

function add_orfs!(segment_data::Dict{String, IncompleteSegmentData}, accession_protein_map::Dict{String, Protein}, io::IO)
    n_updates = 0
    proteinbuffer = ReferenceProtein[]
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
            push!(proteinbuffer, ReferenceProtein(variant, unwrap(orfs)))
        end
        if !is_bad
            resize!(data.proteins, length(proteinbuffer))
            copy!(data.proteins, proteinbuffer)
        end
    end
    return n_updates
end
        
function parse_orf_field(s::Union{String, SubString{String}})::Option{Vector{UnitRange{UInt32}}}
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

function parse_range(s::AbstractString)::UnitRange{UInt32}
    p2 = findfirst(isequal('-'), s)
    # If it's not a range, it must be a single number
    range = if p2 === nothing
        n = parse(UInt32, s)
        n:n
    else
        start = parse(UInt32, @view s[1:p2-1])
        stop = parse(UInt32, s[p2+1:ncodeunits(s)])
        start:stop
    end
    iszero(first(range)) && error("ORF range \"$s\" cannot contain the position zero.")
    return range
end

"""
    strip_false_termini!(data::Dict{String, IncompleteSegmentData})

Some sequences given by NCBI are too long, containing linker or primer contamination.
True influenza sequences begin with AGCAAAAGCAGG and ends with CTTGTTTCTACT. These
termini are highly conserved, but may have a single substitution.
We can allow sequences that do not include these noncoding ends, but we cannot
allow artifacts flanking the true termini.

This function must be called on completely initialized `IncompleteSegmentData`
with both ORFs and sequence, since it modifies both the sequence and ORFs.
"""
function strip_false_termini!(data::Dict{String, IncompleteSegmentData})
    toremove = Set{String}() # we remove these keys
    n_stripped = 0
    for (id, segment_data) in data
        result = strip_false_termini!(segment_data)
        n_stripped += @unwrap_or result begin
            push!(toremove, id)
            continue
        end
    end
    if !iszero(length(toremove))
        filter!(data) do (k, v)
            !in(k, toremove)
        end
    end
    return (n_stripped, length(toremove))
end

# Returns some(true/false), whether data was stripped off, and none if
# the ORF goes into the stripped termini, and the sequence must be discarded
function strip_false_termini!(data::IncompleteSegmentData)::Option{Bool}
    is_error(data.seq) && return some(false)
    false_termini = false_termini_length(data)
    false_termini === nothing && return some(false)
    trim5, trim3 = false_termini
    iszero(trim5) && iszero(trim3) && return some(false)
    seqlen = UInt32(length(unwrap(data.seq)))

    # For safety, error if there are no ORFs. This function should not be called
    # on segments without ORFs, since the orfs need to be modified by this funcion
    isempty(data.proteins) && error("Termini must be stripped after adding ORFs to data")

    # Trim sequence and ORFs.
    data.seq = some(unwrap(data.seq)[trim5+1:seqlen-trim3])
    for (proti, protein) in enumerate(data.proteins), (orfi, orf) in enumerate(protein.orfs)
        if trim5 ≥ first(orf) || (seqlen - trim3) < last(orf)
            return none
        else
            data.proteins[proti].orfs[orfi] = (first(orf) - trim5):(last(orf) - trim5)
        end
    end
    return some(true)
end

# We allow to strip up to 100 bp in each end off
function false_termini_length(data::IncompleteSegmentData)::Union{Nothing, Tuple{UInt32, UInt32}}
    seq = @unwrap_or data.seq (return nothing)
    trim5, trim3 = UInt32(0), UInt32(0)
    # Find true beginning of sequence, if it's within first 100 bp
    p = approxsearch(seq, dna"AGCAAAAGCAGG", 1)
    if !isempty(p)
        if first(p) < 100
            trim5 = UInt32(first(p) - 1)
        end
    end
    # Find true end of sequence
    p = approxrsearch(seq, dna"CTTGTTTCTCCT", 1)
    if !isempty(p)
        trim3 = UInt32(lastindex(seq) - last(p))
        if trim3 > 100
            trim3 = UInt32(0)
        end
    end
    return (trim5, trim3)
end
