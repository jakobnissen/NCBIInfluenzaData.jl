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
    fields = Vector{SubString{String}}(undef, 11)
    for line in eachline(io) |> ifilter(ln -> !isempty(strip(ln)))
        nfields = try_split!(fields, line, '\t')
        if unwrap_or(nfields, 0) != 11
            error("Expected 11 fields, got", line)
        end
        accession = first(fields)
        variant = @unwrap_or try_parse(Protein, fields[3]) continue
        haskey(result, accession) && error("Duplicate key $accession")
        result[accession] = variant
    end
    result
end

function add_orfs!(segment_data::Dict{String, IncompleteSegmentData}, accession_protein_map::Dict{String, Protein}, io::IO)
    n_updates = 0
    proteinbuffer = ReferenceProtein[]
    fields = Vector{SubString{String}}(undef, 128)
    for line in eachline(io) |> imap(strip) |> ifilter(!isempty)
        nfields = unwrap(try_split!(fields, line, '\t'))
        isodd(nfields) || error("Must be an odd-numbered fields")
        gb_accession = first(fields)
        haskey(segment_data, gb_accession) || continue
        data = segment_data[gb_accession]
        n_updates += 1
        is_bad = false
        empty!(proteinbuffer)
        for (protein_accession, orf_field) in zip(@view(fields[2:2:nfields]), @view(fields[3:2:nfields]))

            # For some reason, some of these protein accessions are not actually present
            # in the influenza_aa.dat file.
            if !haskey(accession_protein_map, protein_accession)
                is_bad = true
                break
            end
            variant = accession_protein_map[protein_accession]
            m_orfs = parse_orf_field(orf_field)
            if is_error(data.seq) || is_error(m_orfs)
                is_bad = true
                break
            end

            orfs, seq = unwrap(m_orfs), unwrap(data.seq)

            # NCBI includes stop codon in the ORF. I don't. So if the last three bases
            # of the ORF is a stop, remove them from the ORF.
            # We assume last ORF is 3 bp in length
            if (
                !isempty(orfs) &&
                length(last(orfs)) > 3
            )
                range = last(orfs)[end-2:end]
                if (
                    checkbounds(Bool, eachindex(seq), range) &&
                    all(!BioSequences.isambiguous, seq[range]) &&
                    Influenza.is_stop(BioSequences.DNACodon(seq[range]))
                )
                    orfs[end] = orfs[end][1:end-3]
                end
            end
            push!(proteinbuffer, ReferenceProtein(variant, orfs))
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

    orfstrings = split(s2, ", ")
    result = map(parse_range, orfstrings)

    # Return none if it's not divisible by three and thus cant be ORF
    iszero(sum(length, result, init=0) % 3) || return none
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
    strip_false_termini!(data::Dict{String, IncompleteSegmentData}, filter_termini=true)

Some sequences given by NCBI are too long, containing linker or primer contamination.
True influenza sequences begin with AGCAAAAGCAGG and ends with CTTGTTTCTACT. These
termini are highly conserved, but may have a single substitution.
We can allow sequences that do not include these noncoding ends, but we cannot
allow artifacts flanking the true termini.
If `filter_termini`, remove sequences where the termini cannot be found within the
first/last few basepairs using approximate sequence search.

This function must be called on completely initialized `IncompleteSegmentData`
with both ORFs and sequence, since it modifies both the sequence and ORFs.
"""
function strip_false_termini!(
    data::Dict{String, IncompleteSegmentData},
    filter_termini::Bool
)
    fwquery = BioSequences.ApproximateSearchQuery(dna"AGCAAAAGCAGG", :forward)
    rvquery = BioSequences.ApproximateSearchQuery(dna"CTTGTTTCTCCT", :backward)
    good = Set{String}()

    for (id, segment_data) in data
        seq = @unwrap_or segment_data.seq begin
            continue
        end
        fiveprime, threeprime = false_termini_length(seq, fwquery, rvquery)

        # We discard something like 93% of all seqs here
        if fiveprime === nothing || threeprime === nothing
            if filter_termini
                continue
            else
                fiveprime = fiveprime === nothing ? UInt32(0) : fiveprime
                threeprime = threeprime === nothing ? UInt32(0) : threeprime
            end
        end
        isgood = strip_false_termini!(segment_data, fiveprime, threeprime)
        isgood && push!(good, id)
    end

    filter!(data) do (k, _)
        k in good
    end
    return data
end

# Returns whether 
function strip_false_termini!(
    data::IncompleteSegmentData,
    fiveprime::UInt32,
    threeprime::UInt32,
)::Bool
    iszero(fiveprime) && iszero(threeprime) && return true
    seqlen = UInt32(length(unwrap(data.seq)))

    # For safety, error if there are no ORFs. This function should not be called
    # on segments without ORFs, since the orfs need to be modified by this funcion
    isempty(data.proteins) && error("Termini must be stripped after adding ORFs to data")

    # Trim sequence and ORFs.
    data.seq = some(unwrap(data.seq)[fiveprime+1:seqlen-threeprime])
    for (proti, protein) in enumerate(data.proteins), (orfi, orf) in enumerate(protein.orfs)
        if threeprime ≥ first(orf) || (seqlen - threeprime) < last(orf)
            return false
        else
            data.proteins[proti].orfs[orfi] = (first(orf) - fiveprime):(last(orf) - fiveprime)
        end
    end
    return true
end

# We allow to strip up to 100 bp in each end off
"Returns (5', 3'), where each is nothing, if no termini was found,
or an UInt32 if one was found."
function false_termini_length(
    seq::LongDNASeq,
    fwquery::BioSequences.ApproximateSearchQuery,
    bwquery::BioSequences.ApproximateSearchQuery,
)::NTuple{2, Union{Nothing, UInt32}}
    trim5, trim3 = UInt32(0), UInt32(0)
    SEARCHLEN = 100

    # Find true beginning of sequence, if it's within first SEARCHLEN bp
    p = BioSequences.approxsearch(seq, fwquery, 2, 1, SEARCHLEN)
    trim5 = isempty(p) ? nothing : UInt32(first(p) - 1)

    # Find true end of sequence
    p = BioSequences.approxrsearch(seq, bwquery, 2, lastindex(seq), lastindex(seq)-SEARCHLEN)
    trim3 = isempty(p) ? nothing : UInt32(lastindex(seq) - last(p))
    return (trim5, trim3)
end
