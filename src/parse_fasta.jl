function add_sequences!(segment_data::Dict{String, IncompleteSegmentData}, io::IO)
    n_updated = 0
    record = FASTA.Record()
    reader = FASTA.Reader(io)
    seen = Set{String}()
    while !eof(reader)
        read!(reader, record)
        # headers begins with e.g. gi|59292|gb|X53029|Influenza A virus
        gb_accession = split(FASTA.identifier(record)::String, '|')[4]
        data = get(segment_data, gb_accession, nothing)
        data === nothing && continue
        if in(seen, gb_accession)
            error("Duplicate GI accession: \"$gb_accession\"")
        else
            push!(seen, gb_accession)
        end
        data.seq = some(FASTA.sequence(LongDNASeq, record))
        n_updated += 1
    end
    return n_updated
end

# Some reference sequences are too long. True influenza sequences ends with CTTGTTTCTACT
# and begins with AGCAAAAGCAGG. (with perhaps a single substitution.) We can
# allow sequences that do not include these ends, but we cannot allow sequences
# that include artifacts flanking these true termini.
# Example: KJ476677 contains 19 bp 3' of true 3' terminal.
function strip_false_termini!(data::Dict{String, IncompleteSegmentData})
    toremove = Set{String}() # we remove these keys
    n_stripped = 0
    for (id, segment_data) in data
        false_termini =  false_termini_length(segment_data)
        false_termini === nothing && continue
        trim5, trim3 = false_termini
        iszero(trim5) && iszero(trim3) && continue
        seqlen = UInt16(length(unwrap(segment_data.seq[])))
        segment_data.seq[] = some(unwrap(segment_data.seq[])[trim5+1:seqlen-trim3])
        for (proti, protein) in enumerate(segment_data.proteins), (orfi, orf) in enumerate(protein.orfs)
            if trim5 ≥ first(orf) || (seqlen - trim3) < last(orf)
                push!(toremove, id)
                break
            else
                segment_data.proteins[proti].orfs[orfi] = (first(orf) - trim5):(last(orf) - trim5)
                n_stripped += 1
            end
        end
    end
    if !iszero(length(toremove))
        filter!(data) do (k, v)
            !in(k, toremove)
        end
        println("Removed $(length(toremove)) segments due to ORFs in false termini")
    end
    println("Stripped false termini off $(n_stripped) segment data")
end

# Returns whether data was stripped off
function strip_false_termini!(data::IncompleteSegmentData)::Bool
    is_error(data.seq) && return false
    false_termini =  false_termini_length(data)
    false_termini === nothing && return false
    trim5, trim3 = false_termini
    iszero(trim5) && iszero(trim3) && return false
    seqlen = UInt16(length(unwrap(data.seq)))


    toremove = Set{String}() # we remove these keys
    n_stripped = 0
    for (id, segment_data) in data
        false_termini =  false_termini_length(segment_data)
        false_termini === nothing && continue
        trim5, trim3 = false_termini
        iszero(trim5) && iszero(trim3) && continue
        seqlen = UInt16(length(unwrap(segment_data.seq[])))
        segment_data.seq[] = some(unwrap(segment_data.seq[])[trim5+1:seqlen-trim3])
        for (proti, protein) in enumerate(segment_data.proteins), (orfi, orf) in enumerate(protein.orfs)
            if trim5 ≥ first(orf) || (seqlen - trim3) < last(orf)
                push!(toremove, id)
                break
            else
                segment_data.proteins[proti].orfs[orfi] = (first(orf) - trim5):(last(orf) - trim5)
                n_stripped += 1
            end
        end
    end
    if !iszero(length(toremove))
        filter!(data) do (k, v)
            !in(k, toremove)
        end
        println("Removed $(length(toremove)) segments due to ORFs in false termini")
    end
    println("Stripped false termini off $(n_stripped) segment data")
end

# We allow to strip up to 100 bp in each end off
function false_termini_length(data::SegmentData)::Union{Nothing, Tuple{UInt16, UInt16}}
    seq = @unwrap_or data.seq (return nothing)
    trim5, trim3 = UInt16(0), UInt16(0)
    # Find true beginning of sequence, if it's within first 100 bp
    p = approxsearch(seq, dna"AGCAAAAGCAGG", 1)
    if !isempty(p)
        if first(p) < 100
            trim5 = UInt16(first(p) - 1)
        end
    end
    # Find true end of sequence
    p = approxrsearch(seq, dna"CTTGTTTCTCCT", 1)
    if !isempty(p)
        trim3 = UInt16(lastindex(seq) - last(p))
        if trim3 > 100
            trim3 = UInt16(0)
        end
    end
    return (trim5, trim3)
end

"""
    add_sequences!(segment_data::Dict{String, IncompleteSegmentData}, path::AbstractString)

Parse the FASTA file at `path`, and add sequences to each of the `IncompleteSegmentData` 
objects in `segment_data` found in the file.
Returns the number of updated `IncompleteSegmentData` objects.
"""
function add_sequences!(segment_data::Dict{String, IncompleteSegmentData}, path::AbstractString)
    maybe_gzip(path) do io
        add_sequences!(segment_data, io)
    end
end