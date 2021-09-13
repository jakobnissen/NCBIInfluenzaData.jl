"""
    isok_all(::SegmentData)

Check if the `SegmentData` passes all criteria. Checks:
* Segment has a sequence
* Segment has less than 5 ambiguous bases
* `isok_minimum_proteins`
* `isok_seq_length` (inf A only, skipped for B)
* `isok_orf_length` (inf A only, skipped for B)
* `isok_translatable`
"""
function isok_all(data::SegmentData)::Bool
    count(isambiguous, data.seq) < 5 &&
    isok_segment_label(data)
    isok_minimum_proteins(data) &&
    isok_seq_length(data) &&
    isok_orf_length(data) &&
    isok_translatable(data)
end

"""
    isok_segment_label(::SegmentData)

Check whether a `SegmentData`'s DNA and protein variant is consistent, i.e. no PB2
protein in a PB1 segment. This may happen when segments gets mislabeled, e.g. 
KT853292.1 is labeled "segment 2", but is PB2 (segment 1).
"""
function isok_segment_label(data::SegmentData)::Bool
    all(p -> source(p.var) == data.segment, data.proteins)
end

"""
    isok_minimum_proteins(::SegmentData)

Checks whether a segment has the minimum number of ORFS to be considered "complete",
i.e. whether it does not miss any non-auxilliary proteins.
"""
function isok_minimum_proteins(data::SegmentData)::Bool
    segment = data.segment
    serotype = data.serotype
    minimum = if segment == Segments.PB2
        (Proteins.PB2,)
    elseif segment == Segments.PB1
        (Proteins.PB1,)
    elseif segment == Segments.PA
        (Proteins.PA,)
    elseif segment == Segments.HA
        (Proteins.HA,)
    elseif segment == Segments.NP
        (Proteins.NP,)
    elseif segment == Segments.NA
        # NB is auxilliary
        (Proteins.NA,)
    elseif segment == Segments.MP
        is_error(serotype) ? (Proteins.M1, Proteins.BM2) : (Proteins.M1, Proteins.M2)
    elseif segment == Segments.NS
        (Proteins.NS1, Proteins.NEP)
    else
        @assert false "Unreachable!" # This helps inference
    end
    return !any(minimum) do protein
        isnothing(findfirst(i -> i.var == protein, data.proteins))
    end
end

# This has been empirially determined by looking at the distributions
# in order to find outliers
const ACCEPTABLE_SEGMENT_LENGTHS = Dict(
    Segments.NP => 1490:1570,
    Segments.HA => 1680:1780,
    Segments.MP => 980:1030,
    Segments.PB1 => 2270:2345,
    Segments.PA => 2150:2240,
    Segments.NA => 1350:1470,  # very wide distribution??
    Segments.NS => 820:895,
    Segments.PB2 => 2280:2345
)

"""
    isok_seq_length(::SegmentData)

Check whether a segment's sequence length falls within a reasonable length interval.
Always returns `true` for Influenza B segments.
"""
function isok_seq_length(data::SegmentData)::Bool
    # The sizes for influenza B are a little different.
    is_error(data.serotype) && return true
    length(data.seq) in ACCEPTABLE_SEGMENT_LENGTHS[data.segment]
end

# This is pretty lax, so can be tightened later. It's just to remove
# obvious annotation errors.
const ACCEPTABLE_ORF_LENGTHS = Dict(
    Proteins.PB2 => 2250:2400,
    Proteins.PB1 => 2200:2400,
    Proteins.N40 => nothing, # nothings mean too little data to know range
    Proteins.PB1F2 => 0:350, # often truncated in real viruses
    Proteins.PA => 2100:2250,
    Proteins.PAX => 600:800,
    Proteins.HA => 1600:1800,
    Proteins.NP => 1450:1550,
    Proteins.NA => 1300:1450,
    Proteins.NB => nothing,
    Proteins.M1 => 720:780,
    Proteins.M2 => 250:325,
    Proteins.BM2 => nothing,
    Proteins.M42 => nothing,
    Proteins.NS1 => 600:750,
    Proteins.NEP => 325:400,
    Proteins.NS3 => nothing
)

# Yeah, some of these ORFs are _clearly_ not good enough.
"""
    isok_orf_length(::SegmentData)

Check whether a segment's ORF length falls within a reasonable length interval.
If the segment does not have a SeroType, returns `true.`.
"""
function isok_orf_length(data::SegmentData)::Bool
    is_error(data.serotype) && return true
    
    for protein in data.proteins
        len = sum(length, protein.orfs)
        range = ACCEPTABLE_ORF_LENGTHS[protein.var]
        if range !== nothing && !(len in range)
            return false
        end
    end
    return true
end

"""
    isok_translatable(::SegmentData)

Check whether a segment's ORFs are translatable:
* If the segment does not have a sequence or any ORFS, return `false`.
* If any ORF exceed past the segment sequence, return `false`
* If the joined ORFs' length is not divisible by 3, return `false`
* If any stops in sequence, return `false`.
* If sequence is not followed by a stop codon, return false
"""
function isok_translatable(data::SegmentData)::Bool
    nt_sequence = LongDNASeq()
    for protein in data.proteins
        if any(last(orf) > length(data.seq) for orf in protein.orfs)
            return false
        end
        join!(nt_sequence, (@view(data.seq[orf]) for orf in protein.orfs))

        # Must have a nonzerolength divisible by 3
        (length(nt_sequence) > 2 || iszero(length(nt_sequence) % 3)) || return false
        aa_sequence = BioSequences.translate(nt_sequence)

        # Must not contain stop
        findfirst(AA_Term, aa_sequence) === nothing || return false

        # There must be a stop just after the seq
        lastpos = protein.orfs[end][end]
        stoprange = lastpos+1:lastpos+3
        checkbounds(Bool, eachindex(data.seq), stoprange) || return false
        stopseq = data.seq[stoprange]
        (all(!isambiguous, stopseq) && Influenza.is_stop(DNACodon(stopseq))) || return false
    end
    return true
end

function join!(seq::BioSequence, it)
    len = sum(length, it)
    length(seq) != len && resize!(seq, len)
    offset = 0
    for i in it
        seq[offset+1:offset+length(i)] = i
        offset += length(i)
    end
    @assert offset == len
    return seq
end