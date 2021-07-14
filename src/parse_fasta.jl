function add_sequences!(segment_data::Dict{String, SegmentData}, io::IO)
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

"""
    add_sequences!(segment_data::Dict{String, SegmentData}, path::AbstractString)

Parse the FASTA file at `path`, and add sequences to each of the `SegmentData` 
objects in `segment_data` found in the file.
Returns the number of updated `SegmentData` objects.
"""
function add_sequences!(segment_data::Dict{String, SegmentData}, path::AbstractString)
    maybe_gzip(path) do io
        add_sequences!(segment_data, io)
    end
end