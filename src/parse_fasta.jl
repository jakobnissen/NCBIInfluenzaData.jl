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

function parse_annotation_refs(
    io::IO,
    jls_path::AbstractString,
    fna_path::AbstractString
)::Vector{Reference}
    # First load to Assembly
    reader = FASTA.Reader(io)
    asms = Influenza.Assembly[]
    for record in reader
        id = FASTA.identifier(record)
        isnothing(id) && error("No identifier in record")
        # Skip inf C
        if (m = match(r"^([C])-seg([1-7])$", id); m) !== nothing # Inf C
            continue
        elseif (m = match(r"^([AB])-seg([1-8])(?:_([HN]\d+))?$", id); m) !== nothing # A/B
            species = something(m[1])
            segment = Segment(parse(UInt8, something(m[2])) - 1)
            suffix = if (s = m[3]; s) !== nothing
                @assert segment === (first(s) === 'H' ? Segments.HA : Segments.NA)
                s
            else    
                ""
            end
            name = species * suffix * '_' * string(segment)
            push!(asms, Influenza.Assembly(name, FASTA.sequence(LongDNASeq, record)))
        else
            error("Cannot parse header: \"", id, "\"")
        end
    end
    alnasms = map(unwrap, annotate(asms, jls_path, fna_path))
    return map(alnasms) do alnasm
        prots = map(alnasm.proteins) do protein
            Influenza.ReferenceProtein(
                protein.variant,
                unwrap(protein.orfs)
            )
        end::Vector{Influenza.ReferenceProtein}
        Reference(alnasm.assembly.name, alnasm.reference.segment, alnasm.assembly.seq, prots)
    end
end