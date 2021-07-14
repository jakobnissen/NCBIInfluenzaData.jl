
using Pkg
Pkg.activate(".")
#Pkg.instantiate()

using BioSequences # commit 252720879620716f5f33a54e53ffcd9716706636
using CodecZlib
using ErrorTypes
using FASTX # commit e9ff73b65e6638ed740b2fce5a01b254ecb24c55
using Serialization
using Transducers
using UnicodePlots
using Influenza
import Downloads

using Base: RefValue

"Complete the entire workflow. Get an overview of all that happens here."
function main()
    println("Starting with $(Threads.nthreads()) threads.")

    # Get the data from the internet (except if it already exists)
    download_influenza_data(force=false)

    # Clean the data
    isdir("results") || mkdir("results")
    clean_the_data(joinpath("results", "genomeset.clean.dat"), joinpath("download", "genomeset.dat.gz"))

    # Parse genomeset to SegmentData with all basic info
    segment_data = open(parse_cleaned_genomeset, joinpath("results", "genomeset.clean.dat"))
    
    # Load influenza.fna and update the seq field of the SegmentData
    maybe_gzip(joinpath("download", "influenza.fna.gz")) do io
        add_sequences!(segment_data, io)
    end

    # Load influenza_aa.dat.gz to get information on what kind of protein
    # the different ORFs make
    accession_protein_map = maybe_gzip(parse_inf_aa, joinpath("download", "influenza_aa.dat.gz"))

    # Load influenza.dat and update the ORF fields of SegmentData
    maybe_gzip(joinpath("download", "influenza.dat.gz")) do io
        add_orfs!(segment_data, accession_protein_map, io)
    end

    # We remove the human sequences, because we want to add in a set of manually
    # curated human sequences. Also remove all species in we don't need.
    filter!(pair -> last(pair).host ∈ (avian, swine), segment_data)

    # Add extra records - these are annotated by The NCBI Influenza Virus Sequence Annotation Tool
    # and manually curated to have the exact right format of the header.
    # These are the human sequences that should NOT be deduplicated
    add_extra_records!(segment_data, joinpath("raw", "annotation.tbl.gz"), joinpath("raw", "seqs.fna.gz"), false)

    # Do the exact same thing for extra, user-input sequences. I do this in order to have a homemade
    # file where I can add rare, unseen sequences to. Should I perhaps upload them?
    # These ones are deduplicated
    add_extra_records!(segment_data, joinpath("raw", "annotation.extra.tbl"), joinpath("raw", "seqs.extra.fna"), true)

    # Remove sequencing artifacts beyond the true termini of influenza sequences.
    strip_false_termini!(segment_data)

    # Series of filters on the SegmentData
    filter_segment_data!(segment_data)

    # Cluster using CD-hit for 95% identity
    # This is a Dict{Segment, Set{String}}, where String == accessions
    deduplicated = cd_hit_deduplicate(segment_data)

    # Serialize to files
    serialize_segments(segment_data, deduplicated)
end

function download_influenza_data(;force=false)
    isdir("download") && !force && return nothing
    isdir("download") || mkdir("download")
    ftp_address = "https://ftp.ncbi.nih.gov/genomes/INFLUENZA/"
    for filename in [
        "genomeset.dat.gz",
        "influenza.dat.gz",
        "influenza.fna.gz",
        "influenza_aa.dat.gz"
        ]
        Downloads.download(joinpath(ftp_address, filename), joinpath("download", filename))
        println("Downloaded $filename")
    end
end

function clean_the_data(outpath, inpath)
    file = open(inpath)
    outfile = open(outpath, "w")
    decompressed = GzipDecompressorStream(file)

    for line in eachline(decompressed) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        
        # Check correct number of fields
        if length(fields) != 11
            println("\"$line\"")
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
        nn = parse_name(name)
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

"Parses a name in genomeset.dat during cleaning, returning none if it's malformed.
If it's unexpectedly malformed, throw an error"
function parse_name(s::Union{String, SubString})::Option{String}
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

function parse_from_integer(::Type{Segment}, s::AbstractString)::Option{Segment}
    y = tryparse(UInt8, s)
    y === nothing && return none
    (iszero(y) | (y > 0x08)) && return none
    some(reinterpret(Segment, y - 0x01))
end

@enum Species::UInt8 human swine avian other

function Base.parse(::Type{Species}, s::Union{String, SubString})::Species
    lcase = lowercase(s)
    lcase == "human" && return human
    lcase == "swine" && return swine
    lcase == "avian" && return avian
    return other
end

struct ProteinORF
    variant::Protein
    orfs::Vector{UnitRange{UInt16}}
end

struct SegmentData
    id::String
    host::Species
    segment::Segment
    subtype::SubType
    deduplicate::Bool
    clade::Option{String}
    year::Int16
    name::String
    proteins::Vector{ProteinORF}
    seq::RefValue{Option{LongDNASeq}}
end

function parse_cleaned_genomeset(io::IO)::Dict{String, SegmentData}
    result = Dict{String, SegmentData}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        data = @unwrap_or parse(SegmentData, line) continue
        gi = data.id
        @assert !haskey(result, gi) "GB identifier $(gi) not unique"
        result[gi] = data
    end
    result
end

function Base.parse(::Type{SegmentData}, line::Union{String, SubString{String}})::Option{SegmentData}
    fields = split(line, '\t')
    subtype = @? parse(SubType, fields[4])
    year = @? parse_year(fields[6])
    host = parse(Species, fields[2])
    gi = String(fields[1])
    name = String(fields[8])
    segment = @? parse_from_integer(Segment, fields[3])
    some(SegmentData(gi, host, segment, subtype, true, none(String), year, name, ProteinORF[], Ref(none(LongDNASeq))))
end

function parse_year(s::Union{String, SubString})::Option{Int16}
    isempty(s) && return none
    pos_slash_found = findfirst(isequal('/'), s)
    last_byte = pos_slash_found === nothing ? ncodeunits(s) : pos_slash_found - 1
    return some(parse(Int16, SubString(s, 1:last_byte)))
end

function add_sequences!(segment_data::Dict{String, SegmentData}, io::IO)
    n_updated = 0
    for record in FASTA.Reader(io)
        # headers begins with e.g. gi|59292|gb|X53029|Influenza A virus
        gb_accession = split(FASTA.identifier(record)::String, '|')[4]
        haskey(segment_data, gb_accession) || continue
        data = segment_data[gb_accession]
        @assert is_error(data.seq[])
        data.seq[] = some(FASTA.sequence(LongDNASeq, record))
        n_updated += 1
    end
    println("Updated $n_updated/$(length(segment_data)) records")
    segment_data
end

function parse_inf_aa(io::IO)::Dict{String, Protein}
    result = Dict{String, Protein}()
    for line in eachline(io) |> Map(strip) ⨟ Filter(!isempty)
        fields = split(line, '\t')
        accession = first(fields)
        variant = @unwrap_or parse(Protein, fields[3]) continue
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
    println("Updated $n_updates/$(length(segment_data)) records")
    segment_data
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
    if p2 === nothing
        n = parse(UInt16, s)
        n:n
    else
        start = parse(UInt16, @view s[1:p2-1])
        stop = parse(UInt16, s[p2+1:ncodeunits(s)])
        start:stop
    end
end

function filter_segment_data!(segment_data::Dict{String, SegmentData})
    len = length(segment_data)
    n_human = count(v -> v.host == human, values(segment_data))

    # Must have a nucleotide sequence
    filter!(segment_data) do (name, data)
        !is_error(data.seq[])
    end
    println("No nuc sequence: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Must be a minimum number of ORFs depending on segment an inf A/B
    filter!(segment_data) do (name, data)
        contains_minimum_proteins(data)
    end
    println("Minimum ORFs: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Must not be more than 4 ambiguous bases
    filter!(segment_data) do (name, data)
        count(isambiguous, unwrap(data.seq[])) < 5
    end
    println("Ambiguous bases: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Filter for outliers of lengths
    @warn "Not filtering influenza B for segment lengths"
    filter!(segment_data) do (name, data)
        has_acceptable_seq_length(data)
    end
    println("Nucleotide lengths: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # Filter for outliers of lengths of proteins
    filter!(segment_data) do (name, data)
        has_acceptable_orf_length(data)
    end
    println("ORF lengths: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # All ORFs can be translated and has stop exactly at end
    filter!(segment_data) do (name, data)
        is_all_translatable(data)
    end
    println("Translation: Removed $(len - length(segment_data)) sequences")
    len = length(segment_data)

    # All our human sequences should be good quality, so if we filtered away
    # any human sequences, we did something wrong.
    @assert n_human == count(v -> v.host == human, values(segment_data))
    
    return segment_data
end

isalpha(x::SubType) = x.data isa SubTypes.InfluenzaA
isbeta(x::SubType) = x.data isa Union{SubTypes.Victoria, SubTypes.Yamagata}

function contains_minimum_proteins(data::SegmentData)::Bool
    segment = data.segment
    subtype = data.subtype
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
        isalpha(subtype) ? (Proteins.M1, Proteins.M2) : (Proteins.M1, Proteins.BM2)
    elseif segment == Segments.NS
        (Proteins.NS1, Proteins.NEP)
    else
        @assert false "Unreachable!" # This helps inference
    end
    return !any(minimum) do protein
        isnothing(findfirst(i -> i.variant == protein, data.proteins))
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

function has_acceptable_seq_length(data::SegmentData)::Bool
    # The sizes for influenza B are a little different.
    data.subtype.data isa Union{SubTypes.Victoria, SubTypes.Yamagata} && return true
    length(unwrap(data.seq[])) in ACCEPTABLE_SEGMENT_LENGTHS[data.segment]
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
function has_acceptable_orf_length(data::SegmentData)::Bool
    data.subtype.data isa Union{SubTypes.Victoria, SubTypes.Yamagata} && return true
    
    for protein in data.proteins
        len = sum(length, protein.orfs)
        range = ACCEPTABLE_ORF_LENGTHS[protein.variant]
        if range !== nothing && !(len in range)
            return false
        end
    end
    return true
end

function is_all_translatable(data::SegmentData)::Bool
    seq = unwrap(data.seq[])
    aa_sequence = LongAminoAcidSeq()
    nt_sequence = LongDNASeq()
    for protein in data.proteins
        if any(last(orf) > length(seq) for orf in protein.orfs)
            println(data.id)
            error()
        end
        join!(nt_sequence, (@view(seq[orf]) for orf in protein.orfs))

        # Must have a length divisible by 3
        iszero(length(nt_sequence) % 3) || return false
        BioSequences.translate!(aa_sequence, nt_sequence)
        stop_pos = findfirst(AA_Term, aa_sequence)

        # Must have a stop exactly at the end, nowhere else
        stop_pos == lastindex(aa_sequence) || return false
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

# This is pretty primitive, but fuck it
function maybe_gzip(f, path::AbstractString)
    if endswith(path, ".gz")
        io = GzipDecompressorStream(open(path))
        try
            f(io)
        finally
            close(io)
        end
    else
        open(f, path)
    end
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

function fill_orfs_extra_records!(records::Dict{String, SegmentData}, annotpath::String)
    chunks = maybe_gzip(annotpath) do file
        partition_chunks(file)
    end
    filled_ids = Set{String}()
    for chunk in chunks
        parse_annot_segment_data!(records, filled_ids, chunk)
    end
    @assert length(records) == length(filled_ids)
    records
end

# Please have a look at the files in annotation.tbl to understand this
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

##=
# Some reference sequences are too long. True influenza sequences ends with CTTGTTTCTACT
# and begins with AGCAAAAGCAGG. (with perhaps a single substitution.) We can
# allow sequences that do not include these ends, but we cannot allow sequences
# that include artifacts flanking these true termini.
# Example: KJ476677 contains 19 bp 3' of true 3' terminal.
function strip_false_termini!(data::Dict{String, SegmentData})
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

# We allow to strip up to 50 bp in each end off
function false_termini_length(data::SegmentData)::Union{Nothing, Tuple{UInt16, UInt16}}
    seq = @unwrap_or data.seq[] (return nothing)
    trim5, trim3 = UInt16(0), UInt16(0)
    # Find true beginning of sequence, if it's within first 50 bp
    p = approxsearch(seq, dna"AGCAAAAGCAGG", 1)
    if !isempty(p)
        if first(p) < 50
            trim5 = UInt16(first(p) - 1)
        end
    end
    # Find true end of sequence
    p = approxrsearch(seq, dna"CTTGTTTCTCCT", 1)
    if !isempty(p)
        trim3 = UInt16(lastindex(seq) - last(p))
        if trim3 > 50
            trim3 = UInt16(0)
        end
    end
    return (trim5, trim3)
end

"Deduplicate all non-human seqs and add the human ones back"
function cd_hit_deduplicate(segment_data::Dict{String, SegmentData}
)::Dict{Segment, Set{String}}
    n_non_dedup = count(i -> !(i.deduplicate), values(segment_data))
    @warn "Not deduplicating $n_non_dedup sequences marked as non-duplicate"

    subdir = joinpath("results", "cdhit")
    isdir(subdir) || mkdir(subdir)

    duplicatable = filter(segment_data) do (key, data)
        data.deduplicate
    end

    println("Deduplicating sequences...")
    deduplicated = cd_hit_deduplicate(collect(values(duplicatable)), subdir)

    # Add in the ids of non-duplicated back
    for (key, data) in segment_data
        if !data.deduplicate
            push!(deduplicated[data.segment], data.id)
        end
    end

    deduplicated
end

function cd_hit_deduplicate(data::Vector{SegmentData}, dirname::String
)::Dict{Segment, Set{String}}
    # Collect FASTAS by segment
    bysegment = Dict(s => FASTA.Record[] for s in instances(Segment))
    for segmentdata in data
        record = FASTA.Record(segmentdata.id, unwrap(segmentdata.seq[]))
        push!(bysegment[segmentdata.segment], record)
    end

    # Write FASTA paths
    for (segment, seqs) in bysegment
        fasta_path = joinpath(dirname, string(segment) * ".fna")
        open(FASTA.Writer, fasta_path) do writer
            for seq in seqs
                write(writer, seq)
            end
        end
    end

    # Deduplicate
    fasta_paths = [joinpath(dirname, string(segment) * ".fna") for segment in instances(Segment)]
    run_cd_hit(fasta_paths)

    # Load in the deduplicated ones
    result = Dict(s => Set{String}() for s in instances(Segment))
    for segment in instances(Segment)
        result[segment] = open(joinpath(dirname, string(segment) * ".fna.cdhit")) do file
            eachline(file) |>
            Filter(x -> startswith(x, '>')) ⨟
            Map(x -> x[2:end]) |>
            Set
        end
    end

    result
end

function run_cd_hit(paths::Vector{String})
    Threads.@threads for path in paths
        command = `bin/cd-hit-est -i $path -o $path.cdhit -aS 0.9 -c 0.95 -d 32`
        pipe = pipeline(command, stdout="$path.log")
        run(pipe)
    end
end

function serialize_segments(segment_data::Dict{String, SegmentData},
    deduplicated::Dict{Segment, Set{String}}
)
    serialize_subtypes(segment_data, deduplicated)

    serial_dir = joinpath("results/deduplicated")
    isdir(serial_dir) || mkdir(serial_dir)

    # Serialize the ORFs itself
    for (segment, accessions) in deduplicated
        jls_path = joinpath(serial_dir, string(segment) * ".jls")
        serialize_orfs(segment_data, accessions, jls_path)
    end

    # Write the deduplicated FASTA files
    for (segment, accessions) in deduplicated
        fasta_path = joinpath(serial_dir, string(segment) * ".fna")
        open(FASTA.Writer, fasta_path) do writer
            for accession in accessions
                data = segment_data[accession]
                record = FASTA.Record(data.id, unwrap(data.seq[]))
                write(writer, record)
            end
        end
    end
end

# Serialize subtypes into a vector of (id, subtype)
function serialize_subtypes(segment_data::Dict{String, SegmentData},
    deduplicated::Dict{Segment, Set{String}}
)
    v = Tuple{String, SubType}[]
    for (s, ids) in deduplicated
        append!(v, [(s, segment_data[s].subtype) for s in ids])
    end
    open(io -> serialize(io, v), "results/subtypes.jls", "w")
end

function serialize_orfs(segment_data::Dict{String, SegmentData},
    accessions::Set{String}, path::String)
    #                 Protein  Vector of orfs
    ProteinType = Tuple{UInt8, Vector{UnitRange{UInt16}}}
    #                accession  vector of proteins
    SegmentType = Tuple{String, Vector{ProteinType}}
    result = Vector{SegmentType}()
    
    for accession in accessions
        data = segment_data[accession]
        vec = Vector{ProteinType}()
        segment_repr = (accession, vec)
        push!(result, segment_repr)
        for protein in data.proteins
            var_uint8 = reinterpret(UInt8, protein.variant)
            push!(vec, (var_uint8, protein.orfs))
        end
    end
    serialize(path, result)
end


isinteractive() || main()

  