# In this file, code to cluster (deduplicate) and serialize the sequences
# are put.

"""
    cd_hit_deduplicate_all(data::Dict{String, SegmentData})
        -> Tuple{Dict{String, SegmentData}, String}

Deduplicate the values in `data`, on a per-segment basis in parallel, using
CD-HIT. To work, Julia must have `CD-HIT` on its `PATH` variable, i.e.
`run(cd-hit-est)` must work.

Returns a deduplicated dict, and the path to the directory of the CD-HIT.
"""
function cd_hit_deduplicate_all(
    data::Dict{String, SegmentData},
    tmpdir::String=mktempdir()
)
    # Collect by segment
    bysegment = Dict{Segment, Dict{String, SegmentData}}()
    for dat in values(data)
        get!(Dict{String, SegmentData}, bysegment, dat.segment)[dat.id] = dat
    end

    results = Vector{Vector{SegmentData}}(undef, length(bysegment))
    Threads.@threads for (i, dat) in collect(enumerate(values(bysegment)))
        results[i] = cd_hit_deduplicate(dat, tmpdir)
    end

    result = Dict{String, SegmentData}()
    for i in results, dat in i
        result[dat.id] = dat
    end
    return (result, tmpdir)
end

"Deduplicate using CD-hit"
function cd_hit_deduplicate(
    data::Dict{String, SegmentData},
    tmpdir::String=mktempdir()
)::Vector{SegmentData}

    # Create FASTA file
    isdir(tmpdir) || error("Directory not found $tmpdir")
    fasta_path, file_io = mktemp(tmpdir)
    writer = FASTA.Writer(file_io)
    for i in values(data)
        write(writer, FASTA.Record(i.id, i.seq))
    end
    close(writer)

    # Run CD hit
    cd_path = run_cd_hit(fasta_path)

    # Read results to vector
    ids = open(cd_path) do io
        eachline(io) |>
        Filter(x -> startswith(x, '>')) ⨟
        Map(x -> String(strip(x)[2:end])) |>
        collect
    end

    return [data[i] for i in ids]
end

function run_cd_hit(path::AbstractString)
    outfile = path * ".cdhit"
    command = `cd-hit-est -i $path -o $outfile -aS 0.9 -c 0.95 -d 32`
    pipe = pipeline(command, stdout="$path.log")
    run(pipe)
    return outfile
end

"""
    serialize_segments(::Dict{String, SegmentData}, outdir::String)

In `outdir`, for each segment present in the dict, create a `.jls` serialized
file with the ORFs, and a FASTA file with the sequences.
"""
function serialize_segments(
    deduplicated::Dict{String, SegmentData},
    outdir::String
)
    isdir(outdir) || error("Dir not found: \"$outdir\"")

    # Split to segments
    bysegment = Dict{Segment, Vector{SegmentData}}()
    for data in values(deduplicated)
        push!(get!(Vector{SegmentData}, bysegment, data.segment), data)
    end

    # Serialize
    for (segment, datavec) in bysegment
        jls_path = joinpath(outdir, string(segment) * ".jls")
        serialize_orfs(datavec, jls_path)
        fasta_path = joinpath(outdir, string(segment) * ".fna")
        serialize_fasta(datavec, fasta_path)
    end
end

function serialize_orfs(datavec::Vector{SegmentData}, path::String)
    #                 Protein  Vector of orfs
    ProteinType = Tuple{UInt8, Vector{UnitRange{UInt16}}}
    #                accession  vector of proteins
    SegmentType = Tuple{String, Vector{ProteinType}}
    result = Vector{SegmentType}()
    
    for data in datavec
        vec = Vector{ProteinType}()
        segment_repr = (data.id, vec)
        push!(result, segment_repr)
        for protein in data.proteins
            var_uint8 = reinterpret(UInt8, protein.variant)
            push!(vec, (var_uint8, protein.orfs))
        end
    end
    serialize(path, result)
end

function serialize_fasta(datavec::Vector{SegmentData}, fasta_path::String)
    open(FASTA.Writer, fasta_path) do writer
        for data in datavec
            write(writer, FASTA.Record(data.id, data.seq))
        end
    end
end