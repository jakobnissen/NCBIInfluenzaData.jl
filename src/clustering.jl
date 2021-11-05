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
        ifilter(x -> startswith(x, '>')) |>
        imap(x -> String(strip(x)[2:end])) |>
        collect
    end

    return [data[i] for i in ids]
end

function run_cd_hit(path::AbstractString)
    outfile = path * ".cdhit"
    command = `cd-hit-est -i $path -o $outfile -aS 0.95 -n 9 -c 0.95 -d 32 -T 2`
    pipe = pipeline(command, stdout="$path.log")
    run(pipe)
    return outfile
end
