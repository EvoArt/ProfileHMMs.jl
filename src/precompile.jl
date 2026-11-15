# Compile the search path at install time so the first search in a session
# (e.g. in VS Code) does not wait for code generation.
@setup_workload begin
    M = 40
    emit = fill(-log(0.01 / 3), 4, M)
    cons = "ACGTTGCAAGGTCCATGACTGATCGGATCCTAGCTAGGCA"
    for (k, c) in enumerate(cons)
        emit[findfirst(==(c), "ACGT"), k] = -log(0.99)
    end
    trans = repeat(-log.([0.9, 0.05, 0.05, 0.8, 0.2, 0.8, 0.2]), 1, M + 1)
    prof = HMMProfile("precompile", "", "", M, M + 10, "DNA", 1, 1.0, (-10.0, 0.7),
                      (-10.0, 0.7), (-10.0, 0.7), fill(log(4), 4), emit, fill(log(4), 4, M + 1), trans)
    @compile_workload begin
        h = build_plan7(prof)
        search_sequence(h, "NNNN" * cons * "NNNN"; evalue_threshold=1.0)
        path, io = mktemp()
        println(io, ">s\n", "NNNN", cons, "NNNN"); close(io)
        scan_fasta(h, path; progress=false, evalue_threshold=1.0)
        scan_fasta([h], path; progress=true, io=IOBuffer(), evalue_threshold=1.0)
        scan_fasta(h, path; progress=false, evalue_threshold=1.0, scoring=:forward)
        rm(path; force=true)
        sprint(show, MIME"text/plain"(), h)
    end
end
