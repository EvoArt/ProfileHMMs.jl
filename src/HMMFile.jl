"""
    HMMFile.jl

Parse HMMER3 profile HMM files (.hmm) into structured Julia types.
Supports multi-model files (e.g. barrnap's 16S + 23S + 5S).
"""

"""
    HMMProfile

Parsed representation of a single HMMER3 profile HMM.

Fields:
- `name`: model name (e.g. "16S_rRNA")
- `acc`: accession (e.g. "RF00177")
- `desc`: description
- `length`: model length (number of match states)
- `alphabet`: alphabet string ("RNA" or "DNA")
- `nseq`: number of sequences used to build the model
- `stats_msv`, `stats_viterbi`, `stats_forward`: (mu, lambda) for E-value calculation
- `compo`: background composition (log-odds, 4 values for ACGU/ACGT)
- `match_emissions`: Matrix{Float64} (K × M) — match emission scores per node
- `insert_emissions`: Matrix{Float64} (K × M) — insert emission scores per node
- `transitions`: Matrix{Float64} (7 × M+1) — transition scores per node
    Row order: m→m, m→i, m→d, i→m, i→i, d→m, d→d
    Column 1 = node 0 (COMPO), columns 2:end = nodes 1:M
"""
struct HMMProfile
    name::String
    acc::String
    desc::String
    length::Int
    maxl::Int
    alphabet::String
    nseq::Int
    effn::Float64
    stats_msv::Tuple{Float64, Float64}
    stats_viterbi::Tuple{Float64, Float64}
    stats_forward::Tuple{Float64, Float64}
    compo::Vector{Float64}
    match_emissions::Matrix{Float64}   # K × M
    insert_emissions::Matrix{Float64}  # K × (M+1), col 1 = node 0
    transitions::Matrix{Float64}       # 7 × (M+1), col 1 = node 0
end

function Base.show(io::IO, h::HMMProfile)
    print(io, "HMMProfile(\"$(h.name)\", M=$(h.length), alphabet=$(h.alphabet))")
end

"""
    parse_hmm_value(s::AbstractString) -> Float64

Parse a single HMM value. "*" maps to Inf (zero probability in log space).
"""
function parse_hmm_value(s::AbstractString)
    s == "*" ? Inf : parse(Float64, s)
end

"""
    parse_hmm_file(filepath::String) -> Vector{HMMProfile}

Parse a HMMER3 .hmm file, returning one `HMMProfile` per model in the file.
"""
function parse_hmm_file(filepath)
    lines = readlines(filepath)
    profiles = HMMProfile[]
    i = 1
    while i <= length(lines)
        if startswith(lines[i], "HMMER3")
            profile, i = _parse_single_hmm(lines, i)
            push!(profiles, profile)
        else
            i += 1
        end
    end
    isempty(profiles) && throw(ArgumentError("no HMMER3 profiles in $filepath"))
    profiles
end

function _parse_single_hmm(lines, start_idx)
    i = start_idx
    name = ""; acc = ""; desc = ""; model_len = 0; maxl = 0
    alphabet = ""; nseq = 0; effn = 0.0
    stats_msv = (0.0, 0.0); stats_vit = (0.0, 0.0); stats_fwd = (0.0, 0.0)

    # Parse header
    i += 1  # skip HMMER3 line
    while i <= length(lines) && !startswith(lines[i], "HMM")
        line = lines[i]
        if startswith(line, "NAME")
            name = strip(line[6:end])
        elseif startswith(line, "ACC")
            acc = strip(line[5:end])
        elseif startswith(line, "DESC")
            desc = strip(line[6:end])
        elseif startswith(line, "LENG")
            model_len = parse(Int, strip(line[6:end]))
        elseif startswith(line, "MAXL")
            maxl = parse(Int, strip(line[6:end]))
        elseif startswith(line, "ALPH")
            alphabet = strip(line[6:end])
        elseif startswith(line, "NSEQ")
            nseq = parse(Int, strip(line[6:end]))
        elseif startswith(line, "EFFN")
            effn = parse(Float64, strip(line[6:end]))
        elseif startswith(line, "STATS LOCAL MSV")
            parts = split(line)
            stats_msv = (parse(Float64, parts[end-1]), parse(Float64, parts[end]))
        elseif startswith(line, "STATS LOCAL VITERBI")
            parts = split(line)
            stats_vit = (parse(Float64, parts[end-1]), parse(Float64, parts[end]))
        elseif startswith(line, "STATS LOCAL FORWARD")
            parts = split(line)
            stats_fwd = (parse(Float64, parts[end-1]), parse(Float64, parts[end]))
        end
        i += 1
    end

    # Skip "HMM  A  C  G  U" header line
    i += 1
    # Skip "m->m m->i ..." header line
    i += 1

    uppercase(alphabet) in ("DNA", "RNA") || throw(ArgumentError("only DNA/RNA profiles are supported: $alphabet"))
    model_len > 0 || throw(ArgumentError("profile length must be positive"))
    i + 3model_len + 2 <= length(lines) || throw(ArgumentError("truncated HMM profile: $name"))
    K = 4

    # Parse COMPO (node 0)
    # Line 1: COMPO + K emission values (background composition)
    compo_parts = split(lines[i])
    compo = [parse_hmm_value(compo_parts[j]) for j in 2:K+1]
    i += 1

    # Line 2: insert emission for node 0
    insert0_parts = split(lines[i])
    insert0 = [parse_hmm_value(s) for s in insert0_parts[1:K]]
    i += 1

    # Line 3: transitions for node 0
    trans0_parts = split(lines[i])
    trans0 = [parse_hmm_value(s) for s in trans0_parts[1:7]]
    i += 1

    # Allocate arrays
    match_emissions = Matrix{Float64}(undef, K, model_len)
    insert_emissions = Matrix{Float64}(undef, K, model_len + 1)
    transitions = Matrix{Float64}(undef, 7, model_len + 1)

    insert_emissions[:, 1] = insert0
    transitions[:, 1] = trans0

    # Parse nodes 1..M
    for node in 1:model_len
        # Match emission line: node_num + K values + map + consensus + ...
        m_parts = split(lines[i])
        for k in 1:K
            match_emissions[k, node] = parse_hmm_value(m_parts[k + 1])
        end
        i += 1

        # Insert emission line: K values
        ins_parts = split(lines[i])
        for k in 1:K
            insert_emissions[k, node + 1] = parse_hmm_value(ins_parts[k])
        end
        i += 1

        # Transition line: 7 values
        tr_parts = split(lines[i])
        for t in 1:7
            transitions[t, node + 1] = parse_hmm_value(tr_parts[t])
        end
        i += 1
    end

    # Skip "//"
    while i <= length(lines) && strip(lines[i]) != "//"
        i += 1
    end
    i <= length(lines) || throw(ArgumentError("missing // terminator: $name"))
    i += 1  # skip the "//" line itself

    profile = HMMProfile(
        name, acc, desc, model_len, maxl, alphabet,
        nseq, effn, stats_msv, stats_vit, stats_fwd,
        compo, match_emissions, insert_emissions, transitions
    )
    return profile, i
end
