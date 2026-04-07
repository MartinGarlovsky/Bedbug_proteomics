import subprocess
from pathlib import Path

# ========== USER CONFIGURATION ========== #
query_fasta = "cimex_slaps_clean.fasta"
target_fasta = "Drosophila_clean.fasta"
evalue_cutoff = 1e-2
min_identity = 30.0            # Minimum % identity
min_coverage = 50.0            # Minimum coverage of query
threads = 4
out_dir = Path("blast_results")
out_dir.mkdir(exist_ok=True)
# ======================================== #

def run(cmd, desc):
    print(f"\n🔧 {desc} ...")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"❌ ERROR: {desc}\n{result.stderr}")
        exit(1)
    print(f"✅ Done: {desc}")

def make_db(fasta, db_name):
    run(f"makeblastdb -in {fasta} -dbtype prot -parse_seqids -out {db_name}", f"Building BLAST DB for {db_name}")

def parse_blast(blast_file, fasta_lengths):
    """Returns a dict of best hits per query, filtered by identity and coverage"""
    best_hits = {}
    with open(blast_file) as f:
        for line in sorted(f, key=lambda x: float(x.split()[10])):  # sort by e-value
            parts = line.strip().split('\t')
            qid, sid = parts[0], parts[1]
            identity = float(parts[2])
            aln_len = int(parts[3])
            qlen = fasta_lengths.get(qid, 1)
            coverage = (aln_len / qlen) * 100
            if identity >= min_identity and coverage >= min_coverage:
                if qid not in best_hits:
                    best_hits[qid] = sid
    return best_hits

def get_seq_lengths(fasta_file):
    lengths = {}
    with open(fasta_file) as f:
        seq_id = None
        seq = []
        for line in f:
            if line.startswith(">"):
                if seq_id:
                    lengths[seq_id] = len("".join(seq))
                seq_id = line.strip().lstrip(">").split()[0]
                seq = []
            else:
                seq.append(line.strip())
        if seq_id:
            lengths[seq_id] = len("".join(seq))
    return lengths

# STEP 1: Make BLAST DBs
make_db(query_fasta, "query_db")
make_db(target_fasta, "target_db")

# STEP 2: Run Forward BLAST
forward_out = out_dir / "forward.blast"
run(f"blastp -query {query_fasta} -db target_db -evalue {evalue_cutoff} -outfmt 6 -num_threads {threads} -max_target_seqs 1000 -out {forward_out}", "Forward BLAST")

# STEP 3: Extract best forward hits with filtering
query_lengths = get_seq_lengths(query_fasta)
forward_best = parse_blast(forward_out, query_lengths)
print(f"✅ Forward best hits passing filters: {len(forward_best)}")

# STEP 4: Extract target sequences for reverse BLAST
target_ids_file = out_dir / "target_ids.txt"
with open(target_ids_file, "w") as f:
    for sid in set(forward_best.values()):
        f.write(sid + "\n")

target_best_fasta = out_dir / "target_best.fasta"
run(f"blastdbcmd -db target_db -entry_batch {target_ids_file} -out {target_best_fasta}", "Extracting target best hits")

# STEP 5: Run Reverse BLAST
reverse_out = out_dir / "reverse.blast"
run(f"blastp -query {target_best_fasta} -db query_db -evalue {evalue_cutoff} -outfmt 6 -num_threads {threads} -max_target_seqs 1000 -out {reverse_out}", "Reverse BLAST")

# STEP 6: Extract best reverse hits
target_lengths = get_seq_lengths(target_best_fasta)
reverse_best = parse_blast(reverse_out, target_lengths)
print(f"✅ Reverse best hits passing filters: {len(reverse_best)}")

# STEP 7: Identify Reciprocal Best Hits
rbh_file = out_dir / "reciprocal_best_hits.tsv"
with open(rbh_file, "w") as out:
    out.write("Cimex_ID\tDrosophila_ID\n")
    count = 0
    for q, s in forward_best.items():
        if reverse_best.get(s) == q:
            out.write(f"{q}\t{s}\n")
            count += 1
print(f"\n🎯 Reciprocal Best Hits found: {count} / {len(forward_best)} queries")
print(f"📄 Results saved to: {rbh_file}")

# STEP 8: Summary
with open(out_dir / "summary.txt", "w") as f:
    f.write(f"Total query proteins: {len(query_lengths)}\n")
    f.write(f"Forward best hits passing filters: {len(forward_best)}\n")
    f.write(f"Reverse best hits passing filters: {len(reverse_best)}\n")
    f.write(f"Reciprocal best hits: {count}\n")
