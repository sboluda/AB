from Bio import SeqIO
from Bio.Align import substitution_matrices
import sys

from build_index import build_index
from find_seeds import find_seeds
from extend_seeds import extend_seeds
from merge_overlapping import merge_overlapping


if __name__ == "__main__":

    k = 3  # k-mer length
    t = 11  # BLOSUM62 score threshold

    query_file = __  # complete
    db_file = __  # complete

    query = __  # complete
    db_seqs = __  # complete

    print(f"\nPreprocessing database (k={k})...")
    index = build_index(db_seqs, k)

    print(f"\nFinding seeds (threshold t={t})...")
    seeds = find_seeds(query, index, k, t)

    print(f"\nExtending seeds...")
    hits = extend_seeds(seeds, query, db_seqs, k)

    print(f"\nMerging overlapping hits...")
    merged = merge_overlapping(hits, query, db_seqs)

    print(f"\nMerged hits sorted by length (longest first): {len(merged)}")
    merged.sort(key=lambda h: h["query_end"] - h["query_start"], reverse=True)
    for h in merged[:10]:  # print top 10 longest hits
        diagonal = h["db_start"] - h["query_start"]
        print(
            f"  DB[{h['db_iseq']}] diag={diagonal}"
            f" Q[{h['query_start']}:{h['query_end']}]"
            f" -> DB[{h['db_start']}:{h['db_end']}]"
        )
        print(f"    Query: {h['query_seq']}")
        print(f"    DB:    {h['db_seq']}")
