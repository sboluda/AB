from Bio.Align import substitution_matrices
subst_mat = substitution_matrices.load('BLOSUM62')

print(subst_mat["A", "A"])  # usually 4
print(subst_mat["A", "R"])  # usually -1