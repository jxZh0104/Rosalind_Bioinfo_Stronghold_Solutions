import re
import pandas as pd

# for reading files
def read_fasta(file, output = 'dict'):
    '''
    output = 'dict' (default): {Rosalind id: sequence}
    output = 'list': [sequence, ...] in file order
    '''
    with open(file, 'r') as f:
        raw = f.read().split('>')
        raw.pop(0)
    dnas = {}
    for entry in raw:
        entry = entry.replace('\n', '')
        match = re.search(r'Rosalind_\d+', entry)
        dnas[match.group(0)] = entry[(match.span()[1]):]
    if output == 'list':
        return list(dnas.values())
    return dnas

# for formatting outputs:
def index_transform(matched_indices):
    return " ".join([str(i+1) for i in matched_indices]) # python index starts from 0

def convert_tuples_to_strings(lst_of_tuples):
    for tu in lst_of_tuples:
        print(str(tu[0])+' '+str(tu[1]))

def convert_dict_to_strings(dictionary):
    for key in dictionary:
        matched_indices = dictionary[key]
        if len(matched_indices) > 0:
            print(key)
            print(index_transform(matched_indices))

# scoring matrices:

with open('data/scoring_matrices/BLOSUM62.txt', 'r') as f:
        raw_BLOSUM62 = f.read().split()
raw_BLOSUM62 = raw_BLOSUM62[20:]
BLOSUM62 = {}

for i in range(20):
    BLOSUM62[raw_BLOSUM62[21*i]] = [int(x) for x in raw_BLOSUM62[21*i+1:21*i+21]]
BLOSUM62 = pd.DataFrame(BLOSUM62, index = BLOSUM62.keys())

with open('data/scoring_matrices/PAM250.txt', 'r') as f:
        raw_PAM250 = f.read().split()
raw_PAM250 = raw_PAM250[20:]
PAM250 = {}

for i in range(20):
    PAM250[raw_PAM250[21*i]] = [int(x) for x in raw_PAM250[21*i+1:21*i+21]]
PAM250 = pd.DataFrame(PAM250, index = PAM250.keys())

# for HMM problems
def write_matrix_to_file(file, states, emit_chars, trans_mat, emit_mat):
     with open(file, 'w') as f:
        f.write('\t'+'\t'.join(states)+'\n')
        for state in states:
            f.write(state + '\t' + '\t'.join([str(round(x, 3)) if x != int(x) else str(x) for x in trans_mat.loc[state, :]]) + '\n')
        f.write('--------\n')

        f.write('\t'+'\t'.join(emit_chars)+'\n')
        for state in states:
            f.write(state + '\t' + '\t'.join([str(round(x, 3)) if x != int(x) else str(x) for x in emit_mat.loc[state, :]]) + '\n')
