"""Sequence-alignment problems (edit distance, global/local/semiglobal/overlap/multiple alignment).

Rosalind: lcsq, scsp, edit, edta, ctea, glob, gaff, gcon, loca, laff, smgb, oap, sims, mult, mgap, osym.
laff, smgb and oap also have faster Java versions in java/; ksim exists only in Java.

Run from the repository root; `python python/alignment.py` executes the demo calls at the bottom of this file.
"""

from itertools import combinations

import numpy as np

from utils import BLOSUM62, PAM250, read_fasta
from strings import hamming_dist


# http://rosalind.info/problems/lcsq/
# keeps a table of lengths next to a table of strings (see lcs_string_table for the simpler version used by scsp)
def lcs_with_length_table(s, t):
    longest_string = [['']*len(t) for i in range(len(s))]
    longest_string_length = [[0]*len(t) for i in range(len(s))]

    for i in range(len(s)):
        if i == 0:
            for j in range(len(t)):
                if s[i] == t[j]:
                    longest_string_length[i][j] = 1
                    longest_string[i][j] = s[0]
                else:
                    if j > 0:
                        longest_string_length[i][j] = longest_string_length[i][j-1]
                        longest_string[i][j] = longest_string[i][j-1]
        else:
            for j in range(len(t)):
                if j == 0:
                    if s[i] == t[j]:
                        longest_string_length[i][j] = 1
                        longest_string[i][j] = t[0]
                    else:
                        longest_string_length[i][j] = longest_string_length[i-1][j]
                        longest_string[i][j] = longest_string[i-1][j]
                else:
                    candidate_lengths = {(i-1, j-1): 0, (i-1, j): 0, (i, j-1): 0}
                    candidate_strings = {(i-1, j-1): '', (i-1, j): '', (i, j-1): ''}
                    if s[i] == t[j]:
                        candidate_lengths[(i-1, j-1)] = 1
                        candidate_strings[(i-1, j-1)] = s[i]
                    for index_s, index_t in candidate_lengths.keys():
                        candidate_lengths[(index_s, index_t)] += longest_string_length[index_s][index_t]
                        candidate_strings[(index_s, index_t)] = longest_string[index_s][index_t]+candidate_strings[(index_s, index_t)]
                    max_index = max(candidate_lengths.keys(), key = lambda x: candidate_lengths[x])
                    longest_string_length[i][j] = candidate_lengths[max_index]
                    longest_string[i][j] = candidate_strings[max_index]
    return longest_string[-1][-1]


# http://rosalind.info/problems/edit/
def edit_distance(s, t):
    min_edit = [[0]*len(t) for i in range(len(s))]
    for i in range(len(s)):
        if i == 0:
            for j in range(len(t)):
                min_edit[i][j] = 1-int(s[i] in t[:j+1]) + j # if s[i] in t[:j+1], return 0+j, else, return 1+j

        else:
            for j in range(len(t)):
                if j == 0:
                    min_edit[i][j] = 1-int(t[j] in s[:i+1]) + i
                else:
                    if s[i] == t[j]:
                        min_edit[i][j] = min_edit[i-1][j-1]
                    else:
                        min_edit[i][j] = min([min_edit[i-1][j-1]+1,
                                              min_edit[i-1][j]+1,
                                              min_edit[i][j-1]+1
                                             ])
    print(min_edit[-1][-1])
    return min_edit


# http://rosalind.info/problems/edta/
def edit_reconstruct(s, t):
    '''
    it is easier to reconstruct according to the min_edit table than
    reconstructing during creating the min_edit table
    '''
    min_edit = edit_distance(s, t)
    s_reconstruct = ''
    t_reconstruct= ''
    i = len(s) - 1
    j = len(t) - 1
    while i > 0 and j > 0:
        if min_edit[i][j] == min_edit[i-1][j-1] + 1-int(s[i] == t[j]):
            s_reconstruct = s[i] + s_reconstruct
            t_reconstruct = t[j] + t_reconstruct
            i -= 1
            j -= 1
        elif min_edit[i][j] == min_edit[i][j-1] + 1:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j] + t_reconstruct
            j -= 1
        else: # min_edit[i][j] == min_edit[i-1][j] + 1
            s_reconstruct = s[i] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i -= 1
    if i == 0 and j == 0:
        s_reconstruct = s[0] + s_reconstruct
        t_reconstruct = t[0] + t_reconstruct
    elif i == 0:
        if s[0] in t[:j+1]:
            match_index = t[:j+1].index(s[0])
            s_reconstruct = '-'*match_index + s[0] + '-'*(j-match_index) + s_reconstruct
            # insert s[0] at exactly the place where it matches t
        else:
            s_reconstruct = '-'*j + s[0] + s_reconstruct
        t_reconstruct = t[:j+1] + t_reconstruct
    else:
        if t[0] in s[:i+1]:
            match_index = s[:i+1].index(t[0])
            t_reconstruct = '-'*match_index + t[0] + '-'*(i-match_index) + t_reconstruct
        else:
            t_reconstruct = '-'*i + t[0] + t_reconstruct
        s_reconstruct = s[:i+1] + s_reconstruct

    print(s_reconstruct)
    print(t_reconstruct)
    with open('data/outputs/result.txt', 'w+') as f:
        f.writelines([str(min_edit[-1][-1])+'\n', s_reconstruct+'\n', t_reconstruct+'\n'])
    assert (not any([s_reconstruct[i] == '-' and t_reconstruct == '-' for i in range(len(s_reconstruct))]))
    assert (hamming_dist(s_reconstruct, t_reconstruct) == min_edit[-1][-1])


# http://rosalind.info/problems/glob/

def global_alignment_score_linear(s, t, scoring_matrix = BLOSUM62, gap_penalty = 5):
    '''
    linear gap penalty of -5 per gap
    '''
    max_score = [[0]*(len(t)+1) for i in range((len(s)+1))]
    # max_score[i][j] stores the max score when aligning s[:i] against s[:j]
    # hence max_score[0][0] is always 0

    for i in range(len(s)+1):
        max_score[i][0] = - i * gap_penalty

    for j in range(len(t)+1):
        max_score[0][j] = - j * gap_penalty

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            max_score[i][j] = max([max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]],
                                              max_score[i-1][j] - gap_penalty,
                                              max_score[i][j-1] - gap_penalty
                                             ])
    print(max_score[-1][-1])
    return max_score


# http://rosalind.info/problems/gaff/
def global_alignment_affine(s, t, scoring_matrix = BLOSUM62, start_penalty = 11, extend_penalty = 1):
    '''
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
    returns not only score but also alignment results
    '''
    Sn = [[0]*(len(t)+1) for i in range(len(s)+1)] # no gap at both ends
    Ss = [[0]*(len(t)+1) for i in range(len(s)+1)] # gap at end of s
    St = [[0]*(len(t)+1) for i in range(len(s)+1)] # gap at end of t
    # a gap cannot be at both the end of s and the end of t

    Pn = [[0]*(len(t)+1) for i in range(len(s)+1)]
    Ps = [[0]*(len(t)+1) for i in range(len(s)+1)]
    Pt = [[0]*(len(t)+1) for i in range(len(s)+1)]
    # each P matrix traces the matrix of the previous indices that leads to the max score at the current index
    # 0 represents the matrix Sn
    # 1 represents the matrix Ss
    # 2 represents the matrix St

    lst_match_matrices = [Sn, Ss, St]
    lst_prev_matrices = [Pn, Ps, Pt]

    for i in range(1, len(s)+1):
        Sn[i][0] = -float('inf') # we try to match s[:i] with "", can only match xxxxxxx with -------
        Ss[i][0] = -float('inf')
        St[i][0] = -start_penalty -(i-1)*extend_penalty

    for j in range(1, len(t)+1):
        Sn[0][j] = -float('inf')
        Ss[0][j] = -start_penalty -(i-1)*extend_penalty
        St[0][j] = -float('inf')

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            curr_score = scoring_matrix.loc[s[i-1], t[j-1]]
            alignment_scores_no_gap = {0: Sn[i-1][j-1] + curr_score,
                                       1: Ss[i-1][j-1] + curr_score,
                                       2: St[i-1][j-1] + curr_score}
            alignment_scores_gap_s = {0: Sn[i][j-1] - start_penalty,
                                      2: St[i][j-1] - start_penalty,
                                      1: Ss[i][j-1] - extend_penalty}
            alignment_scores_gap_t = {0: Sn[i-1][j] - start_penalty,
                                      1: Ss[i-1][j] - start_penalty,
                                      2: St[i-1][j] - extend_penalty}

            Pn[i][j] = max([0,1,2], key = lambda x: alignment_scores_no_gap[x])
            Ps[i][j] = max([0,1,2], key = lambda x: alignment_scores_gap_s[x])
            Pt[i][j] = max([0,1,2], key = lambda x: alignment_scores_gap_t[x])

            Sn[i][j] = alignment_scores_no_gap[Pn[i][j]]
            Ss[i][j] = alignment_scores_gap_s[Ps[i][j]]
            St[i][j] = alignment_scores_gap_t[Pt[i][j]]

    max_score_s_t = max([Sn[-1][-1], Ss[-1][-1], St[-1][-1]])
    print(max_score_s_t)

    # add trace step

    s_reconstruct = ''
    t_reconstruct= ''
    i = len(s)
    j = len(t)

    curr_matrix_index = max([0,1,2], key = lambda x: lst_match_matrices[x][i][j])
    prev_matrix_index = lst_prev_matrices[curr_matrix_index][i][j]

    while i > 0 and j > 0:
        if curr_matrix_index == 0:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1
        elif curr_matrix_index == 1:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1
        else:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i-=1
        curr_matrix_index = prev_matrix_index
        prev_matrix_index = lst_prev_matrices[curr_matrix_index][i][j]

    if i == 0 and j == 0:
        x = "nothing" # do nothing
    elif i == 0:
        s_reconstruct = '-'*j + s_reconstruct
        t_reconstruct = t[:j] + t_reconstruct
    else: # j == 0
        s_reconstruct = s[:i] + s_reconstruct
        t_reconstruct = '-'*i + t_reconstruct

    print(s_reconstruct)
    print(t_reconstruct)


# http://rosalind.info/problems/gcon/
def global_alignment_score_constant(s, t, scoring_matrix = BLOSUM62, gap_penalty = 5):
    '''
    constant gap penalty of -5 per contiguous stretch of gaps
    constant gap penalty = affine gap penalty with start penalty = constant and extend penalty = 0
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
    '''
    global_alignment_affine(s, t, scoring_matrix, start_penalty=gap_penalty, extend_penalty=0)


# http://rosalind.info/problems/loca/
def local_alignment_linear(s, t, scoring_matrix = PAM250, gap_penalty = 5):
    '''
    linear gap penalty
    Smith–Waterman algorithm
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/local.pdf
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]

    # we want to track the best score and the location where that best score is taken
    best = 0
    opt_indices = (0, 0)

    # the first row and first column have entry 0

    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            max_score[i][j] = max([max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]],
                                              max_score[i-1][j] - gap_penalty,
                                              max_score[i][j-1] - gap_penalty,
                                             0])
            if max_score[i][j] > best:
                best = max_score[i][j]
                opt_indices = (i, j)

    print(best)
    # trace back:
    i, j = opt_indices
    s_reconstruct = ''
    t_reconstruct = ''

    while max_score[i][j] > 0:

        if max_score[i][j] == max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]]:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i -= 1
            j -= 1

        elif max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j -= 1

        else:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i -= 1

    print(s_reconstruct.replace('-', ''))
    print(t_reconstruct.replace('-', ''))
    print(s_reconstruct)
    print(t_reconstruct)


# http://rosalind.info/problems/laff/

def local_alignment_affine(s, t, scoring_matrix = BLOSUM62, start_penalty = 11, extend_penalty = 1):
    '''
    Affine gap penalty
    Smith–Waterman algorithm
    reference: http://cseweb.ucsd.edu/classes/fa09/cse182/slides/L3.pptx.pdf
    This function is slow for proteins above 1000 aa lengths
    Converted this function to Laff.java
    '''
    Sn = np.zeros((len(t)+1)) # no gap
    Ss = np.zeros((len(t)+1)) # gap at end of s
    St = np.zeros((len(t)+1)) # gap at end of t

    backtrack = 3 * np.ones((len(s)+1, len(t)+1))

    # we no longer need to keep prev_matrix for all 3 matrices
    # since gap_at_end_of_s and gap_at_end_of_t are guaranteed no beter than no_gap_at_end for local alignment

    best, opt_indices = 0, (0,0)

    for i in range(1, len(s)+1):
        new_Sn = np.zeros((len(t)+1))
        new_Ss = np.zeros((len(t)+1))
        new_St = np.zeros((len(t)+1))

        # first entries are all 0

        for j in range(1, len(t)+1):

            curr_score = scoring_matrix.loc[s[i-1], t[j-1]]

            # max_score[0][i][j-1] is guaranteed to be better than max_score[2][i][j-1]
            # max_score[0][i-1][j] is guaranteed to be better than max_score[1][i-1][j]

            new_Ss[j] = max([new_Sn[j-1] - start_penalty,
                             new_Ss[j-1] - extend_penalty])
            new_St[j] = max([Sn[j] - start_penalty,
                             St[j] - extend_penalty])
            candidates = [Sn[j-1] + curr_score,
                          new_Ss[j],
                          new_St[j],
                          0]

            new_Sn[j] = max(candidates)
            backtrack[i,j] = candidates.index(new_Sn[j])

            if new_Sn[j] > best:
                best, opt_indices = new_Sn[j], (i, j)

        Sn = new_Sn
        Ss = new_Ss
        St = new_St

    print(int(best))

    # trace back:

    max_i, max_j = opt_indices
    i = max_i
    j = max_j

    while backtrack[i,j] != 3:

        if backtrack[i,j] == 0:
            i -= 1
            j -= 1

        elif backtrack[i,j] == 1:
            j -= 1

        else:
            i -= 1

    print(s[i:max_i])
    print(t[j:max_j])


# http://rosalind.info/problems/smgb/
def semiglobal_alignment_linear_dna(s, t, match_score = 1, substitution_penalty = 1, gap_penalty = 1):
    '''
    1. same initialization as in local alignment to account for starting gap
    2. same iteration as in global alignment
    3. best = max(last row and last column of score matrix) to account for ending gap
    reference: http://www.cs.cmu.edu/~durand/03-711/2015/Lectures/PW_sequence_alignment_2015.pdf
    4. write a Smgb.java equivalent to deal with large (1000-10000 bp) dnas
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty
            max_score[i][j] = max([max_score[i-1][j-1] + one_match,
                                   max_score[i-1][j] - gap_penalty,
                                   max_score[i][j-1] - gap_penalty])

    best = max([max([max_score[i][-1] for i in range(len(s)+1)]),
                max(max_score[-1])])
    print(best)

    # traceback:

    s_reconstruct = ""
    t_reconstruct = ""

    if best == max([max_score[i][-1] for i in range(len(s)+1)]):
        # gap at end of t
        max_i = [max_score[i][-1] for i in range(len(s)+1)].index(best)
        s_reconstruct = s[max_i:]
        t_reconstruct = '-'*(len(s)-max_i)

        i = max_i
        j = len(t)

    else: # best == max(max_score[-1])
        max_j = max_score[-1].index(best)
        t_reconstruct = t[max_j:]
        s_reconstruct = '-'*(len(t)-max_j)

        j = max_j
        i = len(s)

    while i*j > 0:
        one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1

        elif max_score[i][j] == max_score[i-1][j] - gap_penalty:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = "-" + t_reconstruct
            i-=1

        else: # max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = "-" + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1

    if i == 0 and j == 0: # no gap at start
        x = "nothing" # do nothing

    elif i == 0: # gap at start of s
        s_reconstruct = '-'*j + s_reconstruct
        t_reconstruct = t[:j] + t_reconstruct

    else: # j == 0, gap at start of t
        t_reconstruct = '-'*i + t_reconstruct
        s_reconstruct = s[:i] + s_reconstruct

    print(s_reconstruct)
    print(t_reconstruct)


# http://rosalind.info/problems/oap/
def overlap_alignment_linear_dna(s, t, match_score = 1, substitution_penalty = 2, gap_penalty = 2):
    '''
    1. same initialization as in local alignment to account for trivial suffix of s ('') match with ''
    and '' match with trivial prefix of t ('')
    2. same iteration as in global alignment, so that makes sure always prefix of t match with s
    3. best = max(last row of score matrix) to make sure end at i == -1, so that always suffix of s match with t
    4. traceback until j = 0
    5. write a Oap.java equivalent to deal with large (1000-10000 bp) dnas
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty
            max_score[i][j] = max([max_score[i-1][j-1] + one_match,
                                   max_score[i-1][j] - gap_penalty,
                                   max_score[i][j-1] - gap_penalty])

    best = max(max_score[-1])

    print(best)

    # traceback:

    s_reconstruct = ""
    t_reconstruct = ""

    max_j = max_score[-1].index(best)
    j = max_j
    i = len(s)

    while j > 0:
        one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1

        elif max_score[i][j] == max_score[i-1][j] - gap_penalty:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = "-" + t_reconstruct
            i-=1

        else: # max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = "-" + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1

    print(s_reconstruct)
    print(t_reconstruct)


# http://rosalind.info/problems/sims/
def find_modified_motif_linear_dna(s, t, match_score = 1, substitution_penalty = 1, gap_penalty = 1):
    '''
    local alignment for s, global alignment for t
    traceback until j = 0
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]

    for i in range(0, len(s)+1):
        max_score[i][0] = 0 # local alignment for s

    for j in range(0, len(t)+1):
        max_score[0][j] = -j*gap_penalty  # global alignment for t

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty
            max_score[i][j] = max([max_score[i-1][j-1] + one_match,
                                   max_score[i-1][j] - gap_penalty,
                                   max_score[i][j-1] - gap_penalty])

    candidates = [max_score[i][-1] for i in range(len(s)+1)]
    best = max(candidates)

    print(best)

    # traceback:

    s_reconstruct = ""
    t_reconstruct = ""

    max_i = candidates.index(best)
    i = max_i
    j = len(t)

    while j > 0:
        one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1

        elif max_score[i][j] == max_score[i-1][j] - gap_penalty:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = "-" + t_reconstruct
            i-=1

        else: # max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = "-" + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1

    print(s_reconstruct)
    print(t_reconstruct)


# two slightly different approaches for: http://rosalind.info/problems/ctea/

def count_optimal_alignments_memo(s, t):
    '''
    recursion + memoization approach
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    for i in range(0, len(s)+1):
        max_score[i][0] = -i
    for j in range(0, len(t)+1):
        max_score[0][j] = -j
    index_with_tied_matches = {}

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = 0 if s[i-1] == t[j-1] else -1
            '''
            remember edit distance scoring scheme is that matching symbol does not contribute to score,
            but substitutions and gap contributes -1
            '''
            candidates = [max_score[i-1][j-1] + one_match,
                          max_score[i-1][j] - 1,
                          max_score[i][j-1] - 1]
            max_score[i][j] = max(candidates)

            is_max = [cand == max_score[i][j] for cand in candidates]
            candidate_indices = [(i-1, j-1), (i-1, j), (i, j-1)]
            index_with_tied_matches[(i,j)] = [candidate_indices[k] for k in range(3) if is_max[k]]

    count_matrix = [[0]*(len(t)+1) for i in range(len(s)+1)]

    modulus = 2**27 - 1

    def count(indices):
        i, j = indices
        if i == 0 or j == 0:
            count_matrix[i][j] = 1
            return 1
        elif count_matrix[i][j] != 0:
            return count_matrix[i][j]
        else:
            prev_indices = index_with_tied_matches[(i,j)]
            count_matrix[i][j] = sum([count(prev) for prev in prev_indices]) % modulus
            return count_matrix[i][j]

    return count((len(s),len(t)))


def count_optimal_alignments_dp(s, t):
    '''
    dynamic programming approach
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    count_matrix = [[0]*(len(t)+1) for i in range(len(s)+1)]

    modulus = 2**27 - 1

    for i in range(0, len(s)+1):
        max_score[i][0] = -i
        count_matrix[i][0] = 1
    for j in range(0, len(t)+1):
        max_score[0][j] = -j
        count_matrix[0][j] = 1

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = 0 if s[i-1] == t[j-1] else -1
            '''
            remember edit distance scoring scheme is that matching symbol does not contribute to score,
            but substitutions and gap contributes -1
            '''
            candidates = [max_score[i-1][j-1] + one_match,
                          max_score[i-1][j] - 1,
                          max_score[i][j-1] - 1]
            max_score[i][j] = max(candidates)

            is_max = [cand == max_score[i][j] for cand in candidates]
            candidate_indices = [(i-1, j-1), (i-1, j), (i, j-1)]

            for cand_indices in [candidate_indices[k] for k in range(3) if is_max[k]]:
                i0, j0 = cand_indices
                count_matrix[i][j] += count_matrix[i0][j0]
            count_matrix[i][j] = count_matrix[i][j] % modulus

    return count_matrix[len(s)][len(t)]


# http://rosalind.info/problems/mult/
def generalized_needleman_wunsch(dnas):
    '''
    Generalized Needleman-Wunsch algorithm for n dimension alignment
    Iteration rule: we consider all possible prev_index (see http://readiab.org/book/0.1.2/2/3)
    score_to_add for each prev_index = sum of all possible pairs of -int(s[i-1] != t[j-1]) - number of gaps * non_gaps
    '''
    dims = [len(dna)+1 for dna in dnas]
    n = len(dnas)
    indices = list(range(n))
    S = np.zeros(tuple(dims))
    num_entries = 1
    for dim in dims:
        num_entries *= dim
    traceback = np.array([" "*n for i in range(num_entries)]).reshape(tuple(dims))
    # when initlaizing strings, make the length of string equal to the length of string to be put during iteration

    # initialization, we actually only need to initialize S[0,0,...,0] = 0
    '''
    for k in range(n):
        loc = [0] * n
        for i in range(dims[k]):
            loc[k] = i
            S[tuple(loc)] = -(n-1)*i

    '''

    # iteration:

    # generate all combinations of indices
    # (list(...) because NumPy >= 2.0 returns a tuple, which does not support the item assignment below)
    all_possible_indices = list(np.meshgrid(*[list(range(0, dims[k])) for k in range(n)]))
    for k in range(n):
        while len(all_possible_indices[k].shape) > 1:
            all_possible_indices[k] = np.hstack(all_possible_indices[k]) # flatten the indices
    all_possible_indices = [[all_possible_indices[k][i] for k in range(n)] for i in range(len(all_possible_indices[0]))]

    '''
    Note that in the 2-dimensional case, we only care about (i,0) and (0,j) during initialization
    But in higher dimensions, we have multiple borderline hyperplanes to consider:
    i.e. for 4-dimensional case: we have to consider hyperplanes like (i,j,0,0)
    '''

    for loc in all_possible_indices[1:]: # skip (0,0,0,0)
        candidates = {}
        for k in range(n, 0, -1):
            # we get nCk combinations of choosing index to minus 1
            possible_combinations = combinations(indices, k)

            for comb in possible_combinations:

                new_comb = [i for i in range(n) if (loc[i] > 0 and i in comb)]
                # i.e. if we have comb = (0,1,2) and loc = [0,0,1,1], we can only have prev index as [0,0,0,1]
                if len(new_comb) == 0:
                    continue # i.e. if we have comb = (0,1,2) and loc = [0,0,0,1], we can't minus the first 3 entries
                num = len(new_comb)
                #print(loc, comb, new_comb)
                cand = [loc[i]-1*int(i in new_comb) for i in range(n)]
                score_to_add = 0

                for i in new_comb:
                    for j in new_comb:
                        if i < j:
                            score_to_add -= int(dnas[i][cand[i]] != dnas[j][cand[j]])
                            # add all the scores between different pairs
                score_to_add -= (n-num)*num # gap penalty: (n-num) gaps match with num non-gaps

                #print(loc, cand, score_to_add)
                #print(S[tuple(cand)])
                candidates[tuple(cand)] = S[tuple(cand)] + score_to_add

        #print(loc, candidates)
        traceback[tuple(loc)] = "".join([str(s) for s in list(max(candidates.keys(), key = lambda x: candidates[x]))])
        # store the prev_index as string
        S[tuple(loc)] = candidates[max(candidates.keys(), key = lambda x: candidates[x])]

    #print(S)
    score = S[tuple(all_possible_indices[-1])]
    print(int(score))

    loc = all_possible_indices[-1]
    reconstructs = ["" for k in range(n)]

    # follow the pointers all the way back to the origin: if two or more sequences still have a prefix left,
    # their best alignment can contain gaps, so we can't stop at a border cell and just pad
    while any([index > 0 for index in loc]):
        prev_loc = [int(x) for x in traceback[tuple(loc)]] # unpack indices stored as string
        for k in range(n):
            if prev_loc[k] == loc[k] - 1:
                reconstructs[k] = dnas[k][loc[k] - 1] + reconstructs[k]
            else:
                reconstructs[k] = '-' + reconstructs[k]
        loc = prev_loc

    max_score = 0
    for i in range(n):
        for j in range(i+1, n):
            max_score -= sum([reconstructs[i][k] != reconstructs[j][k] for k in range(len(reconstructs[0]))])

    assert max_score == score # sanity check

    for reconstruct in reconstructs:
        print(reconstruct)


# http://rosalind.info/problems/mgap/
def max_gap(s,t):

    '''
    we set the mismatch penalty to be -inf
    the penalty by adding gap to s and adding gap to t is the same,
    so the max_score must lead to the same gap number under optimized score case
    '''

    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    traceback = [[0]*(len(t)+1) for i in range(len(s)+1)]
    for i in range(0, len(s)+1):
        max_score[i][0] = -i
    for j in range(0, len(t)+1):
        max_score[0][j] = -j

    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = 0 if s[i-1] == t[j-1] else -float('inf')
            '''
            remember edit distance scoring scheme is that matching symbol does not contribute to score,
            but substitutions and gap contributes -1
            '''
            candidates = [max_score[i-1][j-1] + one_match,
                          max_score[i-1][j] - 1,
                          max_score[i][j-1] - 1]
            max_score[i][j] = max(candidates)
            traceback[i][j] = candidates.index(max_score[i][j])

    i = len(s)
    j = len(t)
    gap_num = 0

    while i*j > 0:

        one_match = 0 if s[i-1] == t[j-1] else -float('inf')

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            i-=1
            j-=1
        elif max_score[i][j] == max_score[i-1][j] - 1:
            i-=1
            gap_num += 1
        else:
            j-=1
            gap_num += 1

    if i == 0 and j == 0:
        return gap_num
    elif i == 0:
        return gap_num + j
    else:
        return gap_num + i


# http://rosalind.info/problems/osym/
def align_symbols_cost(s, t):

    def score_matrix(s, t):
        max_score = np.zeros((len(s)+1, len(t)+1))

        for i in range(0, len(s)+1):
            max_score[i,0] = -i # local alignment for s

        for j in range(0, len(t)+1):
            max_score[0,j] = -j # global alignment for t

        for i in range(1,len(s)+1):
            for j in range(1, len(t)+1):
                one_match = 1 if s[i-1] == t[j-1] else -1
                max_score[i][j] = max([max_score[i-1,j-1] + one_match,
                                       max_score[i-1,j] - 1,
                                       max_score[i,j-1] - 1])
        return max_score

    prefix_score_matrix = score_matrix(s, t)
    suffix_score_matrix = score_matrix(s[::-1], t[::-1])

    print(int(prefix_score_matrix[-1][-1]))
    total = 0

    for i in range(len(s)):
        for j in range(len(t)):

            score = prefix_score_matrix[i,j] + int(s[i] == t[j]) - int(s[i] != t[j]) + suffix_score_matrix[len(s)-i-1,len(t)-j-1]
            total += int(score)

    print(total)


# http://rosalind.info/problems/scsp/
# keeps a single table of strings; used by shortest_common_supersequence below
def lcs_string_table(s, t):
    subseq_matrix = [["" for i in range(len(t) + 1)] for j in range(len(s)+1)]
    # sometimes, it is useful the print out the matrix to see where went wrong
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            if s[i-1] == t[j-1]:
                subseq_matrix[i][j] = subseq_matrix[i-1][j-1] + s[i-1]
            else:
                subseq_matrix[i][j] = max(subseq_matrix[i][j-1], subseq_matrix[i-1][j], key = len)
    #print(subseq_matrix)
    return subseq_matrix[len(s)][len(t)]


def shortest_common_supersequence(s, t):
    lcs = lcs_string_table(s, t)
    scs = ""

    while len(lcs) > 0:
        if s[0] == lcs[0] and t[0] == lcs[0]:
            scs += s[0]
            lcs = lcs[1:]
            s = s[1:]
            t = t[1:]
        elif s[0] == lcs[0]: # pause s and add t until s[0] == t[0] == lcs[0]
            scs += t[0]
            t = t[1:]
        else:  # pause t and add s until s[0] == lcs[0]
            scs += s[0]
            s = s[1:]

    return scs + s + t


if __name__ == "__main__":
    lcs_with_length_table('AACCTTGG', 'ACACTGTGA')

    edit_distance('PLEASANTLY', 'MEANLY')

    edit_reconstruct('PRETTY', 'PRTTEIN')

    global_alignment_score_linear('PLEASANTLY', 'MEANLY')

    proteins = read_fasta('data/inputs/rosalind_glob.txt', 'list')

    x = global_alignment_score_linear(proteins[0], proteins[1])

    global_alignment_score_constant('PLEASANTLY', 'MEANLY')

    proteins = read_fasta('data/inputs/rosalind_gcon.txt', 'list')

    global_alignment_score_constant(proteins[0], proteins[1])

    local_alignment_linear('MEANLYPRTEINSTRING', 'PLEASANTLYEINSTEIN')

    local_alignment_affine('PLEASANTLY', 'MEANLY')

    semiglobal_alignment_linear_dna('CAGCACTTGGATTCTCGG', 'CAGCGTGG')

    overlap_alignment_linear_dna('CTAAGGGATTCCGGTAATTAGACAG', 'ATAGACCATATGTCAGTGACTGTGTAA')

    find_modified_motif_linear_dna('GCAAACCATAAGCCCTACGTGCCGCCTGTTTAAACTCGCGAACTGAATCTTCTGCTTCACGGTGAAAGTACCACAATGGTATCACACCCCAAGGAAAC',
                                  'GCCGTCAGGCTGGTGTCCG')

    count_optimal_alignments_memo('PLEASANTLY', 'MEANLY')

    dnas = read_fasta('data/inputs/rosalind_ctea.txt', 'list')

    count_optimal_alignments_dp(dnas[0], dnas[1])

    dnas = read_fasta('data/inputs/rosalind_mult.txt', 'list')

    generalized_needleman_wunsch(dnas)

    max_gap('AACGTA', 'ACACCTA')

    align_symbols_cost('ATAGATA', 'ACAGGTA')

    lcs_string_table('AGGTAB', 'GXTXAYB')

    shortest_common_supersequence('ATCTGAT', 'TGCATA')
