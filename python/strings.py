"""String, counting and motif-finding problems.

Rosalind: dna, rna, revc, fib, fibd, gc, hamm, perm, sign, lexf, kmer, iprb, lia, rstr, prot, mrna, orf,
splc, subs, cons, kmp, sseq, revp, lcsm, mprt.

Run from the repository root; `python python/strings.py` executes the demo calls at the bottom of this file.
"""

import re

import Bio.Data.CodonTable
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from scipy.special import comb

from utils import convert_dict_to_strings, convert_tuples_to_strings, index_transform, read_fasta


# http://rosalind.info/problems/dna/
def count_base(dna):
    dna = dna.upper()
    return [dna.count(b) for b in ['A', 'C', 'G', 'T']]


# http://rosalind.info/problems/rna/
def transcribe(dna):
    dna = dna.upper()
    return dna.replace('T', 'U')


# http://rosalind.info/problems/revc/
def reverse_complement(dna):
    base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    dna = dna.upper()
    return ''.join(map(lambda x: base_pair[x], dna))[::-1]


# http://rosalind.info/problems/fib/
def fibonacci(n, k):
    '''
    Given: Positive integers n≤40 and k≤5
    Return: The total number of rabbit pairs that will be present after n months,
            if we begin with 1 pair and in each generation,
            every pair of reproduction-age rabbits (>= 2 months old) produces
            a litter of k rabbit pairs (instead of only 1 pair).
    '''
    total_number = [1]*n
    for i in range(2, n):
        total_number[i] = total_number[i-1] + k * total_number[i-2]
    return total_number[-1]


# http://rosalind.info/problems/fibd/

def mortal_fibonacci(n, m): # use python list, no integer overflow problem
    '''
    Given: Positive integers n≤100 and m≤20
    Return: The total number of pairs of rabbits that will remain after the nth month
            if all rabbits live for m months. (each rabbit reproduces only when >= 2 month old)
    '''
    total_number = []
    for i in range(n):
        lst = []
        for j in range(m+1):
            lst.append(0)
        total_number.append(lst)
    # the entry at jth column and ith row represents the number of pair of rabbit
    # that is j month old at the i+1 th month
    total_number[0][0] = 1
    total_number[1][1] = 1
    for i in range(2, n):
        for j in range(1, min(i, m)+1):
            total_number[i][j] = total_number[i-1][j-1]
            # update the age of rabbit in each month
        total_number[i][0] = sum(total_number[i][2:(m+1)])
        # print(total_number)
        # only the rabbits >= 2 months old can reproduce
    #print(total_number)
    return int(sum(total_number[n-1][:m]))


# http://rosalind.info/problems/gc/
def gc_content(dna):
    return np.mean([s == 'G' or s == 'C' for s in dna if s in ['A', 'T', 'G', 'C']])*100


def get_highest_gc(file):
    dnas = read_fasta(file)
    GCs = {x: gc_content(dnas[x]) for x in dnas}
    winner = max(GCs)
    print(winner+'\n'+str(round(GCs[winner], 7)))


# http://rosalind.info/problems/hamm/
def hamming_dist(dna1, dna2):
    assert len(dna1)==len(dna2)
    return sum([dna1[i] != dna2[i] for i in range(len(dna1))])


# http://rosalind.info/problems/perm/
def rearrange(n):
    '''
    > rearrange(3)
    [[3, 1, 2], [3, 2, 1], [1, 3, 2], [1, 2, 3], [2, 1, 3], [2, 3, 1]]
    '''
    permutations = {1: [[1]], 2: [[1,2], [2,1]]}

    def switch(item, i, j):
        '''
        switch i to j and j to i and keep all other entries as the same
        '''
        if item == i:
            return j
        elif item == j:
            return i
        else:
            return item

    for i in range(3, n+1):
        total_arrangement = []
        temp_permut_minus_1 = [[i] + x for x in permutations[i-1]]
        total_arrangement.extend(temp_permut_minus_1.copy())
        # e.g. for i = 3, this is [3,1,2] and [3,2,1]
        for j in range(1, i):
            switched_permut = [list(map(lambda item: switch(item, i, j), lst)) for lst in temp_permut_minus_1]
            total_arrangement.extend(switched_permut)
            # switch between i and j e.g. for j = 1, this is [1,3,2] and [1,2,3]
        permutations[i] = total_arrangement
    return permutations[n]


def map_lsts_to_strings(lst):
    '''
    > map_lsts_to_strings([[3, 1, 2], [3, 2, 1], [1, 3, 2], [1, 2, 3], [2, 1, 3], [2, 3, 1]])
    6
    3 1 2
    3 2 1
    1 3 2
    1 2 3
    2 1 3
    2 3 1
    '''
    print(len(lst))
    for item in lst:
        print(' '.join([str(x) for x in item]))


# http://rosalind.info/problems/sign/
def signed_permutation(n):
    '''
    >2
    -1 -2
    -1 2
    1 -2
    1 2
    -2 -1
    -2 1
    2 -1
    2 1
    '''
    def map_to_sign(item, i):
        '''
        convert i to -i in a lst
        '''
        if item == i:
            return -i
        else:
            return item

    unsigned_permutation = rearrange(n)
    for i in range(1, n+1):
        # double the current list by adding the possibility of converting each i to -i for all i = 1, 2, 3, ..., n
        new_permutation = [list(map(lambda item: map_to_sign(item, i), lst)) for lst in unsigned_permutation]
        unsigned_permutation.extend(new_permutation)
    return unsigned_permutation


def map_lst_to_strings_in_file(lst, filename):
    with open(filename, 'w+') as f:
        f.writelines([str(len(lst))+'\n'])
        f.writelines([' '.join([str(x) for x in item])+'\n' for item in lst])


# http://rosalind.info/problems/lexf/
def lexicograph(lst, n):
    '''
    Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n(n≤10).
    Return: All strings of length n that can be formed from the alphabet,
            ordered lexicographically (use the standard order of symbols in the English alphabet).
    '''

    lst = [s.upper() for s in lst]
    lst.sort()
    possible_strings = {1:lst.copy()}
    for i in range(2, n+1):
        temp_strings = []
        for s in lst:
            temp_strings.extend([s + string for string in possible_strings[i-1]])
        possible_strings[i] = temp_strings
    for string in possible_strings[n]:
        print(string)
    return possible_strings[n]


# http://rosalind.info/problems/kmer/
def kmer_composition(dna, k):
    kmer_lst = lexicograph(['A', 'T', 'G', 'C'], k)
    freq = [0]*4**k
    for i in range(4**k):
        freq[i] = len(find_motif(dna, kmer_lst[i], 'list'))
    print(" ".join([str(x) for x in freq]))
    return freq


# http://rosalind.info/problems/iprb/
def mendel_prob(k, m, n):
    '''
    Given: Three positive integers k, m, and n , representing a population containing k+m+n organisms:
           k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
    Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele
            (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    '''
    # two dominant: choose(k, 2) p = 1
    # two heterozygous: choose(m, 2) p = 3/4
    # two recessive: choose(n, 2) p = 0
    # dominant + heterozygous: k*m, p = 1
    # dominant + recessive: k*n, p = 1
    # heterozygous + recessive: m*n, p = 1/2
    # total choices: choose(k+m+n, 2)
    return (k*(k-1)/2 * 1 + m*(m-1)/2 * 3/4 + n*(n-1)/2 * 0 + k*m*1 + k*n*1 + m*n/2)/((k+m+n)*(k+m+n-1)/2)


# http://rosalind.info/problems/prot/
def translate(mrna):
    return str(Seq(mrna).translate())[:-1] #the last one is a stop codon


# http://rosalind.info/problems/mrna/
codon_table = Bio.Data.CodonTable.standard_rna_table.forward_table
reverse_codon_table = {}
for codon in codon_table:
    aa = codon_table[codon]
    if aa not in reverse_codon_table:
        reverse_codon_table[aa] = [codon]
    else:
        reverse_codon_table[aa].append(codon)


def infer_mrna(protein):
    '''
    Given: A protein string of length at most 1000 aa.
    Return: The total number of different RNA strings from which the protein could have been translated,
            modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)
    '''
    total_number = 3 # 3 options for stop codon
    for aa in protein:
        total_number = total_number * len(reverse_codon_table[aa]) % 1000000
    return total_number


# http://rosalind.info/problems/subs/
def find_motif(s, t, output = 'string'):
    '''
    Given: Two DNA strings s and t (each of length at most 1 kbp).
    Return: All locations of t as a substring of s.
    '''
    possible_indices = [str(i+1) for i in range(len(s)-len(t)+1) if s[i:i+len(t)] == t]
    if output == 'list':
        return possible_indices
    return ' '.join(possible_indices)


# http://rosalind.info/problems/cons/
def consensus_string(file):
    dnas = read_fasta(file)
    consensus_dict = {'A':[], 'C':[], 'G':[], 'T':[]}
    for i in range(len(dnas[list(dnas.keys())[0]])):
        for base in consensus_dict:
            consensus_dict[base].append(sum([dna[i] == base for dna in dnas.values()]))
    consensus_matrix = pd.DataFrame(consensus_dict).transpose()
    print(''.join(list(consensus_matrix.idxmax())))
    for base in consensus_dict:
        print(base + ': ' + ' '.join([str(count) for count in consensus_dict[base]]))


# http://rosalind.info/problems/kmp/
# kmf algorithm reference: https://www.inf.hs-flensburg.de/lang/algorithmen/pattern/kmpen.htm
def failure_array(dna, filename = False):
    max_lengths = [0]
    for k in range(1, len(dna)):
        max_lengths.append(0)
        # the sup for prefix = suffix length at position k is reached if we can extend the max length at k-1 by 1
        # i.e. if dna[max_lengths[k-1]] = dna[k]
        if dna[max_lengths[k-1]] == dna[k]:
            max_lengths[k] = max_lengths[k-1] + 1
        else:
            for j in range(max_lengths[k-1], 0, -1):
                if dna[:j] == dna[(k-j+1):(k+1)]:
                    max_lengths[k] = j
                    break
    if filename:
        with open(filename, 'w+') as f:
            f.write(' '.join([str(l) for l in max_lengths]))
    return max_lengths


#reference: https://towardsdatascience.com/pattern-search-with-the-knuth-morris-pratt-kmp-algorithm-8562407dba5b


def find_lps(pattern):

    # name lps indicates longest proper prefix which is also suffix.
    prefix = 0
    lps_array = [0]*len(pattern)

    for i in range(1, len(pattern)):

        while pattern[prefix] != pattern[i] and prefix > 0:
            prefix = lps_array[prefix-1]
            # use proof by contradiction to show that there is no possible match between
            # lps_array[prefix-1] and prefix

        if pattern[prefix] == pattern[i]:
            prefix += 1
            lps_array[i] = prefix

    return lps_array


def kmp(pattern, text, unique):

    # Knuth-Morris-Pratt algorithm

    prefix = 0 # use prefix to keep track of matches
    lps_array = find_lps(pattern)
    matched_indices = []

    for i, ch in enumerate(text):

        while pattern[prefix] != ch and prefix > 0:
            prefix = lps_array[prefix-1]

        if pattern[prefix] == ch:

            if prefix == len(pattern) - 1: # full match
                matched_indices.append(i - prefix)

                if unique: # if we know the match is uqniue, we can return the match as soon as we find one
                    return matched_indices

                prefix = lps_array[prefix] # we turn the suffix into prefix and find the next match

            else:
                prefix += 1

    return matched_indices


# http://rosalind.info/problems/splc/
def exon_protein(file):
    '''
    first dna in file is the entire sequence
    the rest dnas are introns

    Return: a protein made by concatanating the codons
    '''
    dnas = read_fasta(file, 'list')
    text = dnas[0]
    patterns = dnas[1:]
    for pattern in patterns:
        start = kmp(pattern, text, True)[0]
        end = start + len(pattern) - 1
        text = text[0:start] + text[end+1:]
    raw_protein = translate(text)
    return raw_protein.replace('*', '') # remove stop codons


# http://rosalind.info/problems/sseq/
def find_spliced_motif(s, t):
    '''
    Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
    Return: One collection of indices of s in which the symbols of t appear as a subsequence of s.
            If multiple solutions exist, you may return any one.
    '''
    s_index = 0
    t_index = 0
    subsequence = []
    while t_index < len(t):
        if s[s_index] == t[t_index]:
            subsequence.append(s_index)
            s_index += 1
            t_index += 1
        else:
            s_index += 1
    return index_transform(subsequence)


# http://rosalind.info/problems/revp/
def check_palindrome(dna):
    if len(dna) % 2 != 0:
        return False
    forward_pointer = 0
    backward_pointer = len(dna)-1
    base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    while forward_pointer < backward_pointer:
        if dna[forward_pointer] == base_pair[dna[backward_pointer]]:
            forward_pointer += 1
            backward_pointer -= 1
        else:
            return False
    return True


def check_palindrome_brute_force(dna):
    return reverse_complement(dna) == dna


def find_palindrome(dna):

    '''
    Return a list of tuples:
    each tuple: (i, l)
    where i is the position of start base in the dna (index starts from 1)
    and l is the length of reverse palindromes
    if there are mutlple l for oen particular i, return the largest l
    '''

    palindrome_list = []

    for i in range(len(dna)):

        match_indices = []

        candidate_palindromes = [] # candidate palindromes of all lengths

        j = i + 2

        while j < len(dna):

            candidate_j = 0

            pattern = reverse_complement(dna[i:j])
            '''
            for example, if the ith index is 'A', and pattern = 'AG',
            the reverse complement of 'AG' is 'CT'
            we want to see if 'CT' is in dna,
            if not, there is no palindrome that starts from 'A'
            if yes, we find the largest length of palindrome based on all the 'CT' matches
            and increase j until there is no palindrome that starts from 'A'
             '''
            #print(dna[i:j], pattern)
            matched_indices = [x for x in kmp(pattern, dna, False) if x > i]
            #print(matched_indices)

            if matched_indices:
                k = len(matched_indices) - 1
                while k >= 0:
                    candidate = dna[i:(matched_indices[k]+(j-i))]
                    if check_palindrome(candidate):

                        candidate_j = len(candidate)
                        #print('yes', len(candidate))
                        break
                    else:
                        k -= 1

            else:
                break

            if candidate_j > 0:
                candidate_palindromes.append(candidate_j)

            #print(candidate_j)

            j += 1
        #print(candidate_palindromes)

        if candidate_palindromes:
            palindrome_list.append((i+1, max(candidate_palindromes)))
            # i+1 to account for python index start from 0

    return palindrome_list


def find_palindrome_brute_force(dna):
    '''
    Return a list of tuples:
    each tuple: (i, l)
    where i is the position of start base in the dna (index starts from 1)
    and l is the length of reverse palindromes
    there could be multiple l for each i
    '''
    palindrome_candidates = []
    for i in range(len(dna)-3):
        for j in range(i+4, i+13):
            if j > len(dna):
                break
            candidate = dna[i:j]
            if check_palindrome_brute_force(candidate):
                palindrome_candidates.append((i+1, len(candidate)))
    return palindrome_candidates


# http://rosalind.info/problems/lcsm/
def dna_lcs(file):
    from suffix_trees import STree  # optional dependency: imported here so that this module imports without it
    # longest common string
    dnas = read_fasta(file, 'list')
    st = STree.STree(dnas)
    return st.lcs()


# http://rosalind.info/problems/mprt/
def glycosylation_motif(file):
    # optional dependencies: imported here so that this module imports without them
    import requests
    from bs4 import BeautifulSoup
    with open(file, 'r+') as f:
        protein_ids = f.read().split('\n')
    returned = {}
    for protein in protein_ids:
        url = "https://www.uniprot.org/uniprot/" + protein + ".fasta"
        page = requests.get(url)
        soup = BeautifulSoup(page.content, 'html.parser')
        seq = ''.join(str(soup).split('\n')[1:-1])
        matches = re.finditer(r'(?=(N[ARNDCQEGHILKMFSTWYV][ST]))', seq) # r'(?=(pattern))' --> find overlapping matches
        # N-glycosylation motif: Asn-X-Serine/Threonine where X is any amino acid other than Proline
        returned[protein] = []
        for m in matches:
            returned[protein].append(m.start())
    return returned


# http://rosalind.info/problems/orf/
def find_orfs(dna):

    mrna = dna.replace('T', 'U')
    complement_mrna = reverse_complement(dna).replace('T', 'U')
    start_strand1 = list(re.finditer(r'(AUG)', mrna))
    start_strand2 = list(re.finditer(r'(AUG)', complement_mrna))
    # both strands (strand2 is the reverse complement of strand1) can be coding strand --> 6 ORs

    proteins = set()

    for starts, seq in [(start_strand1, mrna), (start_strand2, complement_mrna)]:

        for m in starts:

            i = m.start()
            while seq[i:i+3] not in ['UAG', 'UAA', 'UGA'] and i < len(seq):
                i += 3
            if i < len(seq):
                p = translate(seq[m.start():i+3])
                proteins.add(p)

    for p in proteins:

        print(p)


# http://rosalind.info/problems/lia/
def prob_lia(k, N):
    genotype_probs = [0, 1, 0] # AA, Aa, aa for single allele
    i = 0
    while i <= N:
        print(genotype_probs)
        genotype_probs = [genotype_probs[0] * 0.5 + genotype_probs[1] * 0.25,
                          genotype_probs[0] * 0.5 + genotype_probs[1] * 0.5 + genotype_probs[2] * 0.5,
                          genotype_probs[2] * 0.5 + genotype_probs[1] * 0.25]
        i += 1

    p = genotype_probs[1] ** 2 # prob for AaBb, in fact, this is always 0.25

    if N > 2 ** (k-1):
        return sum([comb(2**k, i) * p**i * (1-p)**(2**k-i) for i in range(N, 2**k + 1)])
    else:
        return 1 - sum([comb(2**k, i) * p**i * (1-p)**(2**k-i) for i in range(N)])


# http://rosalind.info/problems/rstr/
def prob_rstr(N, GC, dna):
    prob = 1
    for b in dna:
        if b in 'GC':
            prob *= GC / 2
        else:
            prob *= (1-GC)/2

    return 1 - (1-prob)**N


if __name__ == "__main__":
    count_base('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC')

    transcribe('GATGGAACTTGACTACGTAAATT')

    reverse_complement('AAAACCCGGT')

    fibonacci(5,3)

    gc_content('CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT')

    get_highest_gc('data/inputs/rosalind_gc.txt')

    map_lsts_to_strings(rearrange(3))

    map_lst_to_strings_in_file(signed_permutation(3), 'data/outputs/test.txt')

    lexicograph(['A','G', 'T', 'C'], 2)

    mendel_prob(2,2,2)

    infer_mrna('MA')

    find_motif('GATATATGCATATACTT',
    'ATAT')

    consensus_string('data/inputs/rosalind_cons.txt')

    failure_array('CAGCATGGTATCACAGCAGAG')  # failure array = lps

    exon_protein('data/inputs/rosalind_splc.txt')

    find_spliced_motif('ACGTACGTGACG', 'GTA')

    find_palindrome('TCAATGCATGCGGGTCTATATGCAT')

    convert_tuples_to_strings(find_palindrome('TCAATGCATGCGGGTCTATATGCAT'))

    convert_dict_to_strings(glycosylation_motif('data/inputs/rosalind_mprt.txt'))

    find_orfs('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')

    prob_lia(7, 34)

    prob_rstr(82895, 0.555535, 'ACCCTCTCGC')
