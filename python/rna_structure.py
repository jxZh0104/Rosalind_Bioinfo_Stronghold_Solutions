"""RNA secondary-structure counting problems.

Rosalind: pmch, mmch, cat, motz, rnas.

Run from the repository root; `python python/rna_structure.py` executes the demo calls at the bottom of this file.
"""

from scipy.special import factorial


# http://rosalind.info/problems/pmch/
def perfect_matching(rna):
    '''
    return total number of perfect matching basepair edges
    '''
    assert rna.count('A') == rna.count('U') and rna.count('G') == rna.count('C')
    return factorial(rna.count('A'), exact = True) * factorial(rna.count('G'), exact = True)
    # for i = 1, 2, ... rna.count('A'):
    # for each Ui, we have n-i+1 ways to choose Aj to pair with Ui


# http://rosalind.info/problems/mmch/
def maximum_matching(rna):

    '''
    return total number of maximum matching basepair edges
    '''
    A = rna.count('A')
    U = rna.count('U')
    G = rna.count('G')
    C = rna.count('C')

    AU_match = factorial(max(A,U), exact=True) // factorial(abs(A-U), exact=True) # essentially a permutation
    GC_match = factorial(max(G,C), exact=True) // factorial(abs(G-C), exact=True)

    return AU_match * GC_match


# http://rosalind.info/problems/cat/
def noncrossing(rna, memo):
    '''
    return total number of noncrossing perfect matching basepair edges
    '''
    if rna not in memo:
        if rna.count('A') != rna.count('U') or rna.count('G') != rna.count('C'):
            memo[rna] = 0
        else:
            memo[rna] = sum([noncrossing(rna[0]+rna[i], memo)*(noncrossing(rna[1:i], memo))*(noncrossing(rna[i+1:], memo)) for i in range(1, len(rna), 2)])
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna]) % 1000000


# http://rosalind.info/problems/motz/
def noncrossing_all(rna, memo):
    '''
    return total number of noncrossing matching basepair (could be non-perfect) edges
    '''
    if rna not in memo:
        if len(rna) <= 1:
            memo[rna] = 1
        else:
            memo[rna] = noncrossing_all(rna[1:], memo) # base 1 is not involved in any base pair
            base_pair = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}
            for i in range(1, len(rna)):
                if base_pair[rna[0]] == rna[i]:
                    memo[rna] += noncrossing_all(rna[1:i], memo) * noncrossing_all(rna[i+1:], memo)
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna]) % 1000000


# http://rosalind.info/problems/rnas/
def noncrossing_with_wobble(rna, memo):
    '''
    return total number of noncrossing matching basepair (could be non-perfect) edges (including U-G basepair)
    two bases can pair only if they are at least 4 bases apart
    '''
    if rna not in memo:
        if len(rna) <= 3:
            memo[rna] = 1
        else:
            memo[rna] = noncrossing_with_wobble(rna[1:], memo) # base 1 is not involved in any base pair
            base_pair = {'A':['U'], 'U':['A', 'G'], 'C':['G'], 'G':['U','C']}
            for i in range(4, len(rna)):
                if rna[i] in base_pair[rna[0]]:
                    memo[rna] += noncrossing_with_wobble(rna[1:i], memo) * noncrossing_with_wobble(rna[i+1:], memo)
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna])


if __name__ == "__main__":
    perfect_matching('AGCUAGUCAU')

    maximum_matching('AUGCUUC')

    total_numbers = {'':1,'A':0, 'G':0, 'U':0, 'C':0, 'AA':0, 'GG':0, 'UU':0, 'CC':0,
                        'AG':0, 'AC':0, 'CA':0, 'CU':0, 'GA':0, 'GU':0, 'UC':0, 'UG':0,
                        'AU':1, 'UA':1, 'CG':1, 'GC':1}
    # we need to add '' to account for the edge that connects the first base to the last base

    noncrossing('AUAU', total_numbers)

    noncrossing_all('AUAU', {})

    noncrossing_with_wobble('GCGGUUCGAGUGCCGUCAUCCGAGGAAUCGUAUGCAUAGUAGUAGGAUUGGACCUGCAAGUAGCAGCGUUUCAACCGUUGUUACCGGAGUAUACUGGAGGGGGUCGCCACUCGGACGUGCAUCCAGUAUUUCAGUUGUUAGCAGUGAUGAUCGUUUGGAUUUAUGGUCUUUUUGUU', {})
