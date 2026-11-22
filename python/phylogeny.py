"""Phylogeny, tree and set problems.

Rosalind: pdst, nwck, nkew, tree, inod, seto.

Run from the repository root; `python python/phylogeny.py` executes the demo calls at the bottom of this file.
"""

from io import StringIO
import re

from Bio import Phylo
import numpy as np

from utils import read_fasta


# http://rosalind.info/problems/tree/
def find_set_index(lst_of_sets, elem):
    '''
    Given a list of sets, return the index of the set in which the element is in
    '''
    for i in range(len(lst_of_sets)):
        if elem in lst_of_sets[i]:
            return i
    print(lst_of_sets, elem)


def compress(lst_of_sets, index1, index2):
    '''
    Given a list of sets, combine the set at index1 and the set at index2 into one set
    '''
    lst_of_sets[index1] = lst_of_sets[index1].union(lst_of_sets[index2])
    lst_of_sets.pop(index2)
    return lst_of_sets


def connected_components(file):
    '''
    return the number of connected components in a graph
    '''
    with open(file, 'r+') as f:
        texts = f.readlines()
    texts = [re.sub('\n', '', text) for text in texts]
    #print(texts)
    n = int(texts[0])
    connected = [{i} for i in range(1, n+1)]
    for text in texts[1:]:
        t1, t2 = re.findall(r'(\d+)\s(\d+)', text)[0]
        index1 = find_set_index(connected, int(t1))
        index2 = find_set_index(connected, int(t2))
        connected = compress(connected, index1, index2)
    return len(connected)


def edges_to_add(file):
    '''
    The minimum number of edges to be added so that a graph becomes a tree (a connected graph)
    is the number of connected components - 1
    '''
    return connected_components(file) - 1


# an even simpler way is to realize that a mimimum of n-1 edges can turn n vertices into a connected graph
# and once the graph is connected, adding additional edges will create cycles

# http://rosalind.info/problems/inod/
def internal_nodes(n):
    '''
    Return: the number of internal nodes (k) in an unrooted binary tree with n leaves
    Each leaf has degree 1
    All internal nodes in an unrooted binary tree have degree 3
    2*edge number = 3k + n
    k+n total vertices --> k+n-1 total edges
    2*(k+n-1) = 3k+n --> k = n - 2
    '''
    return n-2


# http://rosalind.info/problems/pdst/
def p_dist_matrix(file):

    dnas = read_fasta(file, 'list')
    mat = np.zeros((len(dnas), len(dnas)))
    for i in range(len(dnas)):
        dna_i = dnas[i]
        for j in range(len(dnas)):
            dna_j = dnas[j]
            mat[i, j] = np.mean([dnas[i][k] != dnas[j][k] for k in range(len(dna_i))])
        print(' '.join([str(mat[i, j]) for j in range(len(dnas))]))

    return mat


# http://rosalind.info/problems/nwck/
# http://rosalind.info/problems/nkew/
def get_dist(string, node1, node2):
    tree = Phylo.read(StringIO(string), 'newick')
    path1 = tree.get_path(node1)
    path2 = tree.get_path(node2)
    i = 0
    while i < min(len(path1), len(path2)):
        if path1[i] == path2[i]:
            i+=1
        else:
            break
    # identify most recent ancestor
    # and calculate path lengths from there
    return sum([sum([int(path[k].branch_length) if path[k].branch_length else 1 for k in range(i, len(path))]) for path in [path1, path2]])


def output_dist(file):
    dists = []
    with open(file, 'r') as f:
        lst = f.read().split()
    i = 0
    while i < len(lst):
        s = lst[i]
        n1 = lst[i+1]
        n2 = lst[i+2]
        tree = Phylo.read(StringIO(s), 'newick')
        dists.append(str(get_dist(s, n1, n2)))
        i += 3
    print(' '.join(dists))
    return dists


# http://rosalind.info/problems/seto/

# https://www.geeksforgeeks.org/python-set-operations-union-intersection-difference-symmetric-difference/
def set_operation(n, A, B):
    with open('data/outputs/seto_rslt.txt', 'w') as f:
        f.writelines([str(x)+'\n' for x in [A | B,
                    A & B,
                    A - B,
                    B - A,
                    set(range(1,n+1)) - A,
                    set(range(1,n+1)) - B]])


def output_set_op(file):
    with open(file, 'r') as f:
        lst = [x[:-1] for x in f.readlines()] # '\n' counts as one character, strip '\n'
    n = int(lst[0])
    A = set([int(re.sub(r'[{}]', '', x)) for x in lst[1].split(', ')])
    B = set([int(re.sub(r'[{}]', '', x)) for x in lst[2].split(', ')])
    set_operation(n, A, B)


if __name__ == "__main__":
    p_dist_matrix('data/inputs/rosalind_pdst.txt')

    output_dist('data/inputs/rosalind_nwck.txt')

    output_dist('data/inputs/rosalind_nkew.txt')

    output_set_op('data/inputs/rosalind_seto.txt')
