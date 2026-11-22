"""Genome-assembly problems (overlap graphs, de Bruijn graphs, shortest superstring).

Rosalind: grph, long, corr, dbru, pcov, gasm, grep, asmq.

Run from the repository root; `python python/assembly.py` executes the demo calls at the bottom of this file.
"""

import copy

import numpy as np

from utils import read_fasta
from strings import hamming_dist, reverse_complement


# http://rosalind.info/problems/grph/
def overlap_graph(file):
    '''
    For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph Ok
    in which each string is represented by a node, and string s is connected to string t with a directed edge
    when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t
    to prevent directed loops in the overlap graph (although directed cycles may be present).

    Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
    Return: The adjacency list corresponding to O3. You may return edges in any order.
    '''
    dnas = read_fasta(file)
    overlap_dnas = set()
    for dna_id1 in dnas:
        for dna_id2 in dnas:
            if dna_id1 != dna_id2:
                if dnas[dna_id1][-3:] == dnas[dna_id2][:3]:
                    overlap_dnas.add(dna_id1 + ' ' + dna_id2)
                elif dnas[dna_id2][-3:] == dnas[dna_id1][:3]:
                    overlap_dnas.add(dna_id2 + ' ' + dna_id1)
    for item in overlap_dnas:
        print(item)


# https://towardsdatascience.com/solving-tsp-using-dynamic-programming-2c77da86610d

def overlap_dist(dna1, dna2):
    candidate_start_index = [i for i,val in enumerate(dna1) if val==dna2[0]]
    for i in candidate_start_index:
        j1 = i+1
        j2 = 1
        while j1 < len(dna1):
            if dna1[j1] == dna2[j2]:
                j1+=1
                j2+=1
            else:
                break
        if j1 == len(dna1):
            return j2
    return 1


def distance_array(file, output = 'list'):
    dnas = read_fasta(file, 'list')
    if output == 'numpy':
        dist = np.zeros((len(dnas), len(dnas)))
        for i in range(len(dnas)):
            for j in range(len(dnas)):
                dist[i,j] = overlap_dist(dnas[i], dnas[j])
    else:
        dist = [[0]*len(dnas) for i in range(len(dnas))]
        for i in range(len(dnas)):
            for j in range(len(dnas)):
                dist[i][j] = overlap_dist(dnas[i], dnas[j])
    return (dist, dnas)


def shortest_common_superstring(file): # needs an approximate solution, complexity: https://math.stackexchange.com/questions/190313/time-complexity-of-the-travelling-salesman-problem
    dist, dnas = distance_array(file)
    memo = {(tuple([i]), i): tuple([None, 0]) for i in range(len(dnas))}
    # key: prev_visited_indices, prev_last_index; value: next_to_last_index, total dist thus far
    queue = [(tuple([i]), i) for i in range(len(dnas))]
    all_string_indices = set(range(len(dnas)))

    while queue:
        prev_visited_indices, prev_last_index = queue.pop(0)
        _, prev_dist = memo[(prev_visited_indices, prev_last_index)]
        to_visit = all_string_indices.difference(set(prev_visited_indices))

        for next_last_index in to_visit:
            new_visited_indices = tuple(sorted(list(prev_visited_indices)+[next_last_index]))
            new_dist = prev_dist + dist[prev_last_index][next_last_index]

            if (new_visited_indices, next_last_index) not in memo:
                memo[(new_visited_indices, next_last_index)] = (prev_last_index, new_dist)
                queue.append((new_visited_indices, next_last_index))
            elif new_dist > memo[(new_visited_indices, next_last_index)][1]:
                # we want maximum overlap and thus greatest distance between each string
                memo[(new_visited_indices, next_last_index)] = (prev_last_index, new_dist)
                # we do not need to append this one to queue, because as in BFS,
                # we always trace all paths to a smaller subset before heading toward a larger subset
                # so (new_visited_indices, next_last_index) should still be in the queue

    shortest_superstring = retrace_optimal_path(memo, dnas)
    return shortest_superstring


def retrace_optimal_path(memo, dnas):
    superstring = ""
    points_to_retrace = tuple(range(len(dnas)))
    full_path_memo = dict((k,v) for k,v in memo.items() if k[0] == points_to_retrace)
    # only includes keys with full paths
    optimal_choice = max(full_path_memo.keys(), key = lambda x: full_path_memo[x][1])
    last_index = optimal_choice[1]
    superstring += dnas[last_index]
    next_to_last_index, curr_dist = memo[optimal_choice]
    points_to_retrace = tuple(sorted(set(points_to_retrace).difference({last_index})))

    while next_to_last_index is not None:
        last_index = next_to_last_index
        next_to_last_index, prev_dist = memo[(points_to_retrace, last_index)]
        overlap = curr_dist - prev_dist
        superstring = dnas[last_index][:-overlap] + superstring
        curr_dist = prev_dist
        points_to_retrace = tuple(sorted(set(points_to_retrace).difference({last_index})))

    return superstring


# http://rosalind.info/problems/long/ :The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

# create a graph
# topologically sort the graph
# hamiltonian path

def create_graph(file):
    dist, dnas = distance_array(file)
    adjacency_list = {i: [] for i in range(len(dnas))}
    # key: vertex; value: each vertex where there is a directed edge from key to that vertex
    reverse_adjacency_list = {i: [] for i in range(len(dnas))}
    # key: vertex; value: each vertex where there is a directed edge from that vertex to key

    for i in range(len(dnas)):
        for j in range(len(dnas)):
            if i != j and dist[i][j] >= max(len(dnas[i]), len(dnas[j])) / 2:
                adjacency_list[i].append(j)
                reverse_adjacency_list[j].append(i)
    return (adjacency_list, reverse_adjacency_list, dist, dnas)


def topological_sort(adjacency_list, reverse_adjacency_list):

    visited = []
    order = []
    indegree_0_vertices = [v for v in reverse_adjacency_list.keys() if not reverse_adjacency_list[v]]

    def dfs(v):
        visited.append(v)
        for v2 in adjacency_list[v]:
            if v2 not in visited:
                dfs(v2)
        order.append(v)

    for v in indegree_0_vertices:
        dfs(v)

    order = order[::-1]
    return order


def scs(file):

    adjacency_list, reverse_adjacency_list, dist, dnas = create_graph(file)
    order = topological_sort(adjacency_list, reverse_adjacency_list)
    prev_vertex = order[0]
    superstring = dnas[prev_vertex]

    i = 1
    while i < len(dnas):
        next_vertex = order[i]
        overlap = dist[prev_vertex][next_vertex]
        superstring += dnas[next_vertex][overlap:]
        prev_vertex = next_vertex
        i += 1

    return superstring


# http://rosalind.info/problems/corr/
def correct_reads(file):
    dnas = read_fasta(file, 'list')
    dna_counts = {}
    reads = {'correct': [], 'incorrect': []}
    changes = []
    for dna in set(dnas):
        rc_dna = reverse_complement(dna)
        if dna not in dna_counts and rc_dna not in dna_counts:
            dna_counts[dna] = dnas.count(dna) + dnas.count(rc_dna)
            if dna_counts[dna] > 1:
                reads['correct'].append(dna)
            else:
                reads['incorrect'].append(dna)

    for error_dna in reads['incorrect']:
        for correct_dna in reads['correct']:
            if hamming_dist(error_dna, correct_dna) == 1:
                changes.append(error_dna + '->' + correct_dna)
    for change in changes:
        print(change)


# http://rosalind.info/problems/dbru/
def de_bruijn(file):
    with open(file, 'r') as f:
        S = set(f.read().split())
        Src = set([reverse_complement(dna) for dna in S])
    with open('data/outputs/dbru_rslt.txt', 'w') as f:
        for dna in S | Src:
            f.write('(' + dna[:-1] + ', ' + dna[1:] + ')\n')


# http://rosalind.info/problems/pcov/
def construct_perfect_cycle(file):
    with open(file, 'r') as f:
        S = f.read().split()
    adjacency_list = {}
    for dna in S:
        adjacency_list[dna[:-1]] = dna[1:]

    dna = S[0][:-1]
    seq = dna
    while len(seq) < len(adjacency_list):
        dna = adjacency_list[dna]
        seq += dna[-1]
    return seq


# http://rosalind.info/problems/gasm/
def kmers(reads, k):
    adjacency_list = {}
    for read in set(reads) | set([reverse_complement(dna) for dna in reads]):
        i = 0
        while i < len(read) - k:
            adjacency_list[read[i:i+k]] = read[i+1:i+k+1]
            i+=1
    return adjacency_list


def construct_from_reads(file):
    with open(file, 'r') as f:
        reads = f.read().split()

    for k in range(len(reads[0]) - 1, 1, -1):

        adjacency_list = kmers(reads, k)
        kmer = list(adjacency_list.keys())[0]
        seq = kmer
        stop = False

        while seq == kmer or seq[-k:] != seq[:k]:
            if kmer in adjacency_list:
                kmer = adjacency_list[kmer]
                seq += kmer[-1]
            else:
                stop = True
                break

        if stop: # if the cycle does not complete, we need to find the next k
            continue
        else:
            return seq[:-k]


# http://rosalind.info/problems/grep/
def kmers_repeat(file):
    '''
    reads are k + 1 mers
    '''
    with open(file, 'r') as f:
        reads = f.read().split()
    adjacency_list = {}
    for read in reads:
        if read[:-1] in adjacency_list:
            adjacency_list[read[:-1]].append(read[1:])
        else:
            adjacency_list[read[:-1]] = [read[1:]]
    return (adjacency_list, reads[0])
    # we also return read[0] because we want all assembled sequence to start with the first (k+1) mer


def construct_all_seqs(curr_adjacency_lists, curr_kmers, curr_seqs):

    if len(curr_adjacency_lists) == 0 or len(curr_adjacency_lists[0]) == 0:
        # the first case should not happen
        # the second case is when all edges have been traversed
        return curr_seqs

    new_curr_seqs = []
    new_curr_adjacency_lists = []
    new_curr_kmers = []

    for i in range(len(curr_kmers)):

        if curr_kmers[i] in curr_adjacency_lists[i]:

            nxt_kmers = curr_adjacency_lists[i][curr_kmers[i]]

            for nxt_kmer in nxt_kmers:
                if curr_seqs[i]+nxt_kmer[-1] in new_curr_seqs:
                    # we do not want duplicate seqs in our lists
                    continue

                new_curr_kmers.append(nxt_kmer)
                new_curr_seqs.append(curr_seqs[i]+nxt_kmer[-1])

                adjacency_list_copy = copy.deepcopy(curr_adjacency_lists[i])
                adjacency_list_copy[curr_kmers[i]].remove(nxt_kmer)

                if adjacency_list_copy[curr_kmers[i]] == []:
                    adjacency_list_copy.pop(curr_kmers[i]) # for simplicity remove kmer with no out edges

                new_curr_adjacency_lists.append(adjacency_list_copy)
        else:
            continue

    return construct_all_seqs(new_curr_adjacency_lists, new_curr_kmers, new_curr_seqs)


def output_seqs(file):
    adjacency_list, read = kmers_repeat(file)
    adjacency_list[read[:-1]].remove(read[1:])
    if adjacency_list[read[:-1]] == []:
         adjacency_list.pop(read[:-1])
    seqs = construct_all_seqs([adjacency_list], [read[1:]], [read])
    #  we want all assembled sequence to start with the first (k+1) mer

    with open('data/outputs/grep_rslt.txt', 'w') as f:
        for seq in seqs:
            f.write(seq[:-len(read)+1] + '\n')
    # after all edges have been traversed, there will be an overlap between the first kmer in each seq
    # and the last kmer in each seq, and we need to remove that

    return [seq[:-len(read)+1] for seq in seqs]


# http://rosalind.info/problems/asmq/
def nstat(file, prop):
    with open(file, 'r') as f:
        reads = f.read().split()
    lengths = [len(read) for read in reads]
    lengths.sort(reverse = True)
    i = 0
    sums = 0
    while sums < sum(lengths) * prop:
        sums += lengths[i]
        i+=1
    return lengths[i-1]


def output_nstats(file):
    print(nstat(file, 0.5), nstat(file, 0.75))


if __name__ == "__main__":
    shortest_common_superstring('data/inputs/rosalind_long-2.txt') == scs('data/inputs/rosalind_long-2.txt')

    correct_reads('data/inputs/rosalind_corr.txt')

    de_bruijn('data/inputs/rosalind_dbru.txt')

    construct_perfect_cycle('data/inputs/rosalind_pcov.txt')

    construct_from_reads('data/inputs/rosalind_gasm.txt')

    output_seqs('data/inputs/rosalind_grep.txt')

    output_nstats('data/inputs/rosalind_asmq.txt')
