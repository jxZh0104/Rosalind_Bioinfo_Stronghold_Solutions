"""Hidden Markov model problems (Bioinformatics Textbook Track).

Rosalind: ba10a, ba10b, ba10c, ba10d, ba10e, ba10f, ba10g, ba10h, ba10i, ba10j, ba10k.

Run from the repository root; `python python/hmm.py` executes the demo calls at the bottom of this file.
"""

import re

import numpy as np
import pandas as pd

from utils import write_matrix_to_file


# http://rosalind.info/problems/ba10a/
def markov_prob(file):
    with open(file, 'r') as f:
        lst = f.readlines()
    path = lst[0][:-1]
    mat = np.zeros((2,2))
    mat[0,] = [float(x) for x in re.split(r'\s', lst[-2]) if x not in 'A']
    mat[1,] = [float(x) for x in re.split(r'\s', lst[-1]) if x not in 'B']
    mat = pd.DataFrame(mat, index = ['A', 'B'], columns = ['A', 'B'])
    print(mat)

    prob = 1/2 # initial probs are equal
    i = 0
    while i < len(path) - 1:
        prob *= mat.loc[path[i], path[i+1]]
        i += 1
    return prob


# http://rosalind.info/problems/ba10b/
def markov_prob_2(file):
    with open(file, 'r') as f:
        lst = f.readlines()
    string = lst[0][:-1]
    path = lst[4][:-1]
    mat = np.zeros((2,3))
    mat[0,] = [float(x) for x in re.split(r'\s', lst[-2]) if x not in 'A']
    mat[1,] = [float(x) for x in re.split(r'\s', lst[-1]) if x not in 'B']
    mat = pd.DataFrame(mat, index = ['A', 'B'], columns = ['x', 'y', 'z'])
    print(mat)

    prob = 1
    i = 0
    while i < len(path):
        prob *= mat.loc[path[i], string[i]]
        i += 1
    return prob


# 10c, d, and k reference to https://web.stanford.edu/~jurafsky/slp3/A.pdf and https://people.cs.umass.edu/~mccallum/courses/inlp2004a/lect10-hmm2.pdf

# http://rosalind.info/problems/ba10c/
def precompute_hmm_info(file):

    # create transition probability matrix, emission probability matrix, observation string, and states
    with open(file, 'r') as f:
        lst = f.readlines()
    obs = lst[0][:-1]
    emit_chars = [x for x in re.split(r'\s', lst[2]) if x != '']
    states = [x for x in re.split(r'\s', lst[4]) if x != '']
    trans_mat = np.zeros((len(states), len(states)))
    emit_mat = np.zeros((len(states), len(emit_chars)))
    for i in range(len(states)):
        trans_mat[i,] = [float(x) for x in re.split(r'\s', lst[7+i]) if x not in states[i]]
    trans_mat = pd.DataFrame(trans_mat, index = states, columns = states)
    for i in range(len(states)):
        emit_mat[i,] = [float(x) for x in re.split(r'\s', lst[7+len(states)+2+i]) if x not in states[i]]
    emit_mat = pd.DataFrame(emit_mat, index = states, columns = emit_chars)
    print(trans_mat)
    print(emit_mat)

    return (trans_mat, emit_mat, obs, states)


def viterbi(trans_mat, emit_mat, obs, states, printPath = True):

    #trans_mat, emit_mat, obs, states = precompute_hmm_info(file)

    initial_probs = [1/len(states) for state in states] # assumes all states are equally likely at the start

    path_probs = np.zeros((len(states), len(obs)))
    traceback = np.zeros((len(states), len(obs)))

    # initialization
    for i in range(len(states)):
        path_probs[i, 0] = initial_probs[i] * emit_mat.loc[states[i], obs[0]]
        traceback[i, 0] = 0

    # recursion
    for t in range(1, len(obs)):
        for i in range(len(states)):
            candidates = [path_probs[j, t-1] * trans_mat.iloc[j, i] * emit_mat.loc[states[i], obs[t]] for j in range(len(states))]
            path_probs[i, t] = max(candidates)
            traceback[i, t] = candidates.index(max(candidates))

    # termination
    best_path_prob = max(path_probs[:,-1])

    # traceback
    best_loc = np.argmax(path_probs[:,-1])
    t = len(obs) - 1
    path = [states[best_loc]]

    while t > 0:
        best_loc = int(traceback[best_loc, t])
        path = [states[best_loc]] + path
        t -= 1
    if printPath:
        print(''.join(path))
    return path


# http://rosalind.info/problems/ba10d/
def hmm_likelihood(file):

    '''
    apply forward algorithm
    '''

    trans_mat, emit_mat, obs, states = precompute_hmm_info(file)

    initial_probs = [1/len(states) for state in states] # assumes all states are equally likely at the start

    forward_probs = np.zeros((len(states), len(obs)))

    # initialization
    for i in range(len(states)):
        forward_probs[i, 0] = initial_probs[i] * emit_mat.loc[states[i], obs[0]]

    # recursion
    for t in range(1, len(obs)):
        for j in range(len(states)):
            candidates = [forward_probs[i, t-1] * trans_mat.iloc[i, j] * emit_mat.loc[states[j], obs[t]] for i in range(len(states))]
            forward_probs[j, t] = sum(candidates)

    # termination
    return sum(forward_probs[:,-1])


# http://rosalind.info/problems/ba10j/

def forward_algo(trans_mat, emit_mat, obs, states):

    initial_probs = [1/len(states) for state in states] # assumes all states are equally likely at the start

    forward_probs = np.zeros((len(states), len(obs)))

    # initialization
    for i in range(len(states)):
        forward_probs[i, 0] = initial_probs[i] * emit_mat.loc[states[i], obs[0]]

    # recursion
    for t in range(1, len(obs)):
        for j in range(len(states)):
            candidates = [forward_probs[i, t-1] * trans_mat.iloc[i, j] * emit_mat.loc[states[j], obs[t]] for i in range(len(states))]
            forward_probs[j, t] = sum(candidates)

    return forward_probs


def backward_algo(trans_mat, emit_mat, obs, states):

    backward_probs = np.zeros((len(states), len(obs)))

    # initialization
    for i in range(len(states)):
        backward_probs[i, -1] = 1

    # recursion
    for t in range(len(obs)-2, -1, -1):
        for i in range(len(states)):
            candidates = [backward_probs[j, t+1] * trans_mat.iloc[i, j] * emit_mat.loc[states[j], obs[t+1]] for j in range(len(states))]
            backward_probs[i, t] = sum(candidates)

    return backward_probs


def soft_decoding(file = 'data/inputs/rosalind_ba10j.txt'):
    '''
    we basically want the conditional probability of each state at each time t given the emitted string
    which is equal to the forward prob * backward prob of each state at each time t
    divided by the likelihood of tHe emitted string
    see: https://www.inf.ed.ac.uk/teaching/courses/asr/2008-9/asr-hmm.pdf and search for 'soft'
    '''

    trans_mat, emit_mat, obs, states = precompute_hmm_info(file)

    forward_probs = forward_algo(trans_mat, emit_mat, obs, states)
    backward_probs = backward_algo(trans_mat, emit_mat, obs, states)

    # print(sum(forward_probs[:, -1]), sum(backward_probs[:, 0]*emit_mat.loc[:,obs[0]])/len(states))

    assert np.isclose(sum(forward_probs[:, -1]), sum(backward_probs[:, 0]*emit_mat.loc[:,obs[0]])/len(states))
    # check that the likelihood of observation is the same when computed by forward and backward probs

    likelihood = sum(forward_probs[:, -1])

    cond_probs = pd.DataFrame(np.zeros((len(obs), len(states))), columns = states, index = list(range(len(obs))))

    for i in range(len(states)):
        for t in range(len(obs)):
            cond_probs.iloc[t, i] = forward_probs[i, t] * backward_probs[i, t] / likelihood

    with open('data/outputs/ba10j_rslt.txt', 'w') as f:
        f.write('\t'.join(states)+'\n')
        for t in range(len(obs)):
            f.write('\t'.join([str(round(x, 3)) if x != int(x) else str(x) for x in cond_probs.iloc[t, :]]) + '\n')


# http://rosalind.info/problems/ba10k/
def get_hmm_learning_info(file):
    with open(file, 'r') as f:
        lst = f.readlines()
    iter_lim = int(lst[0][:-1])
    obs = lst[2][:-1]
    emit_chars = [x for x in re.split(r'\s', lst[4]) if x != '']
    states = [x for x in re.split(r'\s', lst[6]) if x != '']
    trans_mat = np.zeros((len(states), len(states)))
    emit_mat = np.zeros((len(states), len(emit_chars)))
    for i in range(len(states)):
        trans_mat[i,:] = [float(x) for x in re.split(r'\s', lst[9+i][1:]) if x != ''] # treat the matrices posiitonally, some row labels are mislabeled
    trans_mat = pd.DataFrame(trans_mat, index = states, columns = states)
    for i in range(len(states)):
        emit_mat[i,:] = [float(x) for x in re.split(r'\s', lst[9+len(states)+2+i][1:]) if x != '']
    emit_mat = pd.DataFrame(emit_mat, index = states, columns = emit_chars)

    return (iter_lim, obs, emit_chars, states, trans_mat, emit_mat)


def baum_welch(file):

    iter_lim, obs, emit_chars, states, trans_mat, emit_mat = get_hmm_learning_info(file)

    num_iter = 0
    print(iter_lim)

    while num_iter < iter_lim:

        forward_probs = forward_algo(trans_mat, emit_mat, obs, states)
        backward_probs = backward_algo(trans_mat, emit_mat, obs, states)

        '''
        assert np.isclose(sum(forward_probs[:, -1]), sum(backward_probs[:, 0]*emit_mat.iloc[:,0])/len(states))
        the above is in fact incorrect
        the correct version is:
        '''
        assert np.isclose(sum(forward_probs[:, -1]), sum(backward_probs[:, 0]*emit_mat.loc[:,obs[0]])/len(states))
        # check that the likelihood of observation is the same when computed by forward and backward probs

        likelihood = sum(forward_probs[:, -1])

        # E-step
        gamma_mat = np.zeros((len(states), len(obs)))
        for t in range(len(obs)):
            for j in range(len(states)):
                gamma_mat[j, t] = forward_probs[j, t] * backward_probs[j, t]/likelihood

        eta_mat = np.zeros((len(states), len(states), len(obs)))
        for t in range(len(obs)-1):
            for i in range(len(states)):
                for j in range(len(states)):
                    eta_mat[i, j, t] = forward_probs[i, t] * trans_mat.iloc[i, j] * emit_mat.loc[states[j], obs[t+1]] * backward_probs[j, t+1]/likelihood

        # M-step
        new_trans_mat = pd.DataFrame(np.zeros((len(states), len(states))), index = states, columns = states)
        for i in range(len(states)):
            for j in range(len(states)):
                new_trans_mat.iloc[i, j] = sum([eta_mat[i, j, t] for t in range(len(obs)-1)]) / np.sum(eta_mat[i, :, :-1])

        new_emit_mat = pd.DataFrame(np.zeros((len(states), len(emit_chars))), index = states, columns = emit_chars)
        for j in range(len(states)):
            for k in range(len(emit_chars)):
                new_emit_mat.iloc[j, k] = sum([gamma_mat[j, t] for t in range(len(obs)) if obs[t] == emit_chars[k]])/np.sum(gamma_mat[j, :])

        trans_mat = new_trans_mat
        emit_mat = new_emit_mat
        print(num_iter)
        num_iter+=1

    write_matrix_to_file('data/outputs/ba10k_rslt.txt', states, emit_chars, trans_mat, emit_mat)

    return (trans_mat, emit_mat)


# http://rosalind.info/problems/ba10h/
def get_hmm_path_obs(file):
    with open(file, 'r') as f:
        lst = f.readlines()
    obs = lst[0][:-1]
    emit_chars = [x for x in re.split(r'\s', lst[2]) if x != '']
    path = lst[4][:-1]
    states = [x for x in re.split(r'\s', lst[6]) if x != '']
    return (path, obs, states, emit_chars)


def hmm_param_est(path, obs, states, emit_chars, file = True):
    '''
    see http://cs229.stanford.edu/section/cs229-hmm.pdf
    use Lagrange multiplier to solve for MLE
    '''
    assert len(path) == len(obs)
    trans_mat = pd.DataFrame(np.zeros((len(states), len(states))), index = states, columns = states)
    emit_mat = pd.DataFrame(np.zeros((len(states), len(emit_chars))), index = states, columns = emit_chars)

    for i in range(len(obs)-1):
        trans_mat.loc[path[i], path[i+1]] += 1
        emit_mat.loc[path[i], obs[i]] += 1
    emit_mat.loc[path[-1], obs[-1]] += 1

    for state in states:
        if sum(trans_mat.loc[state, :]) > 0:
            trans_mat.loc[state, :] /= sum(trans_mat.loc[state, :])
        if sum(emit_mat.loc[state, :]) > 0:
            emit_mat.loc[state, :] /= sum(emit_mat.loc[state, :])
    if file:
        write_matrix_to_file('data/outputs/ba10h_rslt.txt', states, emit_chars, trans_mat, emit_mat)
    return (trans_mat, emit_mat)


# http://rosalind.info/problems/ba10i/
def viterbi_learning(file):

    iter_lim, obs, emit_chars, states, trans_mat, emit_mat = get_hmm_learning_info(file)

    num_iter = 0
    print(iter_lim)
    while num_iter < iter_lim:

        print(num_iter)

        path = viterbi(trans_mat, emit_mat, obs, states, False) # compute best hidden path
        trans_mat, emit_mat = hmm_param_est(path, obs, states, emit_chars, file = False) # new estimation
        num_iter += 1

    write_matrix_to_file('data/outputs/ba10i_rslt.txt', states, emit_chars, trans_mat, emit_mat)


# http://rosalind.info/problems/ba10e/
# http://rosalind.info/problems/ba10f/

def profile_hmm_info(file, pseudocount = False):

    with open(file, 'r') as f:
        lst = f.readlines()

    emit_chars = [x for x in re.split(r'\s', lst[2]) if x != '']
    obs = [x[:-1] for x in lst[4:]]

    if not pseudocount:
        thres = float(lst[0][:-1])
        return (thres, emit_chars, obs)
    else:
        figures = [x for x in re.split(r'\s', lst[0]) if x != '']
        thres = float(figures[0])
        pseudo = float(figures[1])
        return (thres, pseudo, emit_chars, obs)


def profile_hmm_matrix(file, pseudocount = False, write_file = True):

    if not pseudocount:
        thres, emit_chars, alignments = profile_hmm_info(file)
    else:
        thres, pseudo, emit_chars, alignments = profile_hmm_info(file, pseudocount)

    cols_to_keep = []
    num_cols = len(alignments[0])

    for i in range(num_cols):
        col = [x[i] for x in alignments]
        if np.mean([x == '-' for x in col]) < thres: # strictly less than
            cols_to_keep.append(i)

    # Determine all possible states

    states = ['S', 'I0']
    for i in range(1, len(cols_to_keep)+1):
        states.extend([x+str(i) for x in ['M', 'D', 'I']])
    states.append('E')

    trans_mat = np.zeros((len(states), len(states)))
    emit_mat = np.zeros((len(states), len(emit_chars)))

    trans_mat = pd.DataFrame(trans_mat, index = states, columns = states)
    emit_mat = pd.DataFrame(emit_mat, index = states, columns = emit_chars)

    for seq in alignments:

        state = 'S'
        col_idx = 0

        while state != 'E':

            if state == 'S':
                curr_idx = 0
            else:
                curr_idx = int(state[1:])

            if col_idx == len(seq):
                nxt_state = 'E'

            elif col_idx in cols_to_keep:

                # a standard seq has len(col_to_keep) characters
                # so the position in col_to_keep cannot generate an insertion

                if seq[col_idx] != '-':
                    nxt_state = 'M' + str(curr_idx+1)
                else:
                    nxt_state = 'D' + str(curr_idx+1)

            elif col_idx not in cols_to_keep and seq[col_idx] != '-':
                nxt_state = 'I' + str(curr_idx) # if col_idx is not conserved --> insertion

            else: # col_idx not in cols_to_keep and seq[col_idx] == '-'
                nxt_state = state

            print(seq, state, nxt_state, col_idx)

            # Mi cannot transit to itself, but Ii can, and we need to make sure the next one is not a '-' to make the insertion valid
            if nxt_state != state or (state[0] == 'I' and seq[col_idx] != '-'):
                trans_mat.loc[state, nxt_state] += 1
            if state[0] in 'MI' and seq[col_idx-1] != '-': # make sure the current one is not a '-'
                emit_mat.loc[state, seq[col_idx-1]] += 1

            state = nxt_state
            col_idx += 1

    for i in range(len(states)):
        if sum(trans_mat.iloc[i,:]) != 0:
            trans_mat.iloc[i,:] /= sum(trans_mat.iloc[i,:])
        if sum(emit_mat.iloc[i,:]) != 0:
            emit_mat.iloc[i,:] /= sum(emit_mat.iloc[i,:])

    if pseudocount:

        # consider all valid transitions, and we will later add pseudo to each valid cases and normalize

        '''
        consider a HMM with two states, S, I0, M1, I1, D1, M2, I2, D2, E

        we can have S --> I0, M1, D1
        we can have I0 --> I0, M1, D1
        we can have I1 --> I1, M2, D2
        we can have I2 --> I2, S
        we can have M1, D1 --> I1, M2, D2
        we can have M2, D2 --> I2, S

        Thus, we need to separate the case of S,I0 from the other cases
        and the case of I2,M2,D2 from the other cases
        '''

        valid_transitions = []
        valid_transitions.extend([('S', x) for x in ['M1', 'I0', 'D1']] + [('I0', x) for x in ['M1', 'I0', 'D1']])
        for s in ['M', 'I', 'D']:
            for i in range(1, len(cols_to_keep)):
                valid_transitions.extend([(s+str(i), x) for x in ['I'+str(i), 'M'+str(i+1), 'D'+str(i+1)]])
            valid_transitions.extend([(s+str(len(cols_to_keep)), x) for x in ['I'+str(len(cols_to_keep)), 'E']])

        print(valid_transitions)

        for state1 in states:
            pseudo_added = False
            for state2 in states:
                if (state1, state2) in valid_transitions:
                    trans_mat.loc[state1, state2] += pseudo
                    pseudo_added = True
            if pseudo_added:
                trans_mat.loc[state1,:] /= sum(trans_mat.loc[state1,:]) # normalize

        for state1 in states:
            if state1[0] in 'MI':
                for char in emit_chars:
                    emit_mat.loc[state1, char] += pseudo
                emit_mat.loc[state1,:] /= sum(emit_mat.loc[state1,:]) # normalize

    '''
    the output file has to be in the exact format as follows (see write_matrix_to_file function in utils):
    '''
    if write_file:  # msa_hmm (ba10g) calls this with write_file = False
        write_matrix_to_file('data/outputs/ba10f_rslt.txt', states, emit_chars, trans_mat, emit_mat)

    return (trans_mat, emit_mat)


# http://rosalind.info/problems/ba10g/
class ViterbiGraph: # it is easier to use graph structure rather than matrix to represent profile HMM
    '''
    adapted from https://github.com/egeulgen/Bioinformatics_Textbook_Track/blob/master/solutions/BA10G.py
    the adapted Viterbi algorithm pseudocode: https://profs.info.uaic.ro/~ciortuz/SLIDES/profileHMM.pdf
    '''
    def __init__(self):
        self.all_nodes = {}
        self.all_edges = []
        self.root = self.add_node(label = 'S', score = 1)

    class node:

        def __init__(self, label = None, score = -1e4):
            self.label = label
            self.edges = []
            self.parents = []
            self.children = []
            self.score = score
            self.traceback = None

    class edge:

        def __init__(self, from_node = None, to_node = None):
            self.from_node = from_node
            self.to_node = to_node

    def add_node(self, label, score = -1e4):

        if label in self.all_nodes.keys():
            return self.all_nodes[label]

        newNode = ViterbiGraph.node(label = label, score = score)

        self.all_nodes[label] = newNode

        return newNode

    def add_edge(self, from_node, to_node):

        newEdge = ViterbiGraph.edge(from_node, to_node)
        #print(newEdge.from_node.label, newEdge.to_node.label)
        from_node.edges.append(newEdge)
        from_node.children.append(to_node)
        to_node.parents.append(from_node)
        self.all_edges.append(newEdge)

    def construct_graph(self, obs, trans_mat):

        states = trans_mat.columns
        l = int(states[-2][1:])

        valid_nodes = ['S']
        valid_nodes.extend([(0, s) for s in states if s[0] == 'D']) # only D can match empty string
        # the first item indicates the posItIon in obs, where we adopt 1-indexing instead of 0-indexing
        for i in range(1, len(obs)+1):
            valid_nodes.extend([(i, s) for s in states if s not in 'SE'])
        valid_nodes.append('E')

        '''
        valid transitions include and only include S to (0,D1), (1,I0), (1,M1)
        (t,Dj),(t,Ij),(t,Mj) to (t,Dj+1), (t+1,Ij), (t+1,Mj+1)
        and (T,Dl),(T,Il),(T,Ml) to E where T == len(obs) and l == len(cols_to_keep) in profile_hmm_matrix
        '''

        valid_transitions = []
        for node in valid_nodes:
            if node == 'S':
                valid_transitions.extend([('S', x) for x in [(0,'D1'), (1,'I0'), (1,'M1')]])
                # we need to put D first, because we need to consider the position 0 of string before position 1
                # this sequence is important in BFS traversal
            elif node == 'E':
                continue
            else:
                curr_idx = node[0]
                curr_state = node[1]
                state = curr_state[0]
                state_idx = int(curr_state[1:])

                if curr_idx == len(obs) and state_idx == l:
                    valid_transitions.append((node, 'E'))
                elif state_idx == l:
                    valid_transitions.append((node, (curr_idx+1, 'I'+str(l))))
                elif curr_idx == len(obs):
                    # the whole string is consumed: only trailing deletions (D) can still lead to E
                    valid_transitions.append((node, (curr_idx, 'D'+str(state_idx+1))))
                else:
                    next_nodes = [(curr_idx, 'D'+str(state_idx+1)),
                                  (curr_idx+1, 'I'+str(state_idx)),
                                  (curr_idx+1, 'M'+str(state_idx+1)),
                                  ]
                    valid_transitions.extend([(node, x) for x in next_nodes])

        #for node1, node2 in valid_transitions:
            #print(str(node1) + ' --> ' + str(node2))

        for node in valid_nodes:
            self.add_node(node)

        # print(self.all_nodes.keys())

        for state1, state2 in valid_transitions:
            self.add_edge(self.all_nodes[state1], self.all_nodes[state2])

        return

    def bfs(self, obs, trans_mat, emit_mat):

        visited = {node: False for node in self.all_nodes.keys() if node != 'E' and node != 'S'}

        queue = [edge.to_node.label for edge in self.all_nodes['S'].edges]

        while queue:

            s = queue.pop(0)
            # print(visited[s])
            curr_node = self.all_nodes[s]
            visited[s] = True

            string_idx = s[0] - 1
            curr_state = s[1]

            for parent in curr_node.parents:
                p = parent.label
                if p == 'S':
                    prev_state = 'S'
                else:
                    prev_state = p[1]

                if curr_state.startswith('D'):
                    tmp_score = parent.score * trans_mat.loc[prev_state, curr_state]
                else:
                    tmp_score = parent.score * trans_mat.loc[prev_state, curr_state] * emit_mat.loc[curr_state, obs[string_idx]]
                if tmp_score > curr_node.score:
                    curr_node.score = tmp_score
                    curr_node.traceback = parent

            for nxt_node in curr_node.children:
                if nxt_node.label != 'E': # isolate the case of E
                    queue.append(nxt_node.label)

            queue = [s for s in queue if not visited[s]] # make sure all nodes in queue have not been visited yet

        return

    def alternative_bfs(self, obs, trans_mat, emit_mat):
        '''
        we can also traverse by string_idx
        '''
        for string_idx in range(len(obs)+1):

            # find nodes in current column
            curr_column = []
            for node in self.all_nodes.values():
                if node.label[0] == string_idx:
                    curr_column.append(node)

            # calculate scores for each node in column
            for curr_node in curr_column:

                if curr_node.label != "S":
                    curr_state = curr_node.label[1]

                    for parent in curr_node.parents:
                        #print(parent_node.label)
                        prev_state = 'S' if parent.label == 'S' else parent.label[1]
                        #print(parent_state, current_state)
                        if curr_state.startswith('D'):
                            tmp_score = parent.score * trans_mat.loc[prev_state, curr_state]
                        else:
                            tmp_score = parent.score * trans_mat.loc[prev_state, curr_state] * emit_mat.loc[curr_state, obs[string_idx-1]]

                        if tmp_score > curr_node.score:
                            curr_node.score = tmp_score
                            curr_node.traceback = parent
        return

    def find_optimal_path(self, obs, trans_mat, emit_mat):

        self.construct_graph(obs, trans_mat)

        self.bfs(obs, trans_mat, emit_mat)

        #self.alternative_bfs(obs, trans_mat, emit_mat)

        # termination
        curr_node = self.all_nodes['E']
        for parent in curr_node.parents:
            p = parent.label
            prev_state = p[1]
            curr_state = 'E'
            # print(prev_state, curr_state, trans_mat.loc[prev_state, curr_state], parent.score)
            tmp_score = parent.score * trans_mat.loc[prev_state, curr_state]
            if tmp_score > curr_node.score:
                curr_node.score = tmp_score
                curr_node.traceback = parent

        last_node = self.all_nodes['E'].traceback
        print(self.all_nodes['E'].score)

        optimal_hidden_path = []
        # print(self)

        while last_node.label != 'S':
            # print(current_node.label)
            optimal_hidden_path = [last_node.label[1]] + optimal_hidden_path
            last_node = last_node.traceback
        return optimal_hidden_path

    def __str__(self):
        ''' Custom __str__ method
        For printing the Viterbi Graph as an adjacency list
        '''
        # header
        res = 'From --> To: score(target)\n\n'
        for edge in self.all_edges:
            res += str(edge.from_node.label) + ' --> ' + str(edge.to_node.label) + ': ' + str(edge.to_node.score) + '\n'
        return res

    def __repr__(self):
        return '<Viterbi Graph adjacency list representation>'


def msa_hmm(file):

    with open(file, 'r') as f:
        lst = f.readlines()

    obs = lst[0][:-1]

    emit_chars = [x for x in re.split(r'\s', lst[4]) if x != '']
    alignments = [x[:-1] for x in lst[6:]]

    figures = [x for x in re.split(r'\s', lst[2]) if x != '']
    thres = float(figures[0])
    pseudo = float(figures[1])

    # we then write a file for profile_hmm_matrix
    with open('data/outputs/ba10g_input1.txt', 'w') as f:
        f.write(str(thres) + '\t' + str(pseudo)+'\n')
        f.write('--------\n')
        f.write('\t'.join(emit_chars)+'\n')
        f.write('--------\n')
        for seq in alignments:
            f.write(seq + '\n')

    transition_matrix, emission_matrix = profile_hmm_matrix('data/outputs/ba10g_input1.txt', True, False)
    viterbi_graph = ViterbiGraph()
    #viterbi_graph.construct_graph(obs, transition_matrix)
    optimal_hidden_path = viterbi_graph.find_optimal_path(obs, transition_matrix, emission_matrix)
    print(' '.join(optimal_hidden_path))
    # print(viterbi_graph)
    return optimal_hidden_path


if __name__ == "__main__":
    markov_prob('data/inputs/rosalind_ba10a.txt')

    markov_prob_2('data/inputs/rosalind_ba10b.txt')

    viterbi(*precompute_hmm_info('data/inputs/rosalind_ba10c.txt'))

    hmm_likelihood('data/inputs/rosalind_ba10d.txt')

    soft_decoding()

    baum_welch('data/inputs/rosalind_ba10k.txt')

    hmm_param_est(*get_hmm_path_obs('data/inputs/rosalind_ba10h.txt'), True)

    viterbi_learning('data/inputs/rosalind_ba10i.txt')

    profile_hmm_matrix('data/inputs/rosalind_ba10e.txt')

    profile_hmm_matrix('data/inputs/rosalind_ba10f.txt', True)

    msa_hmm('data/inputs/rosalind_ba10g.txt')
