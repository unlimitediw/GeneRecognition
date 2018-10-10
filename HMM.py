# PART A: Describe and/or draw the HMM for this problem. Define the states, transition probabilities,
# and the emission probabilities inferred from the training data above.

# HMM description:
# In this problem, we have two states: N(non-coding) and G(coding-gene).
# In this Hidden Markov Models: the transition probabilities are defined by us and the emission probabilities are given.
# We need to use these values and chain rule in HMM to find the probability of N and G in each position of fasta
# sequence.


# PART OF INITIALIZATION: STATES, TRANSITION PROBABILITY, EMISSION PROBABILITY AND FASTA SEQUENCE
# definition of states, 'N' is non-coding, 'G' is coding.

# transition probability of N to N, N to G, G to G, G to N (with normalization)
pNN = 0.9
pNG = 0.1
pGG = 0.9
pGN = 0.1

# emission probability for node-coding N (with normalization)
non_code_em_pro = 1.

# emission probability for coding gene G（with normalization)
code_em_pro = [[[77 * 4 / 184, 66 * 4 / 183, 37 * 4 / 184, 4 * 4 / 184],
                [3 * 4 / 98, 63 * 4 / 98, 13 * 4 / 98, 19 * 4 / 98],
                [1 * 4 / 28, 23 * 4 / 28, 1 * 4 / 28, 3 * 4 / 28],
                [1 * 4 / 188, 98 * 4 / 188, 60 * 4 / 188, 29 * 4 / 188]],
               [[15 * 4 / 116, 23 * 4 / 116, 73 * 4 / 116, 5 * 4 / 116],
                [11 * 4 / 76, 1 * 4 / 76, 55 * 4 / 76, 9 * 4 / 76],
                [1 * 4 / 137, 46 * 4 / 137, 1 * 4 / 137, 89 * 4 / 137],
                [1 * 4 / 171, 18 * 4 / 171, 141 * 4 / 171, 11 * 4 / 171]],
               [[147 * 4 / 338, 85 * 4 / 338, 46 * 4 / 338, 60 * 4 / 338],
                [30 * 4 / 128, 19 * 4 / 128, 49 * 4 / 128, 30 * 4 / 128],
                [1 * 4 / 131, 47 * 4 / 131, 5 * 4 / 131, 78 * 4 / 131],
                [34 * 4 / 144, 21 * 4 / 144, 34 * 4 / 144, 55 * 4 / 144]],
               [[0, 38 * 4 / 56, 0, 18 * 4 / 56],
                [2 * 4 / 77, 38 * 4 / 77, 5 * 4 / 77, 32 * 4 / 77],
                [0, 5 * 4 / 18, 8 * 4 / 18, 5 * 4 / 18],
                [2 * 4 / 69, 44 * 4 / 69, 8 * 4 / 69, 15 * 4 / 69]]]

# given fasta sequence
fasta1 = open("1000_2000_1.txt", "r").read()
fasta2 = open("1000_5000.txt", "r").read()
fasta3 = open("1_39675.txt", "r").read()


# implement a,c,g,t in 0,1,2,3 respectively
def fasta_to_fasta_list(fasta):
    fasta_list = []
    for idx in range(len(fasta)):
        if fasta[idx] == 'A':
            fasta_list.append(0)
        elif fasta[idx] == 'C':
            fasta_list.append(1)
        elif fasta[idx] == 'G':
            fasta_list.append(2)
        elif fasta[idx] == 'T':
            fasta_list.append(3)
    return fasta_list


fasta_list_1 = fasta_to_fasta_list(fasta1)
fasta_list_2 = fasta_to_fasta_list(fasta2)
fasta_list_3 = fasta_to_fasta_list(fasta3)


# PART B: show i) the most likely sequence of states for the genomic sequence provided,
# ii) the joint probability for the path and sequence (path score), and
# iii) the number and the start-end coordinates of the predicted genes in the sequence provided.

# i): Find Max(P(X1, E1,..., Xt, Et))
# For Markov Models, by joint distribution rule we know that:
# P(X1, E1,..., Xt, Et) = P(X1)P(E1|X1)Pi(2 to T)P(Xt|Xt-1)P(Et|Xt)
# Without any optimization, the complexity of this problem is O(2^n)
# With the Viterbi Algorithm, the complexity is O(n)
# As we know, there are two states: N and G
# At first position, we know pN0 = 1, pG0 = 0
# At the second position pN1 = max(pN0 * pNN, pG0 * pGN) * code_em_pro(first position state is N),
# pG1 = max(pN0 * pNG, pG0 * pGG) * non_code_em_pro(first position N).
# Moreover, for t position: pNt = max(pNt-1 * pNN, pGt-1 * pGN) * code_em_pro(sequence of states in t-1 position)
# pGt = max(pNt-1 * pNG, pGt-1 * pGG) * non_code_em_pro(sequence of states in t-1 position)
# Since pNt and pGt both are direct proportion to pNt-1 and pGt-1 thus this dynamic programing is appropriate to
# find the most likely sequence of states of fasta
def solution(fasta_list):
    print("fasta_list(A-0, C-1, G-2, T-3):", fasta_list)
    # add 'AA' at the start position of fasta
    fasta_list.insert(0, 0)
    fasta_list.insert(0, 0)
    # the state of first position is N
    states_N = ['N']
    states_G = ['N']
    probabilities_N = [1.]
    probabilities_G = [0.]
    # HMM with Viterbi Algorithm(dynamic programming)
    for t in range(2, len(fasta_list)):
        # at position t, the probability of N
        if probabilities_N[-1] * pNN > probabilities_G[-1] * pGN:
            new_state_N = states_N[:]
            probabilities_N.append(
                probabilities_N[-1] * pNN * non_code_em_pro)
        else:
            new_state_N = states_G[:]
            probabilities_N.append(
                probabilities_G[-1] * pGN * non_code_em_pro)
        new_state_N.append('N')

        # at position t, the probability of G
        if probabilities_G[-1] * pGG > probabilities_N[-1] * pNG:
            new_state_G = states_G[:]
            probabilities_G.append(
                probabilities_G[-1] * pGG * code_em_pro[fasta_list[t - 2]][fasta_list[t - 1]][fasta_list[t]])
        else:
            new_state_G = states_N[:]
            probabilities_G.append(
                probabilities_N[-1] * pNG * code_em_pro[fasta_list[t - 2]][fasta_list[t - 1]][fasta_list[t]])

        probabilities_N[-1] *= 1.1
        probabilities_G[-1] *= 1.1
        new_state_G.append('G')
        states_N = new_state_N
        states_G = new_state_G

    # use to find total number of G and its start end coordinates
    def coordinate(states_list):
        coordinate_list = []
        cur = 'N'
        count = 0
        start = 0
        for s_idx in range(len(states_list)):
            if states_list[s_idx] == 'G':
                count += 1
                if cur == 'N':
                    start = s_idx
                    cur = 'G'
            elif cur == 'G':
                end = s_idx - 1
                cur = 'N'
                coordinate_list.append((start, end))
        return count, coordinate_list

    if probabilities_N[-1] > probabilities_G[-1]:
        return states_N, probabilities_N[-1], coordinate(states_N)
    else:
        return states_G, probabilities_G[-1], coordinate(states_G)


fasta_test = fasta_list_3
# s is the most likely states of fasta, p is the probability of this states sequence, c is the number of gene and
# coordinates of genes starts and ends
s, p, c = solution(fasta_test)
print("Correspond states sequence(only use the alpha3 part to represent):", s)
print("Probability of this states sequence:", p)
print("Total number of gene code:", c[0])
print("coordinates of gene code start and end(after insertion of 'AA'):", c[1])
# this will show the part of gene in fasta
for i in c[1]:
    print(fasta_test[i[0]:i[1] + 3])
