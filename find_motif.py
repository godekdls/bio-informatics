def get_symbol_count_dictionary(motifs):
    # init count dictionary
    count = {}
    motif_length = len(motifs[0])  # suppose every string in motifs has same length
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(motif_length):
            count[symbol].append(0)

    motif_count = len(motifs)
    for i in range(motif_count):
        for j in range(motif_length):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    return count


def get_profile(motifs):
    count = get_symbol_count_dictionary(motifs)
    motif_count = len(motifs)
    motif_length = len(motifs[0])

    profile = {}
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(motif_length):
            profile[symbol].append(count[symbol][j] / motif_count)

    return profile


# get the most frequent nucleotide from each bp and finally gather them
def get_candidate_regulatory_motif(motifs):
    count = get_symbol_count_dictionary(motifs)
    motif_length = len(motifs[0])

    candidate = ""
    for j in range(motif_length):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        candidate += frequentSymbol
    return candidate


def get_mismatch_score(motifs):
    candidate = get_candidate_regulatory_motif(motifs)
    motif_count = len(motifs)
    motif_length = len(motifs[0])

    score = 0
    for i in range(motif_count):
        for j in range(motif_length):
            if motifs[i][j] != candidate[j]:
                score += 1
    return score


def get_probability_of_motif(motif, profile):
    probability = 1
    for i in range(len(motif)):
        if motif[i] == 'A':
            probability = probability * profile['A'][i]
        elif motif[i] == 'C':
            probability = probability * profile['C'][i]
        elif motif[i] == 'G':
            probability = probability * profile['G'][i]
        else:
            probability = probability * profile['T'][i]
    return probability


# calculate probability using profile and pick the most likely candidate
# k : k-mer
def get_candidate_by_profile(gene, k, profile):
    max_probability = -1
    final_candidate = ''
    for i in range(len(gene) - k + 1):
        candidate = gene[i:i + k]
        probability = 1
        for j in range(k):
            if candidate[j] == 'A':
                probability = probability * profile['A'][j]
            elif candidate[j] == 'C':
                probability = probability * profile['C'][j]
            elif candidate[j] == 'G':
                probability = probability * profile['G'][j]
            else:
                probability = probability * profile['T'][j]
        if probability > max_probability:
            final_candidate = candidate
            max_probability = probability
    return final_candidate


# k : k-mer
def get_best_motifs_by_greedy_algorithm(genes, k):
    best_motifs = []
    for i in range(0, len(genes)):
        best_motifs.append(genes[i][0:k])
    for i in range(len(genes[0]) - k + 1):
        motifs = []
        motifs.append(genes[0][i:i + k])
        for j in range(1, len(genes)):
            profile = get_profile(motifs[0:j])
            candidate = get_candidate_by_profile(genes[j], k, profile)
            motifs.append(candidate)
        if get_mismatch_score(motifs) < get_mismatch_score(best_motifs):
            best_motifs = motifs
    return best_motifs


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
print(get_symbol_count_dictionary(motifs))
print(get_profile(motifs))
print(get_candidate_regulatory_motif(motifs))
print(get_mismatch_score(motifs))

profile = {
    'A': [0.2, 0.2, 0.3, 0.2, 0.3],
    'C': [0.4, 0.3, 0.1, 0.5, 0.1],
    'G': [0.3, 0.3, 0.5, 0.2, 0.4],
    'T': [0.1, 0.2, 0.1, 0.1, 0.2]
}
print(get_probability_of_motif('GGTAC', profile))
print(get_candidate_by_profile('ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT', 5, profile))

genes = [
    "GGCGTTCAGGCA",
    "AAGAATCAGTCA",
    "CAAGGAGTTCGC",
    "CACGTCAATCAC",
    "CAATAATATTCG"
]
print(get_best_motifs_by_greedy_algorithm(genes, 3))

genes = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
         "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
         "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
         "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
         "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
         "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
         "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
         "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
         "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
         "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
print(get_best_motifs_by_greedy_algorithm(genes, 15))