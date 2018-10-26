import random


def get_nucleotide_count_map(motifs):
    count = {}
    motif_length = len(motifs[0])  # suppose every string in motifs has same length
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(motif_length):
            count[symbol].append(1)

    motif_count = len(motifs)
    for i in range(motif_count):
        for j in range(motif_length):
            symbol = motifs[i][j]
            count[symbol][j] += 1

    return count


def get_profile_width_pseudo_counts(motifs):
    count = get_nucleotide_count_map(motifs)
    motif_count = len(motifs)
    motif_length = len(motifs[0])

    profile = {}
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(motif_length):
            profile[symbol].append(count[symbol][j] / (motif_count + 4))  # add 4 because of pseudo count

    return profile


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


# get the most frequent nucleotide from each bp and finally gather them
def get_candidate_regulatory_motif(motifs):
    count = get_nucleotide_count_map(motifs)
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


def greedy_motif_search_with_pseudo_counts(genes, k, t):
    best_motifs = []
    for i in range(0, len(genes)):
        best_motifs.append(genes[i][0:k])
    for i in range(len(genes[0]) - k + 1):
        motifs = []
        motifs.append(genes[0][i:i + k])
        for j in range(1, len(genes)):
            profile = get_profile_width_pseudo_counts(motifs[0:j])
            candidate = get_candidate_by_profile(genes[j], k, profile)
            motifs.append(candidate)
        if get_mismatch_score(motifs) < get_mismatch_score(best_motifs):
            best_motifs = motifs
    return best_motifs


def get_motifs_by_profile(profile, dnas):
    k = len(profile['A'])  # k-mer
    result = []
    for i in range(len(dnas)):
        max_probability = -1
        candidate = ''
        for j in range(len(dnas[i]) - k + 1):
            probability = 1
            for l in range(k):
                bp = dnas[i][j + l]
                probability = probability * profile[bp][l]

            if (max_probability < probability):
                max_probability = probability
                candidate = dnas[i][j:j + k]
        result.append(candidate)
    return result


def get_motif_by_profile(dna, profile, k):
    n = len(dna)
    probabilities = {}
    for i in range(0, n - k + 1):
        probabilities[dna[i:i + k]] = get_probability_of_motif(dna[i:i + k], profile)
    probabilities = normalize(probabilities)
    return roll_a_dice_with_probabilities(probabilities)


# dnas:  A list of strings Dna
# k -> k-mer
def get_random_motifs(dnas, k):
    result = []
    for i in range(len(dnas)):
        index = random.randint(0, len(dnas[0]) - k)  # inclusively
        result.append(dnas[i][index:index + k])
    return result


# k : k-mer
def get_motifs_by_randomized_search_algorithm(dnas, k):
    best_motifs = get_random_motifs(dnas, k)
    while True:
        profile = get_profile_width_pseudo_counts(best_motifs)
        motifs = get_motifs_by_profile(profile, dnas)
        if get_mismatch_score(motifs) < get_mismatch_score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs


# normalize probabilities so that sum of probabilities would be 1
def normalize(probabilities):
    sum = 0
    for key, value in probabilities.items():
        sum = sum + value
    normalized = {}
    for key, value in probabilities.items():
        normalized[key] = value / sum
    return normalized


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def roll_a_dice_with_probabilities(probabilities):
    decimal = random.uniform(0, 1)
    cur = 0
    for key, value in probabilities.items():
        cur = cur + value
        if cur >= decimal:
            return key
    return ''


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


# k : k-mer
# n : iteration number
def get_motifs_by_gibbs_sampler_algorithm(dnas, k, n):
    best_motifs = get_motifs_by_randomized_search_algorithm(dnas, k)
    for j in range(n):
        index = random.randint(0, len(dnas) - 1)  # inclusively
        motifs = []
        for i in range(len(dnas)):
            if i != index:
                motifs.append(best_motifs[i])
        profile = get_profile_width_pseudo_counts(motifs)
        motifs = get_motifs_by_profile(profile, dnas)
        if get_mismatch_score(motifs) < get_mismatch_score(best_motifs):
            best_motifs = motifs

    return best_motifs


motifs = ['AACGTA', 'CCCGTT', 'CACCTT', 'GGATTA', 'TTCCGG']
print(get_nucleotide_count_map(motifs))
print(get_profile_width_pseudo_counts(motifs))
dnas = ['GGCGTTCAGGCA', 'AAGAATCAGTCA', 'CAAGGAGTTCGC', 'CACGTCAATCAC', 'CAATAATATTCG']
print(greedy_motif_search_with_pseudo_counts(dnas, 3, 5))

profile = {'A': [0.8, 0.0, 0.0, 0.2], 'C': [0.0, 0.6, 0.2, 0.0],
           'G': [0.2, 0.2, 0.8, 0.0], 'T': [0.0, 0.2, 0.0, 0.8]}
dnas = ['TTACCTTAAC', 'GATGTCTGTC', 'ACGGCGTTAG', 'CCCTAACGAG', 'CGTCAGAGGT']
print(get_motifs_by_profile(profile, dnas))
print(get_random_motifs(dnas, 3))
dnas = ['GCGTAATAAGCTAATGGTAAAGCGTA', 'ATAAGCTAATGCCTTAGTAAAGCGTA', 'TATAGTCTATCGACGCCCCGGAAGCT',
        'TCCATACACTCGGAACGAAGCACCTC', 'TTGAGGAAAACACAGTCCTCTTGTGG', 'GGCCTCTCATGTAAACGCATCGATAT',
        'ATTCTAGATTTACAGACCGGCCCCGG', 'ATTTTATCCCAAGTTATAGGCCTTGA']
print(get_motifs_by_randomized_search_algorithm(dnas, 5))

dnas = [
    'GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC',
    'CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG',
    'ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC',
    'GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC',
    'GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG',
    'CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA',
    'GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA',
    'GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG',
    'GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG',
    'TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC']

probabilities = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
print(normalize(probabilities))
print(roll_a_dice_with_probabilities(probabilities))
profile = {'A': [0.5, 0.1], 'C': [0.3, 0.2], 'G': [0.2, 0.4], 'T': [0.0, 0.3]}
print(get_motif_by_profile('AAACCCAAACCC', profile, 2))

dna = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
       'TAGTACCGAGACCGAAAGAAGTATACAGGCGT', 'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']
k = 8
n = 100
best_motifs = []
for i in range(20):
    best_motifs = get_motifs_by_gibbs_sampler_algorithm(dna, k, n)
print(best_motifs)
