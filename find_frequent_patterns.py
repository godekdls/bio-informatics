def get_pattern_count(Text, Pattern):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if Text[i:i + len(Pattern)] == Pattern:
            count += 1
    return count


def get_pattern_matching_positions(Pattern, Genome):
    positions = []
    k = len(Pattern)
    n = len(Genome)
    for i in range(n - k + 1):
        kmer = Genome[i:i + k]
        if Pattern == kmer:
            positions.append(i)
    return positions


def get_frequent_patterns(Text, k):
    words = []
    freq = get_frequency_map(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def get_frequency_map(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        freq[Pattern] = 0
    for Pattern, value in freq.items():
        for i in range(len(Text) - len(Pattern) + 1):
            if Text[i:i + k] == Pattern:
                freq[Pattern] += 1
    return freq


def get_reverse_complement(Pattern):
    Pattern = get_reverse_str(Pattern)
    return get_complement(Pattern)


def get_reverse_str(Pattern):
    result = '';
    for char in Pattern:
        result = char + result
    return result;


def get_complement(Pattern):
    result = ''
    for char in Pattern:
        if char == 'A':
            result = result + 'T'
        elif char == 'T':
            result = result + 'A'
        elif char == 'C':
            result = result + 'G'
        else:
            result = result + 'C'
    return result


def get_hamming_distance(p, q):  # Two strings of equal length
    hamming_distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            hamming_distance += 1
    return hamming_distance


def get_approximate_pattern_matching(Text, Pattern, d):
    positions = []
    pattern_len = len(Pattern)
    for i in range(len(Text) - pattern_len + 1):
        candidate = Text[i:i + pattern_len]
        hamming_distance = get_hamming_distance(candidate, Pattern)
        if (hamming_distance <= d):
            positions.append(i)
    return positions

def get_approximate_pattern_count(Pattern, Text, d):
    pattern_matching = get_approximate_pattern_matching(Text, Pattern, d)
    return len(pattern_matching)


Text = 'CGCCTAAATAGCCTCGCGGAGCCTTATGTCATACTCGTCCT'
k = 3
frequent_patterns = get_frequent_patterns(Text, k)
print(frequent_patterns)
print(get_reverse_complement('CCAGATC'))
print(get_hamming_distance('CTTGAAGTGGACCTCTAGTTCCTCTACAAAGAACAGGTTGACCTGTCGCGAAG', 'ATGCCTTACCTAGATGCAATGACGGACGTATTCCTTTTGCCTCAACGGCTCCT'))
