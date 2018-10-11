def faster_symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    extended_genome = genome + genome[0:n // 2]

    # look at the first half of Genome to compute first array value
    array[0] = get_pattern_count(symbol, genome[0:n // 2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i - 1]

        # the current array value can differ from the previous array value by at most 1
        if extended_genome[i - 1] == symbol:
            array[i] = array[i] - 1
        if extended_genome[i + (n // 2) - 1] == symbol:
            array[i] = array[i] + 1
    return array


def get_symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    ExtendedGenome = genome + genome[0:n // 2]
    for i in range(n):
        window = ExtendedGenome[i:i + (n // 2)]
        array[i] = get_pattern_count(symbol, window)
    return array


def get_pattern_count(pattern, text):
    count = 0
    for i in range(len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            count += 1
    return count


def get_minimum_skew(genome):
    positions = []
    skew = get_skew_array(genome)
    minimum = 999999999
    for i in range(len(skew)):
        if skew[i] == minimum:
            positions.append(i)
        elif skew[i] < minimum:
            minimum = skew[i]
            del positions[:]
            positions.append(i)
    return positions


def get_skew_array(genome):
    skew = {}
    count = len(genome)
    if genome[0] == 'G':
        skew[0] = 1
    elif genome[0] == 'C':
        skew[0] = -1
    else:
        skew[0] = 0
    for i in range(1, count):
        if genome[i] == 'G':
            skew[i] = skew[i - 1] + 1
        elif genome[i] == 'C':
            skew[i] = skew[i - 1] - 1
        else:
            skew[i] = skew[i - 1]
    return skew


file = open('./resources/E_coli.txt', 'r')
genome = file.read()
file.close()
symbol = 'C'
# print(faster_symbol_array(genome, symbol))
print(get_minimum_skew('CATTCCAGTACTTCGATGATGGCGTGAAGA'))
