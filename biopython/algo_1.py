import matplotlib.pyplot as plt
#------------------------------------------------------------

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()
            seq = fh.readline().rstrip()
            fh.readline()
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def phred33ToQ(qual):
    return ord(qual) - 33

def QToPhred33(qual):
    return chr(qual+33)

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences


def naive_2mm(p, t):
    occurrences = []
    mm = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mm +=1
            if mm > 2:
                mm = 0
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
            mm = 0
    return occurrences


def naive_with_rc(p, t):
    occurrences = []
    r = reverseComplement(p)
    if r == p:
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
    else:
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
        for i in range(len(t) - len(r) + 1):  # loop over alignments
            match = True
            for j in range(len(r)):  # loop over characters
                if t[i+j] != r[j]:  # compare characters
                    match = False
                    break
            if match:
                occurrences.append(i)  # all chars matched; record
    return occurrences


def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = phred33ToQ(phred)
            hist[q] +=1
    plt.bar(range(len(hist)), hist)
    plt.show()


#------------------------------------------------------------

# p = 'CCC'
# ten_as = 'AAAAAAAAAA'
# t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
# occurrences = naive_with_rc(p, t)
# print(occurrences)
# #------------------------------------------------------------
# p = 'CGCG'
# t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
# occurrences = naive_with_rc(p, t)
# print(occurrences)
# #------------------------------------------------------------
# phix_genome = readGenome('phix.fa')
# occurrences = naive_with_rc('ATTA', phix_genome)
# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))
#------------------------------------------------------------
# lambda_genome = readGenome('lambda_virus.fa')
# occurrences = naive_with_rc('AGTCGA', lambda_genome)
# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))
#------------------------------------------------------------
# p = 'CTGT'
# ten_as = 'AAAAAAAAAA'
# t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as
# occurrences = naive_2mm(p, t)
# print(occurrences)
#------------------------------------------------------------
# phix_genome = readGenome('phix.fa')
# occurrences = naive_2mm('GATTACA', phix_genome)

# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))

#------------------------------------------------------------
# lambda_genome = readGenome('lambda_virus.fa')
# occurrences = naive_2mm('AGGAGGTT', lambda_genome)
# print('offset of leftmost occurrence: %d' % min(occurrences))
# print('# occurrences: %d' % len(occurrences))
#------------------------------------------------------------
# print(phred33ToQ(";"))
# print(QToPhred33(26))
#------------------------------------------------------------

# seqs, qualities = readFastq("ERR037900_1.first1000.fastq")
# # createHist(qualities)
# cycles = []
# for qual in qualities:
#     cycles.append(qual.find(min(qual)))
        
# print(len(cycles))
# plt.bar(range(len(cycles)), cycles)
# plt.show()

#------------------------------------------------------------

