import itertools
from algo_1 import readFastq

#------------------------------------------------------------
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match
#------------------------------------------------------------
# def scs(ss):
#     """ Returns shortest common superstring of given
#         strings, which must be the same length """
#     shortest_sup = None
#     count = 0
#     scs_l = []
#     for ssperm in itertools.permutations(ss):
#         # print(ssperm)
#         sup = ssperm[0]  # superstring starts as first string
#         for i in range(len(ss)-1):
#             # overlap adjacent strings A and B in the permutation
#             olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
#             # add non-overlapping portion of B to superstring
#             sup += ssperm[i+1][olen:]
#         if shortest_sup is not None and len(sup) == len(shortest_sup):
#             scs_l.append(sup)
#             count +=1
#         if shortest_sup is None or len(sup) < len(shortest_sup):
#             scs_l.append(sup)
#             shortest_sup = sup  # found shorter superstring
#     m = len(min(scs_l, key=len))
#     filtered = []
#     for e in scs_l:
#         if len(e) <= m:
#             filtered.append(e)
        
#     print(count)
#     return sorted(filtered)  # return shortest
#------------------------------------------------------------
def make_kmers_dict(reads, kmer_len):
    kmers = {}
    for read in reads:
        i = 0
        while (i + kmer_len) < len(read)+1:
            kmer = read[i:i+kmer_len]
            if kmer not in kmers:
                kmers[kmer] = set()
            kmers[kmer].add(read)
            i+=1
    return kmers

def find_overlap(reads, k):
    kmers = make_kmers_dict(reads, k)
    overlapping= []

    for r in reads:
        suffix = r[-k:]
        if len(kmers[suffix]) > 1:
            for b in kmers[suffix]:
                # print(b, kmers[suffix])
                if b != r:
                    overlapping.append((r, b))
    
    return overlapping

def pick_max_overlap(reads, k):
    read_a, read_b = None, None
    best_olen = 0
    for a, b in find_overlap(reads, k):
        olen = overlap(a,b, min_length=k)
        if olen > best_olen:
            read_a, read_b = a, b
            best_olen = olen
    return read_a, read_b, best_olen
#------------------------------------------------------------
def greedy_scs(reads, k):
    read_a, read_b, olen = pick_max_overlap(reads, k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_max_overlap(reads, k)
    return "".join(reads)
#------------------------------------------------------------
# strings = ['ABC', 'BCA', 'CAB']
# print(greedy_scs(strings, 2))
# strings = ['ABCD', 'CDBC', 'BCDA']
# print(greedy_scs(strings, 1))
# print(scs(strings))
#------------------------------------------------------------
# strings = ['ABC', 'BCA', 'CAB']
# strings = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]
# out = scs(strings)
# print(len(out))
# print(out)
#------------------------------------------------------------
# strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']
# print(scs(strings))
#------------------------------------------------------------
seqs, _ = readFastq("data/ads1_week4_reads.fq")
out = greedy_scs(seqs, 50)
print(len(out))
print()
final_A = out.count("A")
final_T = out.count("T")
print(final_A)
print(final_T)
#------------------------------------------------------------
