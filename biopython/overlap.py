from algo_1 import readFastq

#------------------------------------------------------------
def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

#------------------------------------------------------------
#call overlap only on those we know exist
#every kmer has its own set of reads where that kmer occurs
#make dict to conatin these sets

def make_kmers_dict(reads, kmer_len):
    kmers = {}
    for read in reads:
        # print("read", read)
        
        i = 0
        while (i + kmer_len) < len(read)+1:
            kmer = read[i:i+kmer_len]
            # print(kmer)
            if kmer not in kmers:
                kmers[kmer] = set()
            kmers[kmer].add(read)
            i+=1

    return kmers
#------------------------------------------------------------
# find suffix of a as kmer and lookup reads in dict and call overlap on each
def find_overlap(reads, k):
    kmers = make_kmers_dict(reads, k)
    overlapping= []

    for r in reads:
        suffix = r[-k:]
        if len(kmers[suffix]) > 1:
            for b in kmers[suffix]:
                # print(b, kmers[suffix])
                if b != r:
                    if overlap(r, b, min_length=k) > 0:
                        overlapping.append((r, b))
    return overlapping

#------------------------------------------------------------
# reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']
# print(find_overlap(reads, 3))
# print(find_overlap(reads, 4))
#------------------------------------------------------------
# reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
# print(find_overlap(reads, 4))
# print(find_overlap(reads, 5))
#------------------------------------------------------------
reads, _ = readFastq("data/ERR266411_1.for_asm.fastq")
tuples = find_overlap(reads, 30)
print(len(tuples))

unique_firsts = set(tup[0] for tup in tuples)
count = 0
for entry in unique_firsts:
    if entry in reads:
        count +=1
print(count)