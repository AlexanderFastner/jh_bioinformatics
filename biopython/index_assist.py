from kmer_index import Index
from algo_1 import readGenome
from bm_preproc import SubseqIndex

from algo_1 import naive_2mm
# #------------------------------------------------------------
def indexed_2mm(p, t, index_obj):
    segmnent_length = int(round(len(p)/(2+1)))
    all_matches = set()
    for i in range(2+1):
        start  = i * segmnent_length
        end = min((i + 1)*segmnent_length, len(p))
        # print(start, end)
        # print(p[start:end])

        matches = index_obj.query(p[start:end])
        print(len(matches))

        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mm = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mm += 1
                    if mm > 2:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mm += 1
                    if mm > 2:
                        break
            if mm <=2:
                all_matches.add(m -start)
    return list(all_matches)
# #------------------------------------------------------------
def subseq_indexed_2mm(p, t, index_obj):
    segmnent_length = int(round(len(p)/(2+1)))
    num_index_hits = 0
    all_matches = set()
    for i in range(2+1):
        # print(p[i:len(p):3])

        matches = index_obj.query(p[i:len(p):3])
        num_index_hits += len(matches)
        # print(matches)

        start  = i * segmnent_length
        end = min((i + 1)*segmnent_length, len(p))
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            mm = 0
            for j in range(0, start):
                if not p[j] == t[m-start+j]:
                    mm += 1
                    if mm > 2:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mm += 1
                    if mm > 2:
                        break
            if mm <=2:
                all_matches.add(m -start)
    return list(all_matches), num_index_hits
# #------------------------------------------------------------
# p = "AACTTG"
# t = "CACTTAATTTG"
# print(indexed_2mm(p, t, Index(t, 2)))
# #------------------------------------------------------------
print("regular index")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
t  = readGenome('chr1.GRCh38.excerpt.fasta')
ind = Index(t, 8)
hits = indexed_2mm(p, t, ind)   
print(sorted(hits))
print(len(hits))
print()
hits = naive_2mm(p, t)
print(sorted(hits))
print(len(hits))
# #------------------------------------------------------------
# t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
# p = 'to-morrow and to-morrow '
# subseq_ind = SubseqIndex(t, 8, 3)
# # print(subseq_ind.index)
# occurrences, num_index_hits = subseq_indexed_2mm(p, t, subseq_ind)
# print()
# print(occurrences)
# print(num_index_hits)
# #------------------------------------------------------------
# t = open('1110.txt.utf-8').read()
# p = 'English measure backward'
# subseq_ind = SubseqIndex(t, 8, 3)
# occurrences, num_index_hits = subseq_indexed_2mm(p, t, subseq_ind)
# print()
# print(occurrences)
# print(num_index_hits)
# #------------------------------------------------------------


# #------------------------------------------------------------
print("SubseqIndex")
p = "GGCGCGGTGGCTCACGCCTGTAAT"
t  = readGenome('chr1.GRCh38.excerpt.fasta')
subseq_ind = SubseqIndex(t, 8, 3)
occurrences, num_index_hits = subseq_indexed_2mm(p, t, subseq_ind)
print()
print(occurrences)
print(len(occurrences))
print(num_index_hits)