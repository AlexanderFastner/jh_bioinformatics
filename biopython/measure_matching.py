from bm_preproc import BoyerMoore
from algo_1 import readGenome
#------------------------------------------------------------
#naive and bm that return
#1 num of character comparisons
#2 num of alignments tried
#------------------------------------------------------------
def naive_with_counts(p, t):
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        num_alignments +=1
        match = True
        for j in range(len(p)):  # loop over characters
            num_character_comparisons+=1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences, num_alignments, num_character_comparisons
#------------------------------------------------------------
# p = 'word'
# t = 'there would have been a time for such a word'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# #------------------------------------------------------------
# p = 'needle'
# t = 'needle need noodle needle'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
#------------------------------------------------------------
#------------------------------------------------------------
#TODO add counts
def boyer_moore_with_counts(p, p_bm, t):
    i = 0
    num_alignments = 0
    num_character_comparisons = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        num_alignments += 1
        for j in range(len(p)-1, -1, -1):
            num_character_comparisons += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, num_alignments, num_character_comparisons

#------------------------------------------------------------
# print("NAIVE with Counts")
# p = 'word'
# t = 'there would have been a time for such a word'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# #------------------------------------------------------------
# p = 'needle'
# t = 'needle need noodle needle'
# occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
# print(occurrences, num_alignments, num_character_comparisons)
# #------------------------------------------------------------
# print()
# print("BM with Counts")
# p = 'word'
# t = 'there would have been a time for such a word'
# lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
# # #------------------------------------------------------------
# p = 'needle'
# t = 'needle need noodle needle'
# p_bm = BoyerMoore(p, lowercase_alphabet)
# occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
# print(occurrences, num_alignments, num_character_comparisons)
#------------------------------------------------------------
print()
print("Actual with counts")

p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
t  = readGenome('chr1.GRCh38.excerpt.fasta')
occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
print(occurrences, num_alignments, num_character_comparisons)
#------------------------------------------------------------
p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
t  = readGenome('chr1.GRCh38.excerpt.fasta')
p_bm = BoyerMoore(p, "ACGT")
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)