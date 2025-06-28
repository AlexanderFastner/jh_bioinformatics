from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
#------------------------------------------------------------
#(1) How many records are in the file?
records = list(SeqIO.parse("dna2.fasta", "fasta"))
# for record in records:
#     print(record.id)
print("Number of records in Fasta", len(records))
# print(records[6])
print("------------------------------------------------------------")

#------------------------------------------------------------
#(2) What are the lengths of the sequences in the file? 
#What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers?
lengths_records = {}
for record in records:
    lengths_records[record.id] = len(record.seq)
    # print(len(record.seq))
print("Lengths of all sequences: ", lengths_records)

# lengths_records = {"a": 1,"b": 1, "t": 6, "u": 7, "k": 7, "o": 2}
max_value = max(lengths_records.values())
max_items = [(key, value) for key, value in lengths_records.items() if value == max_value]
min_value = min(lengths_records.values())
min_items = [(key, value) for key, value in lengths_records.items() if value == min_value]

if len(max_items) > 1:
    print("Multiple Longest:")
    print(max_items)
else:
    print("Longest sequence is:", max_items)
if len(min_items) > 1:
    print("Multiple Shortest:")
    print(min_items)
else:
    print("Shortest sequence is:", min_items)

print("------------------------------------------------------------")
#------------------------------------------------------------
#3 identify all ORFs present in each sequence of the FASTA file
#What is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? 
#For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it? 

start = "ATG"
stop = ["TAA", "TAG", "TGA"]
all_ORFs = {}
longest_orf = None

#!!!only positives so far * need to reverse commplement first before input into this method
def find_ORFS(seq):
    # return tuple of reading Frame and the ORF itself
    i = 0
    while i < 3:
        # print("i:", i)
        # print(seq[i:])
        working_seq = seq[i:]
        ORFs = []
        starting = findstart(working_seq, 0)
        j = 0
        
        while starting < len(seq):
        #while either starting or stopping position is not past len(working_seq) keep going
            if starting != -1:
                #find stop
                stopping = findstop(working_seq, starting)
                if stopping != -1:
                    #Have both start and stop
                    #make and save orf
                    ORF = (i, working_seq[starting: stopping], starting)
                    # print(ORF)
                    ORFs.append(ORF)
                    #look for next start and restart the loop
                    starting = findstart(working_seq, stopping)
                    # print("NEXT START", starting)

                else:
                    # print("NO STOP") 
                    break
            else:
                # print("NO START")
                break
        i+=1
    return ORFs

#given a sequence and a start position find the next start codon
def findstart(seq, starting):
    offset = starting
    if seq[starting:].find(start):
        # print("START FOUND:", seq[starting:].find(start)+ offset)
        return seq[starting:].find(start)+ offset
    else: 
        return -1

#given a sequence and a start position find the next stop codon
def findstop(seq, starting):
    for j in range(starting + 3, len(seq), 3):
        if seq[j:j+3] in stop:
            # print("STOP FOUND:", j, seq[j:j+3])
            return j + 3
        else:
            continue
    return -1


for record in records:
    #check all 3 ORFs for start AND stop codons then return ORFs for each
    all_ORFs[record.id] = find_ORFS(record.seq)
    print(record.id, "|", len(all_ORFs[record.id]))
    
print("------------------------------------------------------------")
longest_orf = ("None", 0, 0)
for id, seq in all_ORFs.items():
    for orfs in seq:
        if len(orfs[1]) > longest_orf[1]:
            longest_orf = (id, len(orfs[1]), orfs[2])
            
#What is the length of the longest ORF in the file?
print("Length of the longest ORF in the file:", longest_orf[1])
#What is the identifier of the sequence containing the longest ORF? 
print("ID for longest ORF in the file:", longest_orf[0])
#For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it?

def longest_Orf_in_id(id):
    longest_orf = ("None", 0, 0) 
    for orf in all_ORFs[id]:
        # print(orf, len(orf[1]))
        if len(orf[1]) > longest_orf[1]:
            longest_orf = (id, len(orf[1]), orf[2])
    print("Longest Orf in" , longest_orf[0], "is", longest_orf[1], "starts at", longest_orf[2])
longest_Orf_in_id("gi|142022655|gb|EQ086233.1|16")
print("------------------------------------------------------------")

#------------------------------------------------------------
#4 Given a length n, your program should be able to identify all repeats of length n in all sequences in the FASTA file. 
#Your program should also determine how many times each repeat occurs in the file, and which is the most frequent repeat of a given length.
import itertools
def find_repeats_of_length_n(length):
    #returns a list of all repeats of length n that occur at least 1 time

    #make all possible subset combinations to try against
    #given ACTG make a list of all possible products of length n
    c = list(itertools.product("ACTG", repeat= length))
    combinations = {}
    for combination in c:
        combinations[("".join(combination))] = 0

    for record in records:
        for c in combinations.keys():
            combinations[c] +=record.seq.count_overlap(c)
        filtered = {key : value for key, value in combinations.items() if value > 1}
        # print(filtered)

    return filtered


repeats = find_repeats_of_length_n(3)
print(repeats)
max_item = max(repeats.items(), key=lambda item: item[1])
print("Most Frequent repeat:", max_item)
print("------------------------------------------------------------")