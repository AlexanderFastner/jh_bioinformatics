from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
# print(dir(Bio))
print("-------")

# sequence_data = open("blast.fasta").read() 
# result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)

# # blast_record = NCBIXML.read(result_handle)
# # print(len(blast_record.alignments))

# with open('results.xml', 'w') as save_file: 
#     blast_results = result_handle.read() 
#     save_file.write(blast_results)


# E_VALUE_THRESH = 1e-20 
# for record in NCBIXML.parse(open("results.xml")): 
#     if record.alignments: 
#        print("\n") 
#        print("query: %s" % record.query[:100]) 
#        for align in record.alignments: 
#           for hsp in align.hsps: 
#              if hsp.expect < E_VALUE_THRESH: 
#                 print("match: %s " % align.title[:100])


# seq = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")
# print(seq.translate())

