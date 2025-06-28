from algo_1 import readGenome
#------------------------------------------------------------
#TODO return the edit distance of the match between P and T with the fewest edits.
#TODO return smalled value in bottom row, not the last 
def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = 0
    # for i in range(len(y)+1):
    #     D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    # print("     ", end = "")
    # for l in y:
    #     print(f"{l:2}", end = "  ")
    # print()
    # for row in D:
    #     for num in row:
    #         print(f"{num:2}", end="  ")
    #     print()
    return min(D[-1])
#------------------------------------------------------------
p = "GCGTATGC"
p = "GCGTATGC"
t = "TATTGGCTATACGGTT"
# print(editDistance(p, t))
# print()
#------------------------------------------------------------
genome = readGenome('data/chr1.GRCh38.excerpt.fasta')
p = "GCTGATCGATCGTACG"
print(editDistance(p, genome))
#------------------------------------------------------------
p = "GATTTACCAGATTGAG"
print(editDistance(p, genome))
#------------------------------------------------------------

