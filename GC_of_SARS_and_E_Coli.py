from Bio.SeqUtils import GC
from Bio.Seq import Seq
def rna(DNA):
    RNA = []
    DNAbases = ['A','G','C','T']
    RNAbases = ['U','C','G','A']
    for n in DNA:
        for i in range(4):
            if n == DNAbases[i]:
                RNA.append(RNAbases[i])
    return ''.join(RNA)
file = open("SARS-2020.fasta")
Data = Seq(rna(file.read()))
file.close()
file2 = open("small-genome")
Data2 = Seq(rna(file2.read()))
file2.close()
print("GC of E-Coli: ",GC(Data))
print("GC of SARS Virus: ",GC(Data2))


