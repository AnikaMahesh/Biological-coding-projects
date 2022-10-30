from Bio import SeqIO
# Six Frame Translation W biopython
def rna(DNA):
    RNA = []
    DNAbases = ['A','G','C','T']
    RNAbases = ['U','C','G','A']
    for n in DNA:
        for i in range(4):
            if n == DNAbases[i]:
                RNA.append(RNAbases[i])
    return RNA
def Protiens(sequence):
    starts = []
    stops = []
    index = 0
    sequence = split(sequence)
    while 'AUG' in sequence[index:len(sequence)]:
        index = sequence.index('AUG', index)
        starts.append(index)
        subends = []
        for i in ['UAG','UAA','UGA']:
            if i in sequence[index:len(sequence)]:
                subends.append(sequence.index(i,index))
            else:
                subends.append(len(sequence))
        index = min(subends)
        stops.append(index)
    if not len(starts) == len(stops):
        print("Incomplete Sequence")
        return False
    for i in range(len(starts)):
        protiens = [sequence[starts[i]:stops[i]+3] for i in range(len(starts)) if stops[i] - starts[i] >= 300]
    return protiens
    
def split(sequence):
    tmp = []
    for i in range(int(len(sequence)/3)):
        tmp.append(sequence[(i*3): (i*3)+3])
    sequence = tmp
    return sequence
def translate(segment):
    AminoAcids = []
    segment = ''.join(segment)
    for i in range(int(len(segment)/3)):
        if segment[3*i] == 'G':
            if segment[3*i+1] == 'G':
                AminoAcids.append('Gly')
            elif segment[3*i + 1] == 'U':
                AminoAcids.append('Val')
            elif segment[3*i + 1] == 'A':
                if segment[3*i + 2] == 'G' or segment[3*i + 2] == 'A':
                    AminoAcids.append('Glu')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('Asp')
            elif segment[3*i + 1] == 'C':
                AminoAcids.append('Ala')

        if segment[3*i] == 'U':
            if segment[3*i+1] == 'G':
                if segment[3*i + 2] == 'G':
                    AminoAcids.append('Trp')
                elif segment[3*i + 2] == 'A':
                    AminoAcids.append('STOP')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C': # 'UAG' or 'UAA' or 'UGA' or 'AUG'
                    AminoAcids.append('Asp')                
            elif segment[3*i + 1] == 'U':
                if segment[3*i + 2] == 'G' or segment[3*i + 2] == 'A':
                    AminoAcids.append('Leu')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('Phe')
            elif segment[3*i + 1] == 'A':
                if segment[3*i + 2] == 'G' or segment[3*i + 2] == 'A':
                    AminoAcids.append('STOP')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('Thr')
            elif segment[3*i + 1] == 'C':
                AminoAcids.append('Ser')

        if segment[3*i] == 'A':
            if segment[3*i+1] == 'G':
                if segment[3*i + 2] == 'A' or segment[3*i + 2] == 'G':
                    AminoAcids.append('Arg')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('Ser')                
            elif segment[3*i + 1] == 'U':
                if segment[3*i + 2] == 'G':
                    AminoAcids.append('START')
                else:
                    AminoAcids.append('Ile')
            elif segment[3*i + 1] == 'A':
                if segment[3*i + 2] == 'G' or segment[3*i + 2] == 'A':
                    AminoAcids.append('Lys')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('Asn')
            elif segment[3*i + 1] == 'C':
                AminoAcids.append('Thr') 

        if segment[3*i] == 'C':
            if segment[3*i+1] == 'G':
                AminoAcids.append('Arg')                
            elif segment[3*i + 1] == 'U':
                AminoAcids.append('Leu')
            elif segment[3*i + 1] == 'A':
                if segment[3*i + 2] == 'G' or segment[3*i + 2] == 'A':
                    AminoAcids.append('Gln')
                elif segment[3*i + 2] == 'U' or segment[3*i + 2] == 'C':
                    AminoAcids.append('His')
            elif segment[3*i + 1] == 'C':
                AminoAcids.append('Pro')
    return AminoAcids
def complement(DNA):
    DNAcomplement = []
    DNAbases = ['A','G','C','T']
    complementry = ['T','C','G','A']
    for n in DNA:
        for i in range(4):
            if n == DNAbases[i]:
                DNAcomplement.append(complementry[i])
    return DNAcomplement
file = open("SARS-2020.fasta")
DNA = list(file.read())
file.close()
for i in range(2):
    print("complement",i+1)
    for j in range(3):
        print(" Frame ",j+1,": ")
        LOP = Protiens(''.join(rna(DNA[j:len(DNA)])))
        for k in range(len(LOP)):
            print("     lenght: ",len(LOP[k]), " protien: ", translate(LOP[k]), "string:", ''.join(LOP[k]))
    DNA = complement(DNA)
