def Complement(DNA):
    RNA = []
    DNAbases = ['A','G','C','T']
    RNAbases = ['T','C','G','A']
    for n in DNA:
        for i in range(4):
            if n == DNAbases[i]:
                RNA.append(RNAbases[i])
    return RNA
def Protiens(sequence):
    starts = []
    stops = []
    index = 0
    sequence = list(sequence)

    while '*' in sequence[index:len(sequence)]:
        index = sequence.index('*', index)
        starts.append(index)
        subends = []
        for i in ['M']:
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
        if stops[i] - starts[i] >= 100:
            print("exon: ",sequence[starts[i]:stops[i]+3])
        else:
            print("intron: ",sequence[starts[i]:stops[i]+3])

file = open('chr12.fa')
data = file.read()
data = data.upper()
file.close()
ALDH2 = 'MLRAAARFGPRLGRRLLSAAATQAVPAPNQQPEVFCNQIFINNEWHDAVSRKTFPTVNPSTGEVICQVAEGDKEDVDKAVKAARAAFQLGSPWRRMDASHRGRLLNRLADLIERDRTYLAALETLDNGKPYVISYLVDLDMVLKCLRYYAGWADKYHGKTIPIDGDFFSYTRHEPVGVCGQIIPWNFPLLMQAWKLGPALATGNVVVMKVAEQTPLTALYVANLIKEAGFPPGVVNIVPGFGPTAGAAIASHEDVDKVAFTGSTEIGRVIQVAAGSSNLKRVTLELGGKSPNIIMSDADMDWAVEQAHFALFFNQGQCCCAGSRTFVQEDIYDEFVERSVARAKSRVVGNPFDSKTEQGPQVDETQFKKILGYINTGKQEGAKLLCGGGIAADRGYFIQPTVFGDVQDGMTIAKEEIFGPVMQILKFKTIEEVVGRANNSTYGLAAAVFTKDLDKANYLSQALQAGTVWVNCYDVFGAQSPFGGYKMSGSGRELGEYGLQAYTEVKTVTVKVPQKNS'
print("Has Asian Flush: ",Complement(data[111803962]) == 'G')
Protiens(ALDH2)
