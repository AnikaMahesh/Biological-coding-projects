import operator
#showing similarity of a gene in different organisms
file = open('animals.txt')
data = file.readlines()
file.close()
center = 'Whale'
sequences = []
differences = []
priviousi = 0
for i in range(int(len(data)/2)):
    if center in data[i*2]:
        center = data[i*2 + 1]
    sequences.append((data[i*2],data[i*2 + 1]))
for name,sequence in sequences:
    similarity = 0
    sequence = ''.join(sequence)
    center = ''.join(center)       
    for i in range(min(len(''.join(center)),len(''.join(sequence)))):
        if ''.join(center)[i] == ''.join(sequence)[i]:
            similarity = similarity+1
    differences.append(similarity)
sequences = zip(differences,sequences)
sequences = sorted(sequences, key=operator.itemgetter(0))
print(sequences)
