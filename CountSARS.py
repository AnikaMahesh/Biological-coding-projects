
file = open('SARS-2020.fasta')
reading = file.read()
length = 0
for i in list(reading):
   if i in ["A","C","G","T","U"]:
      length = length+1
file.close()
print(length)

