#this section will contain all the necessary packages that are loaded

library(msa)
library(Biostrings)
library(seqinr)
library(phangorn)
library(tidyr)
library(dplyr)
library(UniprotR)

#aligning the .txt file 
#naming it for convenience

humanalign <- readDNAStringSet("/Users/kamilkawiecki/Downloads/sequences.fasta")

#msa alignment 

hmsaalign<- msa(humanalign)
nchar(hmsaalign)

#this allows me to see the number of nucleotides in each sequence
#sequence 6 has an A nucleotide deletion

print(hmsaalign,show="complete")

#further investigation of differences by printing sequences
#sequence 4 has a C->A point mutation
#sequence 6 has an A deletion (position 1) as well as two point mutations
#one of which is A -> G (postion 47) and the other A->T at about 140
#sequence 10 has 2 point mutations (one is A->T)

alfreq <- alphabetFrequency(hmsaalign)

  


#3 placed sequences into ncbi BLAST and the results is a gene
# for providing the instruction for beta globin which is a component 
# of hemoglobin (2 alpha and 2 beta globins)
# LC121775.1


#translating the sequence to protein first means you have to seperate it from
#the rest of the sequence

class(humanalign)
homo6 <- humanalign$Homo_sapiens_6
print(homo6)

#Biostrings translates this to an amino acid 

aaseq6 <- Biostrings::translate(homo6)
print(aaseq6)

#write AAsequence to a blank file as a FASTA 

write.fasta(Biostrings::translate(homo6),"Homo_sapiens_6","/Users/kamilkawiecki/Documents/Bioinformatics/FASTA",open = "w",nbchar = 60,as.string = FALSE)

#5 A0A0J9YWK4



#6 UniProt results indicate beta Thalassemia
#checked 4 and 10 for good measure but they were related closer to gorillas


