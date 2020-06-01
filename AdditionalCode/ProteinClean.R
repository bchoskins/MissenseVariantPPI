#Protein Data clean up

##############Interactions##################
inBio = read.table("core.txt", header=FALSE, sep='\t')

inBioUnique <- inBio[!duplicated(inBio$V1),]
inBioUnique2 <- inBio[!duplicated(inBio$V2),]

df <- inBio[-c(7:14,16)]

colnames(df) <- c("Protein_A", 'Protein_B',  "AlternateID_A", 'AlternateID_B', "Alias_A", "Alias_B", 'Confidence_Score')

df[] <- lapply(df, function(x) gsub("[|].*", "", x))

df[] <- lapply(df, function(x) gsub(".*[:]", "", x))

df[] <- lapply(df, function(x) gsub("[(].*", "", x))

write.csv(df, file = "Interactions_InBio_Map.csv", row.names = FALSE)

####################Sequences################
seq <- read.table("uniprot_Sequences.txt", sep="\t", header=FALSE,
                  na.strings=".", stringsAsFactors=FALSE,
                  quote="", fill=FALSE)


data <- data.frame("Protein" = (1:74034), "Sequence" = (1:74034), stringsAsFactors = FALSE)

library(tidyverse)
data$Protein <- seq %>% filter(str_detect(V1, ">")) %>% unlist() %>% as.vector()

data$Sequence <- seq %>% filter(str_detect(V1,">", negate = TRUE)) %>% unlist() %>% as.vector()


step1 <- separate(data, Protein, into = c("Drop1", "Protein", "Rest"), sep = "\\|")

step2 <- separate(step1, "Rest", into = c("AlternateID", "rest"), sep = "\\s", extra = "merge")

step3 <- separate(step2, rest, into = c("Drop2", "Drop3", "Drop4", "Alias", "Drop5", "Drop6"), sep = "\\=")

step4 <- separate(step3, Alias, into = c("Alias", "Drop7"), sep = "\\s")

sequenceData <- step4[-c(1,4:6,8:10)]

write.csv(sequenceData, file="Uniprot_SequencesGOOD.csv", row.names = FALSE)

###############Combine sequences to Protein_A and Protein_B in interactions data#############################3

sequences <- read.csv(file="Uniprot_SequencesGOOD.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

interactions <- read.csv(file="InBio_Map_core_2016_09_12/Interactions_InBio_Map.csv", header = TRUE, sep=",", stringsAsFactors = FALSE)

intUniqueA <- interactions[!duplicated(interactions$Protein_A),]
intUniqueB <- interactions[!duplicated(interactions$Protein_B),]

allInt <- dplyr::union(intUniqueA, intUniqueB)

allIntUnique <- allInt[!duplicated(allInt),]
#306262

sequences <- sequences[-c(2,3)]

interactions <- dplyr::left_join(interactions, sequences, by = c("Protein_A" = "Protein"))

colnames(interactions)[8] <- "Sequence_A"

interactions <- dplyr::left_join(interactions, sequences, by = c("Protein_B" = "Protein"))

colnames(interactions)[9] <- "Sequence_B"

complete <- na.omit(interactions)

write.csv(complete, file = "UniprotModelData.csv", row.names = FALSE)

#########################Create negatives#######################
#shuffle the interactions data and then see if any pairs show up in the original
#if not keep and set their conf. score to 0

neg <- read.csv(file="UniprotModelData.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

justAProt <- neg[-c(4,6,7,9)]

shuffle <- transform(justAProt, Shuffled = sample(Protein_B) )

shuffle <- shuffle[-c(2)]

colnames(shuffle)[5] <- c("Protein_B")

#removes rows that have interaction from original data
shuffle <- shuffle[c(1,5,2:4)]
test <- dplyr::setdiff(shuffle, justAProt)

B <- neg[-c(1,3,5,7,8)]

#so I can add Protein B back in matched by Protein name
Bunique <- B[!duplicated(B$Protein_B),]

#check total num of proteins
# Aunique <- justAProt[!duplicated(justAProt$Protein_A),]
# 
# totalUnique <- neg[!duplicated(neg$Protein_A),]

add1 <- dplyr::left_join(test, Bunique, by = "Protein_B")

add1$Confidence_Score <- 0

add1 <- add1[c(1,2,3,6,4,7,5,8,9)]

##NEW: subsample these negatives since we will get other later 
add1 <- sample_n(add1, 200000)

check <- dplyr::intersect(add1, neg)

fullneg <- dplyr::bind_rows(neg, add1)

set.seed(42)

rows <- sample(nrow(fullneg))

fullneg <- fullneg[rows, ]

fullneg <- fullneg[c(1:6,8,9,7)]

#write.csv(fullneg, file="UniprotModelDataWithNegatives.csv")

#added binary column 
fullneg$Binary <- 0

fullneg$Binary <- ifelse(fullneg$Confidence_Score>=0.5, 1, 0)
#negatives = 653946, positives = 158226

rownames(fullneg) <- NULL

write.csv(fullneg, file="UniprotModelDataBinary.csv")

checkUni <- read.csv(file="UniprotModelData.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

check2 <- dplyr::intersect(fullneg, checkUni)



########Need to get all unique proteins and their respective sequence to encode the seqeuences
data = read.csv(file="UniprotModelDataBinary.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
write.csv(data, file="UniprotModelDataBinary.csv")

check3 <- dplyr::setdiff(fullneg, data)

data <- data[-c(1)]

protein <- data[-c(3:6,9,10)]

protein <- protein[c(1,3,2,4)]

proteinA <- protein[c(1,2)]
proteinB <- protein[c(3,4)]

colnames(proteinB) <- c("Protein_A", "Sequence_A")

allProt <- dplyr::union(proteinA, proteinB)
allProtUnique <- allProt[!duplicated(allProt$Protein_A),]

n <- dim(allProtUnique)[1]
allProtUnique <- allProtUnique[1:(n-1),]

write.csv(allProtUnique, file="AllProteins.csv")



