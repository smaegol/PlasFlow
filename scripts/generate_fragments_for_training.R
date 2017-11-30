#######################################################################################
###                                                                                 ###
###     Copyright (C) 2017  Pawel Krawczyk (p.krawczyk@ibb.waw.pl)                  ###
###                                                                                 ###
###     This program is free software: you can redistribute it and/or modify        ###
###     it under the terms of the GNU General Public License as published by        ###
###     the Free Software Foundation, either version 3 of the License, or           ###
###     (at your option) any later version.                                         ###
###                                                                                 ###
###     This program is distributed in the hope that it will be useful,             ###
###     but WITHOUT ANY WARRANTY; without even the implied warranty of              ###
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               ###
###     GNU General Public License for more details.                                ###
###                                                                                 ###
###     You should have received a copy of the GNU General Public License           ###
###     along with this program. If not, see <http://www.gnu.org/licenses/>.        ###
###                                                                                 ###
#######################################################################################

library(Biostrings)
library(data.table)
library(plyr)
library(ggplot2)
library(XML)
library(rentrez)


# define splitting parameters
SAMPLE_LENGTH <- 10000
FRACTION_COVERED_PLASMIDS <- 0.75
FRACTION_COVERED_CHROMOSOMES <- 0.05

make_sample <- function(x, sample_length, fraction_covered) {
  #function generating random sequence fragments of given length
    splitted_set <- DNAStringSet()
    x_length <- width(x)
    print(paste("Sequence of length:", x_length))
    if (x_length <= SAMPLE_LENGTH) {
        splitted_set <- c(splitted_set, x)
    } else {
        no_samples <- ceiling(x_length/sample_length * fraction_covered)
        print(paste("dividing into", no_samples, "samples"))
        sequence_starts <- sample(x_length - SAMPLE_LENGTH - 1, no_samples)
        for (z in sequence_starts) {
            temp_seq <- subseq(x, z, z + SAMPLE_LENGTH - 1)
            splitted_set <- c(splitted_set, temp_seq)
        }
    }
    return(splitted_set)
}



generate_class2 <- function(x) {
    # function generating classes used for training, based on origin (plasmid/chromosome) and taxonomy (phylum)
    # to be used with apply
    element1 <- as.character(x[9])
    plasmid <- x[5]
    element1 <- gsub("[[:space:]]", "", element1)
    element1 <- gsub(",", "", element1)
    element1 <- gsub("-", "", element1)
    if (!is.na(element1)) {
        if (element1 == "") {
            element1 <- "other"
        }
        if (element1 == "undef") {
            element1 <- "other"
        }
    } else {
        element1 <- "other"

    }
    final_class = paste(plasmid, element1, sep = ".")
    print(final_class)
    return(final_class)
}

generate_class3 <- function(x) {
    # classes tax rank 7 - class 4 - plasm
    element1 <- as.character(x[10])
    plasmid <- x[5]
    element1 <- gsub("[[:space:]]", "", element1)
    element1 <- gsub(",", "", element1)
    element1 <- gsub("-", "", element1)
    if (!is.na(element1)) {
        if (element1 == "") {
            element1 <- "other"
        }
        if (element1 == "undef") {
            element1 <- "other"
        }
    } else {
        element1 <- "other"

    }
    final_class = paste(plasmid, element1, sep = ".")
    print(final_class)
    return(final_class)
}


filter_class2 <- function(x) {
    # functions used to filter low count classes

    min_class_count <- 100 #minimum number of sequences to include in class
    class_get <- as.character(x)
    class_return = class_get
    if (levels_count[levels_count$factor == class_get, ]$count < min_class_count) {
        if (grepl("plasmid", class_get)) {
            class_return = "plasmid.other"
        }
        if (grepl("chromosom", class_get)) {
            class_return = "chromosome.other"
        }
    }
    return(class_return)
}




# read previously imported sequences
load("imported_sequences.RData")

#load taxonomy annotation data
# should be manually edited to correct wrong annotations and exclude unwanted sequences
annotation_data<-fread("refseq_annotation.tsv",data.table=F)


accessions_df <- data.frame(accessions = accessions_dedup)
# get_annotation_of_all_sequences in the next step include Archaea
all_seqs_acc_annotated <- merge(annotation_data, accessions_df,
    by.x = "accession", by.y = "accessions")


accessions_mapping <- data.frame()

#DNAStringSet for storing splitted sequences:
all_sequences_split <- DNAStringSet()
#vector storing accessions of splitted fragments
all_sequences_split_acc <- c()


#Perform actual random splitting:
for (y in seq(1, length(all_sequences_dedup))) {
    temp_seq <- all_sequences_dedup[y]
    temp_acc <- accessions_dedup[y]
    if (nrow(all_seqs_acc_annotated[all_seqs_acc_annotated$accession == temp_acc,
        ]) > 0) {
          #check if given sequence is plasmid
        isplasmid <- all_seqs_acc_annotated[all_seqs_acc_annotated$accession == temp_acc,
            ]$plasmid
        if (isplasmid == "chromosome") {
            print(isplasmid)
            splitted_set <- make_sample(temp_seq, SAMPLE_LENGTH, FRACTION_COVERED_CHROMOSOMES)
        } else {
            splitted_set <- make_sample(temp_seq, SAMPLE_LENGTH, FRACTION_COVERED_PLASMIDS)
        }
    }
    print(paste("processed", y, "out of", length(all_sequences_dedup), "sequences"))
    all_sequences_split <- c(all_sequences_split, splitted_set)
    no_splitted_seqs <- length(splitted_set)
    accessions_splitted <- rep(temp_acc, no_splitted_seqs)
    all_sequences_split_acc <- c(all_sequences_split_acc, accessions_splitted)
}

#Calculate kmer counts (kmers from 3 to 7):
kmer3<-trinucleotideFrequency(all_sequences_split)
kmer3_frame<-as.data.frame(kmer3)
kmer3_frame$accessions<-all_sequences_split_acc

kmer4<-oligonucleotideFrequency(all_sequences_split,4)
kmer4_frame<-as.data.frame(kmer4)
kmer4_frame$accessions<-all_sequences_split_acc

kmer5<-oligonucleotideFrequency(all_sequences_split,5)
kmer5_frame<-as.data.frame(kmer5)
kmer5_frame$accessions<-all_sequences_split_acc

kmer6<-oligonucleotideFrequency(all_sequences_split,6)
kmer6_frame<-as.data.frame(kmer6)
kmer6_frame$accessions<-all_sequences_split_acc

kmer7<-oligonucleotideFrequency(all_sequences_split,7)
kmer7_frame<-as.data.frame(kmer7)
kmer7_frame$accessions<-all_sequences_split_acc

#annotate kmer count tables:
kmer3_annotated<-merge(annotation_tax_edited_no_Archaea,kmer3_frame,by.x="accession",by.y="accessions")
kmer4_annotated<-merge(annotation_tax_edited_no_Archaea,kmer4_frame,by.x="accession",by.y="accessions")
kmer5_annotated<-merge(annotation_tax_edited_no_Archaea,kmer5_frame,by.x="accession",by.y="accessions")
kmer6_annotated<-merge(annotation_tax_edited_no_Archaea,kmer6_frame,by.x="accession",by.y="accessions")
kmer7_annotated<-merge(annotation_tax_edited_no_Archaea,kmer7_frame,by.x="accession",by.y="accessions")

#generate class (in form plasmid.phylum) for each sequence
annotation_temp<-kmer3_annotated[,c(1:14)]
tax_classes<-apply(annotation,1,function(x) generate_class2(x))

#get counts for each class
levels_count<-as.data.frame(table(tax_classes))
colnames(levels_count) = c('factor','count')

#filter classes by count
tax_classes_filtered<-sapply(tax_classes,function(x) filter_class2(x))

#convert class names to numeric values (which will be provided as an input to TensorFlow)
tax_classes_numeric<-as.numeric(as.factor(tax_classes_filtered))-1
#create class labels data frame - will be required for prediction
class_labels<-levels(as.factor(tax_classes_filtered))
labels_df_num <- as.numeric(as.factor(class_labels)) - 1
labels_df <- data.frame(id = labels_df_num, label = class_labels)
# write class labels to file:
write.table(labels_df,file="class_labels_phyla_df.tsv",sep="\t",row.names=F)

#add generated classes to kmer counts tables and remove other annotation data:
kmer3_annotated_tax_classes<-cbind(tax_classes_numeric,kmer3_annotated[,-c(1:14)])
colnames(kmer3_annotated_tax_classes)[1]<-'plasmid'

kmer4_annotated_tax_classes<-cbind(tax_classes_numeric,kmer4_annotated[,-c(1:14)])
colnames(kmer4_annotated_tax_classes)[1]<-'plasmid'

kmer5_annotated_tax_classes<-cbind(tax_classes_numeric,kmer5_annotated[,-c(1:14)])
colnames(kmer5_annotated_tax_classes)[1]<-'plasmid'

kmer6_annotated_tax_classes<-cbind(tax_classes_numeric,kmer6_annotated[,-c(1:14)])
colnames(kmer6_annotated_tax_classes)[1]<-'plasmid'

kmer7_annotated_tax_classes<-cbind(tax_classes_numeric,kmer7_annotated[,-c(1:14)])
colnames(kmer7_annotated_tax_classes)[1]<-'plasmid'


#Write kmer count data to files:
write.table(kmer3_annotated_tax_classes,file="kmer3_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv",sep="\t",row.names=F,quote=F)
write.table(kmer4_annotated_tax_classes,file="kmer4_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv",sep="\t",row.names=F,quote=F)
write.table(kmer5_annotated_tax_classes,file="kmer5_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv",sep="\t",row.names=F,quote=F)
write.table(kmer6_annotated_tax_classes,file="kmer6_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv",sep="\t",row.names=F,quote=F)
write.table(kmer7_annotated_tax_classes,file="kmer7_split_raw_counts_tax_classes_numeric_annotated_filtered_tax.tsv",sep="\t",row.names=F,quote=F)
