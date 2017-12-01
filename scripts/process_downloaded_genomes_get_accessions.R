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

# location of downloaded sequences
downloaded_genomes_location <- "downloaded_refseq_genomes"

# create object for storing all sequences
all_sequences <- DNAStringSet()
all_plasmid_sequences <- DNAStringSet()

print("reading_sequences")
# read all sequences from specified folder and append them to all_sequences
# object
for (file in list.files(downloaded_genomes_location, "genomic.fna.gz$", full.names = T)) {
    print(file)
    sequences <- readDNAStringSet(file)
    all_sequences <- c(all_sequences, sequences)
}

sequence_names <- names(all_sequences)

# get accessions of sequences
accessions1_matches <- regexpr("^(.*?)\\.[0-9]", sequence_names)
accessions1 <- regmatches(sequence_names, accessions1_matches)
accessions <- sub("^(ref\\||)(.*)\\.[0-9]", "\\2", accessions1)


# read all sequences of plasmids and append them to all_sequences object
for (file in list.files(".", "plasmid.*genomic.fna.gz$", full.names = T)) {
    print(file)
    sequences <- readDNAStringSet(file)
    all_plasmid_sequences <- c(all_plasmid_sequences, sequences)
}

sequence_names_plasmids <- names(all_plasmid_sequences)

# get accessions of plasmids
accessions1 <- sub("gi(.*?)[:blank:].*", "\\1", sequence_names_plasmids)
accessions2 <- sub(".*ref(.*).*", "\\1", accessions1)
accessions3 <- sub(".(.*[0-9])\\|.*", "\\1", accessions2)
accessions_plasmids <- sub("(.*)..", "\\1", accessions3)

# merge genomic and plasmid sequences
all_sequences <- c(all_sequences, all_plasmid_sequences)
# merge genomic and plasmid accessions
accessions <- c(accessions, accessions_plasmids)

# remove duplicated accessions
check_dup_accessions <- !duplicated(accessions)

# create deduplicated sequences dataset
accessions_dedup <- accessions[check_dup_accessions]
all_sequences_dedup <- all_sequences[check_dup_accessions]

# save imported sequences to use them in further steps
save.image("imported_sequences.RData")

# write a file with accessions (will be used for filtering RefSeq catalog)
write.table(accessions_dedup, file = "accessions.tsv", sep = "\t", col.names = F,
    row.names = F, quote = F)
