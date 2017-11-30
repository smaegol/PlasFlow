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

library(data.table)

# get complete reference or representative genomes

filtered_data <- fread("assembly_summary_merged_complete_genomes_filtered.txt", sep = "\t",
    data.table = F)
uniq_species <- fread("refseq_uniq_species_complete_genomes.txt", sep = "\t")

output_download_links <- c()
for (id in uniq_species$V1) {
    temp_data <- filtered_data[filtered_data$V7 == id, ]
    if (nrow(temp_data[grep("reference genome", temp_data$V5), ]) > 0) {
        output_download_links <- c(output_download_links, temp_data[grep("reference genome",
            temp_data$V5), ]$V20[1])
    } else if (nrow(temp_data[grep("representative genome", temp_data$V5), ]) > 0) {
        output_download_links <- c(output_download_links, temp_data[grep("representative genome",
            temp_data$V5), ]$V20[1])
    } else {
        output_download_links <- c(output_download_links, temp_data$V20[1])
    }
}

output_table <- as.data.frame(output_download_links)
write.table(output_table, file = "refseq_uniq_species_complete_genomes_download_links.txt",
    row.names = F, col.names = F, quote = F)
