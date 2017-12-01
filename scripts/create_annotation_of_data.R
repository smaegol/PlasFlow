
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
library(XML)
library(rentrez)

get_taxonomy_by_taxid <- function(x) {
    # download taxonomy information from NCBI using rentrez

    taxon_id <- x
    taxonomy <- entrez_fetch(db = "taxonomy", id = taxon_id, rettype = "XML")
    taxonomy_xml <- xmlParse(taxonomy)
    taxonomy_df <- xmlToDataFrame(getNodeSet(taxonomy_xml, "//TaxaSet/Taxon/LineageEx/Taxon"))

    taxon_phylum <- as.character(taxonomy_df[taxonomy_df$Rank == "phylum", ]$ScientificName)
    taxon_class <- as.character(taxonomy_df[taxonomy_df$Rank == "class", ]$ScientificName)
    taxon_order <- as.character(taxonomy_df[taxonomy_df$Rank == "order", ]$ScientificName)
    taxon_family <- as.character(taxonomy_df[taxonomy_df$Rank == "family", ]$ScientificName)
    taxon_genus <- as.character(taxonomy_df[taxonomy_df$Rank == "genus", ]$ScientificName)

    taxon_species <- as.character(taxonomy_df[taxonomy_df$Rank == "species", ]$ScientificName)
    if (length(taxon_phylum) == 0) {
        taxon_phylum <- NA
    }
    if (length(taxon_class) == 0) {
        taxon_class <- NA
    }
    if (length(taxon_order) == 0) {
        taxon_order <- NA
    }
    if (length(taxon_family) == 0) {
        taxon_family <- NA
    }
    if (length(taxon_genus) == 0) {
        taxon_genus <- NA
    }
    if (length(taxon_species) == 0) {
        scientific_name <- sub("^.*?<ScientificName>(.*?)</ScientificName>.*", "\\1",
            taxonomy)
        taxon_species <- scientific_name
        # taxon_species<-NA
    }

    print(taxon_phylum)

    return_values <- data.frame(phylum = taxon_phylum, class = taxon_class, order = taxon_order,
        family = taxon_family, genus = taxon_genus, species = taxon_species, stringsAsFactors = FALSE)
    rownames(return_values) <- taxon_id

    return(return_values)

}

# read filtered Refseq catalog
annotation_refseq_filtered <- fread("RefSeq-release82.catalog.filtered", data.table = F)
colnames(annotation_refseq_filtered) <- c("tax_id", "species", "accession", "gi",
    "collection", "refseq_status", "length")

# get uniq taxon ids from Refseq catalog
taxon_ids <- annotation_refseq_filtered$tax_id
taxon_ids <- unique(taxon_ids)

# create refseq taxonomy table
refseq_taxonomy <- data.frame(phylum = character(), class = character(), order = character(),
    family = character(), genus = character(), species = character())

# get taxonomy information using rentrez
i = 0
for (taxon_id in taxon_ids) {
    i = i + 1
    print(i)
    if (length(grep(paste("^", taxon_id, "$", sep = ""), row.names(refseq_taxonomy))) ==
        0) {
        result_taxonomy <- get_taxonomy_by_taxid(taxon_id)
        refseq_taxonomy <- rbind(refseq_taxonomy, result_taxonomy)
    }
}

refseq_taxonomy$tax_id<-rownames(refseq_taxonomy)

taxonomy<-merge(annotation_refseq_filtered,refseq_taxonomy,by.x="tax_id",by.y="tax_id",all=T)

# write annotation data to file:
write.table(refseq_taxonomy, file = "refseq_annotation.tsv", sep = "\t", col.names = T,
    quote = T, row.names = F)

    # write annotation data to file:
    write.table(taxonomy, file = "refseq_annotation2.tsv", sep = "\t", col.names = T,
        quote = T, row.names = F)
