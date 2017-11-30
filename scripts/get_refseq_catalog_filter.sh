#!/bin/bash


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


#definition of files
refseq_catalog_link='ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/RefSeq-release*.catalog.gz'
refseq_catalog_file='Refseq-release82.catalog.gz'
refseq_gunzipped_catalog_file='RefSeq-release82.catalog'
filtered_catalog_file='Refseq-release82.catalog.filtered'
output_folder='downloaded_refseq_genomes'
accessions_file='accessions.tsv'


#paths
filter_large_file_path='/home/smaegol/storage/analyses/plasmids/final_modelling4/new_bacteria/filter_large_file.pl'


#download RefSeq catalog if needed
if [ -e "$refseq_gunzipped_catalog_file" ]; then
	echo "Catalog file exists..."
else
	echo "not";
	wget -O $refseq_catalog_file -nc $refseq_catalog_link
	gunzip $refseq_catalog_file
fi

#read downloaded genomic sequences to get accession numbers
/usr/bin/Rscript process_downloaded_genomes_get_accessions.R

#filter RefSeq by obtained accessions
$filter_large_file_path $accessions_file $refseq_gunzipped_catalog_file $filtered_catalog_file

#create annotation table of downloaded sequences (get taxonomy info)
