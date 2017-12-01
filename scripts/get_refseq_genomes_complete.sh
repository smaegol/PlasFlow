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
refseq_bacteria_assembly_summary_link='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
refseq_archaea_assembly_summary_link='ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'
assembly_summary_file_bacteria='assembly_summary_bacteria.txt'
assembly_summary_file_archaea='assembly_summary_archaea.txt'
assembly_summary_file='assembly_summary_merged.txt'
assembly_summary_complete_genomes='assembly_summary_merged_complete_genomes.txt'
assembly_summary_complete_genomes_filtered='assembly_summary_merged_complete_genomes_filtered.txt'
uniq_species_file='refseq_uniq_species_complete_genomes.txt'
download_links_file='refseq_uniq_species_complete_genomes_download_links.txt'
output_folder='downloaded_refseq_genomes'

#other variables
complete_genome_mark='Complete Genome'

#download required files
wget -O $assembly_summary_file_bacteria -nc $refseq_bacteria_assembly_summary_link
wget -O $assembly_summary_file_archaea -nc $refseq_archaea_assembly_summary_link
#merge bacteria and archaea
cp $assembly_summary_file_bacteria $assembly_summary_file
cat $assembly_summary_file_archaea >> $assembly_summary_file
#get complete genomes
grep "$complete_genome_mark" $assembly_summary_file > $assembly_summary_complete_genomes
#filter genomes excluded from Refseq (based on 21th column in assembly_summary.txt
awk -F "\t" '$21=="" {print $0}' $assembly_summary_complete_genomes > $assembly_summary_complete_genomes_filtered
#get uniq species
cut -f7 $assembly_summary_complete_genomes_filtered | uniq | sort -n > $uniq_species_file

#get download links for each species (prefer reference or representative genomes)
/usr/bin/Rscript get_refseq_genomes_complete_get_download_links.R

#create output folder (if needed)
if [ -d "$output_folder" ] ; then
  echo "output folder $output_folder exisits!!"
else
  mkdir $output_folder
fi


#download all sequences
for next in $(cat $download_links_file); do
  wget -a download_log.txt -P $output_folder -nc "$next"/*genomic.fna.gz;
done

for next in $(cat $download_links_file); do wget -a download_log_gbff.txt -P $output_folder -nc "$next"/*genomic.gbff.gz; done

exit
