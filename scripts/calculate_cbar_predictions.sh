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


#define paths
cbar_path="/home/smaegol/storage/analyses/plasmids/cBar/cBar.1.2/cBar.pl"
#those scripts are expected to be located in the current folder
output_fasta_from_cbar_path="perl get_seqs_by_ids_list_cBar.pl "
filter_sequences_by_length_script='perl filter_sequences_by_length.pl'

#Get variables from command line
while getopts ":i:h:" optname
do
  case "$optname" in
    "h")
      echo "Prediction of plasmids sequences using cBar"
      echo ""
      echo "Options :"
      echo " -i : input contigs (Required)"
      echo " -h : THIS HELP"
      echo ""
      exit 1
      ;;
    "i")
      input_file=$OPTARG
      ;;
    "?")
      echo "Unknown option $OPTARG"
      exit 1
      ;;
    ":")
      echo "No argument value for option $OPTARG"
      exit 1
      ;;
    *)
      # Should not occur
      echo "Unknown error while processing options"
      ;;
  esac
done

# check if proper variables were set in command line
if [ -z $input_file ] ; then
  echo "missing arguments"
  echo "please ensure that all options are specified"
  echo ""
  echo "Options :"
  echo " -i : input contigs (Required)"
  echo " -h : THIS HELP"
  echo ""
  exit 1
fi

scaffolds_file=$input_file
scaffolds_file_1kb=$scaffolds_file"_1kb.fasta"
#filter inut sequences to have minimal length of 1kb
$filter_sequences_by_length_script -input $scaffolds_file -output $scaffolds_file_1kb -thresh 1000
cbar_output_file=$scaffolds_file_1kb".cbar.tsv"
cbar_plasmids=$scaffolds_file_1kb".cbar.tsv_plasmids.fasta"
cbar_chromosomes=$scaffolds_file_1kb".cbar.tsv_chromosomes.fasta"
#perform cBar predition
$cbar_path $scaffolds_file_1kb $cbar_output_file
#filter input fasta file based on predictions
$output_fasta_from_cbar_path -fasta $scaffolds_file_1kb -cBar $cbar_output_file -chromosome $cbar_chromosomes -plasmid $cbar_plasmids
