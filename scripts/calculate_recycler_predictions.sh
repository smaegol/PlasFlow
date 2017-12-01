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
samtools_path="/usr/local/bin/samtools"
recycler_path="/home/smaegol/storage/analyses/plasmids/recycler/Recycler/recycle.py"
make_fasta_from_fastg_path="/home/smaegol/storage/analyses/plasmids/recycler/Recycler/make_fasta_from_fastg.py"
bwa_path="/usr/local/bioinformatics/bwa/bwa"
spades_path="/usr/local/bioinformatics/SPAdes/SPAdes-3.9.1/bin/spades.py"

#define spades options
spades_threads=16


#Get variables from command line
while getopts ":1:2:h:" optname
do
  case "$optname" in
    "h")
      echo "Prediction of plasmids sequences using Recycler"
      echo ""
      echo "Options :"
      echo " -1 : input R1 file (required)"
      echo " -2 : input R2 file (required)"
      echo " -h : THIS HELP"
      echo ""
      exit 1
      ;;
    "1")
      R1=$OPTARG
      ;;
    "2")
      R2=$OPTARG
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
if [ -z $R1 ] || [ -z $R2 ] ; then
  echo "missing arguments"
  echo "please ensure that all options are specified"
  echo ""
  echo "Options :"
  echo " -1 : input R1 reads (Required)"
  echo " -2 : input R2 reads (Required)"
  echo " -h : THIS HELP"
  echo ""
  exit 1
fi



FILENAME_PREFIX=`expr match "$R1" '\(.*\)_R1.*'`
PREFIX_BASENAME=`expr match "$FILENAME_PREFIX" '.*\/\(.*\)'`
SPADES_FOLDER=$PREFIX_BASENAME"_spades" #name of the SPAdes assembly folder

if [ -d $SPADES_FOLDER ]; then
  #if spades folder exists it means that assembly was already done
  echo "spades folder already exists. Skipping assembly..."
else
  #run the assembly using SPAdes on provided sequences. Assembly will be stored in the $SPADES_FOLDER
  $spades_path -1 $R1 -2 $R2 -t $spades_threads -o $SPADES_FOLDER
fi

# the following commands are based on the manual available at https://github.com/Shamir-Lab/Recycler

cd $SPADES_FOLDER
#create fasta from the fastg(assembly graph) using script provided with Recycler
if [ -e "assembly_graph.nodes.fasta" ];
then
  echo "Fasta file from the fastg already prepared..."
else
  python2.7 $make_fasta_from_fastg_path -g assembly_graph.fastg
fi
#index obtained fasta sequences for BWA
if [ -e "assembly_graph.nodes.fasta.sa" ];
then
  echo "BWA index already built..."
else
  $bwa_path index assembly_graph.nodes.fasta
fi

cd ..
#use BWA to map sequences used for assembly to the assembly graph
if [ -e "$SPADES_FOLDER/paired_reads.sam" ];
then
  echo "Skipping samtools part..."
else
  $bwa_path mem -t 12 $SPADES_FOLDER/assembly_graph.nodes.fasta $R1 $R2  > $SPADES_FOLDER/paired_reads.sam
  #use Samtools to filter paired reads
  $samtools_path view -buS $SPADES_FOLDER/paired_reads.sam > $SPADES_FOLDER/reads_pe.bam
  $samtools_path view -bF 0x0800 $SPADES_FOLDER/reads_pe.bam > $SPADES_FOLDER/reads_pe_primary.bam
  $samtools_path sort -o $SPADES_FOLDER/reads_pe_primary.sort.bam $SPADES_FOLDER/reads_pe_primary.bam
  $samtools_path index $SPADES_FOLDER/reads_pe_primary.sort.bam
fi

cd $SPADES_FOLDER

#Get maximum K used by the assembler
if [ ! -z `find . -name "K127"` ]; then
  echo "$SPADES_FOLDER found 127"
  K="127"
elif [ ! -z `find . -name "K99"` ]; then
  echo "$SPADES_FOLDER found 99"
  K="99"
elif [ ! -z `find . -name "K77"` ]; then
  echo "$SPADES_FOLDER found 77"
  K="77"
elif [ ! -z `find . -name "K55"` ]; then
  echo "$SPADES_FOLDER found 55"
  K="55"
elif [ ! -z `find . -name "K33"` ]; then
  echo "$SPADES_FOLDER found 33"
  K="33"
fi

#Run Recycler
python2.7 $recycler_path -g assembly_graph.fastg -k $K -b read_pe_primary.sort.bam

cd ..
