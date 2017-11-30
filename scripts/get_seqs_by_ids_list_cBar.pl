#!/usr/bin/perl

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

use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

#define variables
my $input_file;
my $output_chromosomes;
my $output_plasmids;
my $output_file;
my $names;
my $i;

my $random_number;

#read options from command line
GetOptions(
    'fasta=s'      => \$input_file,
    'cBar=s'       => \$names,
    'chromosome=s' => \$output_chromosomes,
    'plasmid=s'    => \$output_plasmids
);

#create hash tables for storing ids of classified sequences
my %names       = ();
my %plasmids    = ();
my %chromosomes = ();

#open cBar predictions file and get predictions for each sequence
open( LIST, "<$names" );
while (<LIST>) {
    chomp $_;
    my ( $name, $length, $prediction ) = (/(\S+)\s+(\d+)\s+(\S+)/);
    if ( $prediction eq 'Plasmid' ) {
        $plasmids{$name} = 1;
    }
    else {
        $chromosomes{$name} = 1;
    }
}

#open fasta file used for classification by cBar
my $inseq = Bio::SeqIO->new(
    -file   => "<$input_file",
    -format => 'fasta',
);

#open output file for plasmid sequences
my $out_plasmids = Bio::SeqIO->new(
    -file   => ">$output_plasmids",
    -format => 'fasta',
);

#open output file for chromosome sequences
my $out_chromosomes = Bio::SeqIO->new(
    -file   => ">$output_chromosomes",
    -format => 'fasta',
);

my $seq_name;

#iterate through sequences in input fasta file
while ( my $seq = $inseq->next_seq ) {
    $seq_name = $seq->id;

    #store sequences in the proper output files:
    if ( $plasmids{$seq_name} ) {
        $out_plasmids->write_seq($seq);
    }
    else {
        $out_chromosomes->write_seq($seq);
    }
}

exit;
