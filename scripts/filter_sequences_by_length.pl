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

my $input_file;
my $output_file;
my $length_threshold;

# opis uzycia programu
my $USAGE =
"Filters sequences in multifasta file by length\n\nUsage of program:\n\n $0 [options] \n\t-input\t\tmultifasta file with sequences to be filtered\n\t-output\t\toutput file for filtered sequences \n\t-thresh\t\tsequence length threshold \n";

#pobranie zmiennych z linii polecen
unless (@ARGV) {
    print $USAGE;
    exit;
}

GetOptions(
    'input=s'  => \$input_file,
    'thresh=s' => \$length_threshold,
    'output=s' => \$output_file
);

if ( $input_file eq '' ) {
    die "no input file specified\n";
}

if ( $output_file eq '' ) {
    die "no output file specified\n";
}

if ( $length_threshold eq '' ) {
    die "no length threshold provided\n";
}

my $inseq;

#properly treat gzipped files
if ( $input_file =~ /\.gz$/ ) {
    $inseq = Bio::SeqIO->new(
        -file   => "gunzip -c $input_file |",
        -format => 'fasta',
    );
}

else {
    $inseq = Bio::SeqIO->new(
        -file   => "<$input_file",
        -format => 'fasta',
    );

}

#define output fasta file
my $out_file = Bio::SeqIO->new(
    -file   => ">$output_file",
    -format => 'fasta',
);

#counter variables:
my $no_seq      = 0;
my $missed_seqs = 0;

my $seqLength;

#iterate over sequences in the input file
while ( my $seq = $inseq->next_seq ) {
    $seqLength = $seq->length;
    if ( $seqLength >= $length_threshold ) {
        $out_file->write_seq($seq);
    }
    else {
        #		print "Sequence $seqID of length $seqLength will not be included\n";
        $missed_seqs++;
    }
    $no_seq++;
}

#remaining sequences
my $remaining = $no_seq - $missed_seqs;

print
"$missed_seqs out of $no_seq sequences from $input_file were not included in the dataset $output_file.
Remaining $remaining seqeunces will be used for further processing steps
";

exit;
