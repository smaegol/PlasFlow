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
use warnings;

my $ids_file = $ARGV[0];
my $big_file = $ARGV[1];
my $output   = $ARGV[2];

my %my_hash = ();

my $i = 0;
open( IDS, "<$ids_file" );
while ( my $line = <IDS> ) {
    my ($id) = ( $line =~ /^(\S+)/ );
    $my_hash{$id} = 1;
}

close(IDS);

foreach my $a ( sort keys %my_hash ) {
    print "$a\n";
}

my $accession;

open( BIG, "<$big_file" );
open( OUT, ">$output" );
while ( my $line = <BIG> ) {
    $i++;
    if ( $i % 100000 == 0 ) {
        print "$i sequences processed\n";

    }

    ($accession) = ( $line =~ /.+\t.+\t(\S+)\...?\t.+\t.+\t.+\t.+$/ );
    if ( $my_hash{$accession} ) {
        print OUT $line;
        print "$accession found!!\n";
    }
}
close(OUT);
close(BIG);
