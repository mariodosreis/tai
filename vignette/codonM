#! /usr/bin/perl -w
#
# --------------------------------------------------------------------
# CodonM -- Codon Matrix:
# 
# This program generates a matrix of codon frequencies per ORF from
# an input file containing a list of coding sequences in FastA format.
# The output of this program should be useful for calculating a series
# of codon usage bias indexes.
#
# (c) 2003 Mario dos Reis.
#
# This file forms part of codonR, a package for the analysis of codon
# usage in DNA coding sequences.
#
# codonR is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# Mario dos Reis
# m.reis@mail.cryst.bbk.ac.uk
# School of Crystallography
# Birkbeck College
# London WC1E 7HX
#
# --------------------------------------------------------------------
#
# version: 0.4 (Nov 2003)
#
# codonM [INFILE] [OUTFILE]

use strict;
use warnings;

my $version = 0.4;

if(@ARGV < 2) {
    print "\nCodonM version: $version\n",
    "usage: codonM <INFILE> <OUTFILE>\n\n";
    exit;
}

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(FASTA_FILE, $infile) 
    or die "Fatal: cannot open file $infile:\n$!";

open(CODON_MATRIX, ">$outfile")
    or die "Fatal: couldn't open $outfile:\n$!";

my @file = <FASTA_FILE>;

my @tcag = ('T', 'C', 'A', 'G');
my (@codons) = Codons(@tcag);

my @sequences = ReadFasta(@file);

print "Analysing $infile, please be patient ...\n";

foreach my $seq (@sequences) {

    my @count;

    for(my $k = 0; $k < 61; $k++) { $count[$k] = 0 }

    # $i = 3 because we are ignoring START codons
    for(my $i = 3; $i < length($seq); $i += 3) {

	my $codon = substr($seq, $i, 3);

	# the whole list, for optimisation reasons
	# STOP codons have been ignored
	if($codon =~ /TTT/o) { $count[0]++ }
	elsif($codon =~ /TTC/o) { $count[1]++ }
	elsif($codon =~ /TTA/o) { $count[2]++ }
	elsif($codon =~ /TTG/o) { $count[3]++ }
	elsif($codon =~ /TCT/o) { $count[4]++ }
	elsif($codon =~ /TCC/o) { $count[5]++ }
	elsif($codon =~ /TCA/o) { $count[6]++ }
	elsif($codon =~ /TCG/o) { $count[7]++ }
	elsif($codon =~ /TAT/o) { $count[8]++ }
	elsif($codon =~ /TAC/o) { $count[9]++ }
	elsif($codon =~ /TGT/o) { $count[10]++ }
	elsif($codon =~ /TGC/o) { $count[11]++ }
	elsif($codon =~ /TGG/o) { $count[12]++ }
	elsif($codon =~ /CTT/o) { $count[13]++ }
	elsif($codon =~ /CTC/o) { $count[14]++ }
	elsif($codon =~ /CTA/o) { $count[15]++ }
	elsif($codon =~ /CTG/o) { $count[16]++ }
	elsif($codon =~ /CCT/o) { $count[17]++ }
	elsif($codon =~ /CCC/o) { $count[18]++ }
	elsif($codon =~ /CCA/o) { $count[19]++ }
	elsif($codon =~ /CCG/o) { $count[20]++ }
	elsif($codon =~ /CAT/o) { $count[21]++ }
	elsif($codon =~ /CAC/o) { $count[22]++ }
	elsif($codon =~ /CAA/o) { $count[23]++ }
	elsif($codon =~ /CAG/o) { $count[24]++ }
	elsif($codon =~ /CGT/o) { $count[25]++ }
	elsif($codon =~ /CGC/o) { $count[26]++ }
	elsif($codon =~ /CGA/o) { $count[27]++ }
	elsif($codon =~ /CGG/o) { $count[28]++ }
	elsif($codon =~ /ATT/o) { $count[29]++ }
	elsif($codon =~ /ATC/o) { $count[30]++ }
	elsif($codon =~ /ATA/o) { $count[31]++ }
	elsif($codon =~ /ATG/o) { $count[32]++ }
	elsif($codon =~ /ACT/o) { $count[33]++ }
	elsif($codon =~ /ACC/o) { $count[34]++ }
	elsif($codon =~ /ACA/o) { $count[35]++ }
	elsif($codon =~ /ACG/o) { $count[36]++ }
	elsif($codon =~ /AAT/o) { $count[37]++ }
	elsif($codon =~ /AAC/o) { $count[38]++ }
	elsif($codon =~ /AAA/o) { $count[39]++ }
	elsif($codon =~ /AAG/o) { $count[40]++ }
	elsif($codon =~ /AGT/o) { $count[41]++ }
	elsif($codon =~ /AGC/o) { $count[42]++ }
	elsif($codon =~ /AGA/o) { $count[43]++ }
	elsif($codon =~ /AGG/o) { $count[44]++ }
	elsif($codon =~ /GTT/o) { $count[45]++ }
	elsif($codon =~ /GTC/o) { $count[46]++ }
	elsif($codon =~ /GTA/o) { $count[47]++ }
	elsif($codon =~ /GTG/o) { $count[48]++ }
	elsif($codon =~ /GCT/o) { $count[49]++ }
	elsif($codon =~ /GCC/o) { $count[50]++ }
	elsif($codon =~ /GCA/o) { $count[51]++ }
	elsif($codon =~ /GCG/o) { $count[52]++ }
	elsif($codon =~ /GAT/o) { $count[53]++ }
	elsif($codon =~ /GAC/o) { $count[54]++ }
	elsif($codon =~ /GAA/o) { $count[55]++ }
	elsif($codon =~ /GAG/o) { $count[56]++ }
	elsif($codon =~ /GGT/o) { $count[57]++ }
	elsif($codon =~ /GGC/o) { $count[58]++ }
	elsif($codon =~ /GGA/o) { $count[59]++ }
	elsif($codon =~ /GGG/o) { $count[60]++ }
    }
    for(my $j = 0; $j < @count; $j++) {
	print CODON_MATRIX "$count[$j]\t";
    }
    print CODON_MATRIX "\n";
}

print "\tDONE\n";

close FASTA_FILE;
close CODON_MATRIX;


# --------------------------------------------------------------------
# SUBS start here:
# --------------------------------------------------------------------

# Creates the subset of 64 codons
sub Codons {

    my(@tcag) = @_;
    my(@codons);
    
    for (my $i = 0; $i < @tcag; $i++) {

	for (my $j = 0; $j < @tcag; $j++) {

	    for (my $k = 0; $k < @tcag; $k++) {
		my $codon = $tcag[$i].$tcag[$j].$tcag[$k];
		push (@codons, $codon);
	    }
	}
    }

    return @codons;
}

# Reads and pre-proccess array of FastA seqs
sub ReadFasta {

    my (@file) = @_;;
    
    my @seqs;
    
    my $seq = "";
    
    # load sequences
    foreach my $line (@file) {
	
	if ($line =~ />/) {
	    push(@seqs, $seq);
	    $seq = "";
	    next;
	
	} else {

	    chomp($line);
	    $line =~ tr/atcg/ATCG/; # to UPPERCASE
	    $seq = $seq.$line;
	    
	}
    }

    # add last sequence
    push(@seqs, $seq);
    
    # delete first element of the array (empty sequence)
    shift @seqs;

    return @seqs;
}
