#!/usr/bin/perl -w

################################################################################
# Profiling of T4 DNA ligase fidelity and bias via single-molecule sequencing
# Copyright (C) 2018 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;
use Getopt::Long;

### command-line options
my $verbose = 0;
my @regions = ();

GetOptions( 
    "region=s" => \@regions,
    "verbose+" => \$verbose,
    );

@regions = split(/,/, join(",", @regions));

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 input.bam reference.fasta\n";
    print "\n";
    print "options:\n";
    print "  --verbose\n";
    print "  --region\n";
    print "\n";
    exit;
}

### command-line arguments
my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;

my ($rseq, $name) = read_refseq($reffile);

### data structures for parsing CIGAR string
my %digit = (
    "0" => 1,
    "1" => 1,
    "2" => 1,
    "3" => 1,
    "4" => 1,
    "5" => 1,
    "6" => 1,
    "7" => 1,
    "8" => 1,
    "9" => 1,
    );

my %operator = (
    "M" => 1,
    "I" => 1,
    "D" => 1,
    "=" => 1,
    "X" => 1,
    );

### number of regions to parse
my $numreg = scalar(@regions);

### save for future reference
print STDERR join("\n", map { "R$_,$regions[$_-1]" } 1..$numreg), "\n";

### print table headers
print join(",", "Movie", "ZMW", map { "R$_" } 1..$numreg), "\n";

open(BAM, "samtools view $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    ### skip BAM headers
    if( substr($line, 0, 1) eq "@" )
    {
	next;
    }

    my $s1 = ""; # aligned reference
    my $s2 = ""; # aligned read
    my $bq = ""; # base quality
    
    ### parse fields
    my ( $qname,
	 $flag,
	 $rname,
	 $pos,
	 $mapq,
	 $cigar,
	 $rnext,
	 $pnext,
	 $tlen,
	 $seq,
	 $qual,
	 @tag ) = split(/\t/,$line);
    
    ### strict filtering
    if( $flag != 0 )
    {
	next;
    }
    
    # ----- process CIGAR ------------------------------------------------------
    
    my $p1  = $pos - 1; # start of aligned reference
    my $p2  = 0;        # start of aligned read
    my $num = "";       # length of the block to operate
    
    my $read = uc($seq);
    my @cigar = split(//,$cigar);
    
    while( @cigar )
    {
	my $c = shift @cigar;
	
	if( exists $digit{$c} )
	{
	    $num .= $c;
	    next;
	}
	
	if( $c eq "M" || $c eq "=" || $c eq "X" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);
	    
	    $p1 += $num;
	    $p2 += $num;
	}
	elsif( $c eq "D" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= ("-" x $num);
	    $bq .= (" " x $num);
	    
	    $p1 += $num;
	}
	elsif( $c eq "I" )
	{
	    $s1 .= ("-" x $num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);
	    
	    $p2 += $num;
	}
	elsif( $c eq "S" )
	{
	    if( $p2 == 0 )
	    {
		# ================== #
		# reconstruct 5'-end #
		# ================== #

		my ($max,$min) = sort { $b <=> $a } ($p1,$num);

		my $diff = $max - $min;

	    	if( $p1 < $num )
	    	{
		    ### reference
	    	    $s1 .= "-" x $diff;
		    $s1 .= substr($rseq,0,$p1);

		    ### read
	    	    $s2 .= substr($read,$p2,$max);

		    ### qual
		    $bq .= substr($qual,$p2,$max);
	    	}
	    	else
	    	{
		    ### reference
	    	    $s1 .= substr($rseq,0,$max);

		    ### read
		    $s2 .= "-" x $diff;
		    $s2 .= substr($read,$p2,$num);

		    ### qual
		    $bq.= " " x $diff;
		    $bq .= substr($qual,$p2,$num);
	    	}
	    }
	    else
	    {
		# ================== #
		# reconstruct 3'-end #
		# ================== #

		my ($max,$min) = sort { $b <=> $a } (length($rseq)-$p1,$num);

		my $diff = $max - $min;

		if( length($rseq)-$p1 < $num )
		{
		    ### reference
		    $s1 .= substr($rseq,$p1,$num);
		    $s1 .= "-" x $diff;

		    ### read
		    $s2 .= substr($read,$p2,$max);

		    ### qual
		    $bq .= substr($qual,$p2,$max);
		}
		else
		{
		    ### reference
		    $s1 .= substr($rseq,$p1,$max);

		    ### read
		    $s2 .= substr($read,$p2,$num);
		    $s2 .= "-" x $diff;

		    ### qual
		    $bq .= substr($qual,$p2,$num);
		    $bq .= "-" x $diff;
		}
	    }
	    
	    $p2 += $num;
	}
	else
	{
	    print STDERR "[ERROR] unknown operator in CIGAR string '$c'\n\n";
	    exit;
	}
	
	$num = "";
    }

    my @data = ();

    ### extract sepcified regions
    for my $region ( @regions )
    {
	### start, end, trim flag
	my ($p1, $p2, $trim) = split(/[-\/]/, $region);

	$trim = 0 if( ! defined $trim );

	### find position in the alignment corresponding to a particular
	### sequence position
	my $b = seq2aln([split(//, $s1)], $p1);
	my $e = seq2aln([split(//, $s1)], $p2);
	
	my $str = "NA";
	
	### extract required substring
	if( $b != -1 && $e != -1 )
	{
	    $str = substr($s2,$b,$e-$b+1);

	    if( $trim > 0 )
	    {
		### trim if necessary
		$str = substr($str,$trim,length($str)-2*$trim);
	    }
	}

	push @data, $str;
    }

    my ($movie,$zmw,@pth) = split(/\//,$qname);

    ### print combined results
    print join(",",$movie,$zmw,@data), "\n";
}

close(BAM);

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub read_refseq {
    my ($file) = @_;
    
    my $rseq = "";
    my $name = "";

    open(REF,$file) || die "Can't open '$file'";
    
    while( my $line = <REF> )
    {
	chomp( $line );
	
	if( index($line,">") == 0 )
	{
	    $name = substr($line,1);
	    next;
	}	
	
	$rseq .= $line;
    }
    
    close(REF);
    
    ### upcase reference sequence once and for all
    $rseq = uc($rseq);

    return ($rseq,$name);
}

sub seq2aln {
    my ($sequence,$seqpos) = @_;

    my $pos = 0;
    
    for( my $i = 0; $i < @$sequence; $i++ )
    {
	if( $$sequence[$i] ne "-" )
	{
	    $pos++;
	    
	    if( $pos == $seqpos )
	    {
		return $i;
	    }
	}
    }
    
    return -1;
}
