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
my $o_insert_length   =   99;  # default insert length
my $o_insert_min      = 0.75;  # minimum insert length (-25%)
my $o_insert_max      = 1.25;  # maximum insert length (+25%)
my $o_smrtbell_length =   45;  # default SMRTbell length
my $o_smrtbell_min    = 0.75;  # minimum SMRTbell length (-25%)
my $o_smrtbell_max    = 1.25;  # maximum SMRTbell length (+25%)

GetOptions(
    "insert_length=f"   => \$o_insert_length,
    "insert_min=f"      => \$o_insert_min,
    "insert_max=f"      => \$o_insert_max,
    "smrtbell_length=f" => \$o_smrtbell_length,
    "smrtbell_min=f"    => \$o_smrtbell_min,
    "smrtbell_max=f"    => \$o_smrtbell_max,
    );

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 [options] input.bam\n";
    print "\n";
    print "options:\n";
    print " --insert_length\t($o_insert_length)\n";
    print " --insert_min\t\t($o_insert_min)\n";
    print " --insert_max\t\t($o_insert_max)\n";
    print " --smrtbell_length\t($o_smrtbell_length)\n";
    print " --smrtbell_min\t\t($o_smrtbell_min)\n";
    print " --smrtbell_max\t\t($o_smrtbell_max)\n";
    print "\n";
    exit;
}

### command-line arguments
my $bamfile = shift @ARGV;

my $prev = -1;
my @reads = ();

## read PacBio CCS reads
open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    ### skip headers
    if( substr($line, 0, 1) eq "@" )
    {
	next;
    }
    
    ### parse fields
    my @tokens = split(/\t/, $line);
    my $qname = shift @tokens;

    ### parse pacbio qname
    my ($movie,$zmw,$region) = split(/\//, $qname);
    my ($b,$e) = split(/_/, $region);
    
    ### update $prev
    if( $prev == -1 )
    {
	$prev = $zmw;
    }

    ### accumulate read names (per ZMW)
    if( $zmw ne $prev || eof(BAM) )
    {
	### find longest sequence of inserts/adapters
	pbreads(@reads);

	### reset
	$prev = $zmw;
	@reads = ([$qname, $b, $e]);
    }
    else
    {
	push @reads, [$qname, $b, $e];
    }
}

close(BAM);

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub pbreads {
    my (@reads) = @_;

    ### sort reads by qStart
    @reads = sort { $$a[1] <=> $$b[1] } @reads;
    
    ### find contigs
    my $contig_id = 0;
    my %contig_size = ();
    
    ### init
    $reads[0][3] = $contig_id;
    $contig_size{$contig_id}++;
    
    for( my $i = 1; $i < @reads; $i++ )
    {
	my $sep = $reads[$i][1] - $reads[$i-1][2]; # adapter length
	my $len = $reads[$i][2] - $reads[$i][1];   # insert length
	
	### we expect succession of inserts and adapters of specific length.
	### we allow 25% deviation from expected length to account for Pacbio variability.
	if( $sep < $o_smrtbell_length * $o_smrtbell_min ||
	    $sep > $o_smrtbell_length * $o_smrtbell_max ||
	    $len < $o_insert_length   * $o_insert_min ||
	    $len > $o_insert_length   * $o_insert_max )
	{
	    $contig_id++;
	}
	
	### assign contig_id
	$reads[$i][3] = $contig_id;

	### count inserts within contig
	$contig_size{$contig_id}++;
    }

    ### find largest contig
    my @sorted = sort { $contig_size{$b} <=> $contig_size{$a} } keys %contig_size;

    ### print subreads comprising the longest contig
    my $count = 0;
    
    for( my $i = 0; $i < @reads; $i++ )
    {
	if( $reads[$i][3] == $sorted[0] )
	{
	    printf( "%s,%i\n", $reads[$i][0], $count % 2 );
	    $count++;
	}
    }
}
