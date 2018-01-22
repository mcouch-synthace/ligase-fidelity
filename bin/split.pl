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
my $o_fwd = "subreads.fwd.sam";
my $o_rev = "subreads.rev.sam";

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 bamfile [options] clusters.bz2\n";
    print "\n";
    print "options:\n";
    print "  --fwd\tfile name for forward reads ($o_fwd)\n";
    print "  --rev\tfile name for reverse reads ($o_rev)\n";
    print "\n";
    exit;
}

my $bamfile = shift @ARGV;
my $csvfile = shift @ARGV;

### read mapping information
my %csv = ();

open(CSV,"bzip2 -cd $csvfile |") || die "Can't open '$csvfile'";

while( my $line = <CSV> )
{
    chomp($line);
    
    my ($qname,$flag) = split(/,/,$line);

    ### save read direction
    $csv{$qname} = $flag;
}

close(CSV);

### split BAM file
open(FWD,">",$o_fwd) || die "Can't write '$o_fwd'";
open(REV,">",$o_rev) || die "Can't write '$o_rev'";

open(BAM,"samtools view -h $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    if( substr($line,0,1) eq "@" )
    {
	### copy BAM headers to both files
	print FWD $line;
	print REV $line;
    }
    else
    {
	my ($qname,$flag,$rname,$pos,$mapq,@other) = split(/\t/,$line);
	
	if( exists $csv{$qname} )
	{
	    if( $csv{$qname} == 0 )
	    {
		print FWD $line;
	    }
	    else
	    {
		print REV $line;
	    }
	}
    }
}

close(BAM);

close(FWD);
close(REV);
