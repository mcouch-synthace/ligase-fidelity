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
my $o_bcout   = "barcodes.csv";
my $o_ohout   = "overhangs.csv";
my $o_match   = 5;
my $o_np      = 5;
my $o_verbose = 0;
my $o_etype   = "";

GetOptions(
    "np=f"     => \$o_np,
    "match=f"  => \$o_match,
    "etype=s"  => \$o_etype,
    "bcout=s"  => \$o_bcout,
    "ohout=s"  => \$o_ohout,
    "verbose+" => \$o_verbose,
    );

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 sampledir1 sampledir2 sampledir3 ...\n";
    print "\n";
    print "options:\n";
    print "  --etype\tsubstrate type (default '$o_etype')\n";
    print "  --np\t\tminimum number of passes (default '$o_np')\n";
    print "  --match\tminimum number of matchesin the barcode region (default '$o_match')\n";
    print "  --bcout\toutput file for barcodes (default '$o_bcout')\n";
    print "  --ohout\toutput file for overhangs (default '$o_ohout')\n";
    exit;
}

### command-line arguments
my @rundirs = @ARGV;

my %barcodes  = ();
my %overhangs = ();

for my $rundir ( @rundirs )
{
    my $frag_fwd = ();
    my $zmws_fwd = ();

    my $frag_rev = ();
    my $zmws_rev = ();

    if( $o_etype ne "blunt-single" )
    {
	$frag_fwd = read_zmw_data("$rundir/fragments.fwd.csv");
	$zmws_fwd = read_zmw_data("$rundir/zmws.fwd.csv");
	
	$frag_rev = read_zmw_data("$rundir/fragments.rev.csv");
	$zmws_rev = read_zmw_data("$rundir/zmws.rev.csv");
    }
    else
    {
	$frag_fwd = read_zmw_data("$rundir/fragments.csv");
	$zmws_fwd = read_zmw_data("$rundir/zmws.csv");
    }

    for my $movie ( keys %$frag_fwd )
    {
	for my $zmw ( keys %{$$frag_fwd{$movie}} )
	{
	    next if( ! exists $$frag_rev{$movie}{$zmw} && $o_etype ne "blunt-single" );

	    if( $o_etype eq "b3" )
	    {
		b3( $$frag_fwd{$movie}{$zmw},
		    $$frag_rev{$movie}{$zmw},
		    $$zmws_fwd{$movie}{$zmw},
		    $$zmws_rev{$movie}{$zmw} );
	    }
	    elsif( $o_etype eq "b4" )
	    {
		b4( $$frag_fwd{$movie}{$zmw},
		    $$frag_rev{$movie}{$zmw},
		    $$zmws_fwd{$movie}{$zmw},
		    $$zmws_rev{$movie}{$zmw} );
	    }
	    elsif( $o_etype eq "blunt" )
	    {
		blunt( $$frag_fwd{$movie}{$zmw},
		       $$frag_rev{$movie}{$zmw},
		       $$zmws_fwd{$movie}{$zmw},
		       $$zmws_rev{$movie}{$zmw} );
	    }
	    elsif( $o_etype eq "blunt-single" )
	    {
		blunt_single( $$frag_fwd{$movie}{$zmw},
			      $$zmws_fwd{$movie}{$zmw} );
	    }
	    else
	    {
		print "[WARNING] specify correct experiment type with '--etype'\n";
		exit;
	    }
	}
    }
}

write_file($o_bcout,\%barcodes,["Barcode","Count"]);

if( $o_etype eq "blunt-single" )
{
    write_file($o_ohout,\%overhangs,["Blunt","Count"]);
}
else
{
    write_file($o_ohout,\%overhangs,["O1,O2","Count"]);
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub write_file {
    my ($file,$data,$headers) = @_;
    
    open(OUT,">",$file) || die "Can't write '$file'";
    
    print OUT join(",",@$headers), "\n";
    
    for my $k ( sort { $$data{$b} <=> $$data{$a} } keys %$data )
    {
	print OUT join(",",$k,$$data{$k}), "\n";
    }
    
    close(OUT);
}

sub b3 {
    my ($fwd,$rev,$zmws_fwd,$zmws_rev) = @_;

    ### check fwd-BC1 and fwd-BC2
    return 1 if( length($$fwd{"R2"}) != 6 );
    return 1 if( length($$fwd{"R8"}) != 6 );
    
    ### check fwd-OH (strict)
    return 1 if( $$fwd{"R4"} ne "TCG" || length($$fwd{"R5"}) != 3 || $$fwd{"R6"} ne "CGA" || index($$fwd{"R5"},"-") != -1 );
    
    ### check rev-BC1 and rev-BC2
    return 1 if( length($$rev{"R2"}) != 6 );
    return 1 if( length($$rev{"R8"}) != 6 );

    ### check rev-OH (strict)
    return 1 if( $$rev{"R4"} ne "TCG" || length($$rev{"R5"}) != 3 || $$rev{"R6"} ne "CGA" || index($$rev{"R5"},"-") != -1 );
    
    ### self-consistency BC check
    my $m1 = count_matches($$fwd{"R8"},scalar reverse(complement($$rev{"R2"})));
    my $m2 = count_matches($$rev{"R8"},scalar reverse(complement($$fwd{"R2"})));
    return 1 if( $m1 < $o_match || $m2 < $o_match );

    ### require a certain numer of passes for CCS read
    return 1 if( $$zmws_fwd{"NP"} < $o_np || $$zmws_rev{"NP"} < $o_np );

    ### PacBio CCS read are reverse complements. Convert it to template-based sequences
    my $o1 = scalar reverse(complement($$fwd{"R5"}));
    my $o2 = scalar reverse(complement($$rev{"R5"}));

    ### keep overhangs sorted
    ($o1,$o2) = sort($o1,$o2);

    if( $o_verbose )
    {
	print $$fwd{"R2"}, " ", $$fwd{"R5"}, " ", $$fwd{"R8"}, "---", $$rev{"R2"}, " ", $$rev{"R5"}, " ", $$rev{"R8"}, " ", $m1, " ", $m2, "\n";
    }

    $overhangs{"$o1,$o2"}++;
    
    $barcodes{$$fwd{"R8"}}++;
    $barcodes{$$rev{"R8"}}++;
}

sub b4 {
    my ($fwd,$rev,$zmws_fwd,$zmws_rev) = @_;

    ### check fwd-BC1 and fwd-BC2
    return 1 if( length($$fwd{"R2"}) != 6 );
    return 1 if( length($$fwd{"R8"}) != 6 );
    
    ### check fwd-OH (strict)
    return 1 if( $$fwd{"R4"} ne "TCC" || length($$fwd{"R5"}) != 4 || $$fwd{"R6"} ne "GGA" || index($$fwd{"R5"},"-") != -1 );
    
    ### check rev-BC1 and rev-BC2
    return 1 if( length($$rev{"R2"}) != 6 );
    return 1 if( length($$rev{"R8"}) != 6 );

    ### check rev-OH (strict)
    return 1 if( $$rev{"R4"} ne "TCC" || length($$rev{"R5"}) != 4 || $$rev{"R6"} ne "GGA" || index($$rev{"R5"},"-") != -1 );
    
    ### self-consistency BC check
    my $m1 = count_matches($$fwd{"R8"},scalar reverse(complement($$rev{"R2"})));
    my $m2 = count_matches($$rev{"R8"},scalar reverse(complement($$fwd{"R2"})));
    return 1 if( $m1 < $o_match || $m2 < $o_match );

    ### require a certain numer of passes for CCS read
    return 1 if( $$zmws_fwd{"NP"} < $o_np || $$zmws_rev{"NP"} < $o_np );

    ### PacBio CCS read are reverse complements. Convert it to template-based sequences
    my $o1 = scalar reverse(complement($$fwd{"R5"}));
    my $o2 = scalar reverse(complement($$rev{"R5"}));

    ### keep overhangs sorted
    ($o1,$o2) = sort($o1,$o2);

    if( $o_verbose )
    {
	print $$fwd{"R2"}, " ", $$fwd{"R5"}, " ", $$fwd{"R8"}, "---", $$rev{"R2"}, " ", $$rev{"R5"}, " ", $$rev{"R8"}, " ", $m1, " ", $m2, "\n";
    }

    $overhangs{"$o1,$o2"}++;
    
    $barcodes{$$fwd{"R8"}}++;
    $barcodes{$$rev{"R8"}}++;
}

sub blunt {
    my ($fwd,$rev,$zmws_fwd,$zmws_rev) = @_;

    ### check fwd-BC1 and fwd-BC2
    return 1 if( length($$fwd{"R2"}) != 6 );
    return 1 if( length($$fwd{"R8"}) != 6 );
    
    ### check fwd-OH (strict)
    return 1 if( $$fwd{"R4"} ne "TCG" || length($$fwd{"R5"}) != 6 || $$fwd{"R6"} ne "CGA" || index($$fwd{"R5"},"-") != -1 );
    
    ### check rev-BC1 and rev-BC2
    return 1 if( length($$rev{"R2"}) != 6 );
    return 1 if( length($$rev{"R8"}) != 6 );

    ### check rev-OH (strict)
    return 1 if( $$rev{"R4"} ne "TCG" || length($$rev{"R5"}) != 6 || $$rev{"R6"} ne "CGA" || index($$rev{"R5"},"-") != -1 );
    
    ### self-consistency BC check
    my $m1 = count_matches($$fwd{"R8"},scalar reverse(complement($$rev{"R2"})));
    my $m2 = count_matches($$rev{"R8"},scalar reverse(complement($$fwd{"R2"})));
    return 1 if( $m1 < $o_match || $m2 < $o_match );

    ### require a certain numer of passes for CCS read
    return 1 if( $$zmws_fwd{"NP"} < $o_np || $$zmws_rev{"NP"} < $o_np );

    ### PacBio CCS read are reverse complements. Convert it to template-based sequences
    my $o1 = scalar reverse(complement($$fwd{"R5"}));
    my $o2 = scalar reverse(complement($$rev{"R5"}));

    ### for blunt ligation we can have an extra complementarity check
    return 1 if( $o1 ne scalar reverse(complement($o2)));

    # DON'T EVEN TRY TO SORT OVERHANGS FROM BLUNT LIGATION OR ELSE...
    # ### keep overhangs sorted
    # ($o1,$o2) = sort($o1,$o2);

    if( $o_verbose )
    {
	print $$fwd{"R2"}, " ", $$fwd{"R5"}, " ", $$fwd{"R8"}, "---", $$rev{"R2"}, " ", $$rev{"R5"}, " ", $$rev{"R8"}, " ", $m1, " ", $m2, "\n";
    }

    ### note that we keep only one overhang
    ### subsequent analysis involves "doubling", so we do not do it here.
    my $pair = sprintf( "%s,%s", substr($o1,0,3), substr($o1,3,3) );
    $overhangs{$pair}++;
    
    $barcodes{$$fwd{"R8"}}++;
    $barcodes{$$rev{"R8"}}++;
}

sub blunt_single {
    my ($fwd,$zmws_fwd) = @_;

    ### check fwd-BC1 and fwd-BC2
    return 1 if( length($$fwd{"R2"}) != 6 );
    return 1 if( length($$fwd{"R8"}) != 6 );
    
    ### check fwd-OH (strict)
    return 1 if( $$fwd{"R4"} ne "TCG" || length($$fwd{"R5"}) != 6 || $$fwd{"R6"} ne "CGA" || index($$fwd{"R5"},"-") != -1 );
    
    ### require a certain numer of passes for CCS read
    return 1 if( $$zmws_fwd{"NP"} < $o_np );

    ### PacBio CCS read are reverse complements. Convert it to template-based sequences
    my $o1 = $$fwd{"R5"};
    my $o1rc = scalar reverse(complement($o1));
    
    $overhangs{ $o1 }++;
    $overhangs{ $o1rc }++;

    $barcodes{$$fwd{"R8"}}++;
    $barcodes{scalar reverse(complement($$fwd{"R2"}))}++;

    if( $o_verbose )
    {
	print $$fwd{"R2"}, " ", $$fwd{"R5"}, " ", $$fwd{"R8"}, "\n";
	print "      ", " ", $o1rc, " ", "      ", "\n";
    }
}

sub count_matches {
    my ($bc1,$bc2) = @_;
    
    my $count = 0;

    for( my $i = 0; $i < length($bc1); $i++ )
    {
	if( substr($bc1,$i,1) eq substr($bc2,$i,1) && substr($bc1,$i,1) ne "-" )
	{
	    $count++;
	}
    }

    return $count;
}

sub complement {
    my ($seq) = @_;

    my %bases = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
	"a" => "t",
	"c" => "g",
	"g" => "c",
	"t" => "a",
	"-" => "-",
	);

    my $complement = join("",map { $bases{$_} } split(//,$seq));

    return $complement;
}

sub read_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();
    my %head = ();

    if( substr($file,-3) eq ".7z" )
    {
	open(IN,"7za e -so $file |") || die "Can't open '$file'";
    }
    elsif( substr($file,-3) eq ".gz" )
    {
	open(IN,"gzip -cd $file |") || die "Can't open '$file'";
    }
    elsif( substr($file,-4) eq ".bz2" )
    {
	open(IN,"bzip2 -cd $file |") || die "Can't open '$file'";
    }
    else
    {
	open(IN,$file) || die "Can't open '$file'";
    }
    
    while( my $line = <IN> )
    {
	chomp($line);
	
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    @head = @tokens;

	    for( my $i = 0; $i < @head; $i++ )
	    {
		if( uc($head[$i]) eq "ZMW" )
		{
		    $head[$i] = "ZMW";
		}
	    }

	    %head = map { $_ => 1 } @head;
	    
	    if( ! exists $head{"Movie"} )
	    {
		print STDERR "error :: missing 'Movie' column\n";
		exit;
	    }

	    if( ! exists $head{"ZMW"} )
	    {
		print STDERR "error :: missing 'ZMW' column\n";
		exit;
	    }
	}
	else
	{
	    my %entry = ();
	    
	    for( my $i = 0; $i < @head; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    $zmws{$entry{"Movie"}}{$entry{"ZMW"}} = \%entry;
	}
    }

    close(IN);

    return \%zmws;
}
