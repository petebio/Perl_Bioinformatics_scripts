#!/usr/bin/env perl
use Getopt::Long;
use POSIX;
use warnings;
use strict;

########################################################################################################################################################

my$infile;
my$outfile;

my$max_dist = 2;

GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "dist=i" => \$max_dist);

########################################################################################################################################################

open(IN, $infile) || die "Cannot open file $infile : $!.\n";
open(OUT, ">", $outfile);

my@prev_line;
while(my$line = <IN>){
	chomp $line;
	next if $line !~ m/^chr/;

	my@cols = split("\t", $line);

	my$chr = $cols[0];
	my$start = $cols[1];
	my$end = $cols[2];
	my$strand = $cols[5];

	my$pos = ceil(($start + $end) / 2);

	my$pass = 1;
	if(scalar @prev_line != 0){
		my($prev_chr, $prev_pos, $prev_strand) = @prev_line;
		my$dist = abs($pos - $prev_pos);	

		if($chr eq $prev_chr && $strand ne $prev_strand && $dist <= $max_dist){
			$pass = 0;
		}
	}

	if($pass == 1){
		print OUT $chr,"\t",$start,"\t",$end,"\n";
	}

	@prev_line = ($chr, $pos, $strand);
}

close IN;
close OUT;

########################################################################################################################################################

exit;

########################################################################################################################################################
