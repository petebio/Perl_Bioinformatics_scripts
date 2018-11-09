#!/usr/bin/env perl
use Getopt::Long;
use warnings;
use strict;

my$dist = 2000;

GetOptions("d=i" => \$dist);

################################################################################################################

open(PROX, ">", "Proximal.bed");
open(DIST, ">", "Distal.bed");

while(my$line = <STDIN>){
	next if $. == 1;
	
	chomp $line;
	my@fields = split("\t", $line);
	
	my$peak_id = $fields[0];
	my$chr = $fields[1];
	my$start = $fields[2];
	my$end = $fields[3];
	my$strand = $fields[4];
	
	my$tss_dist = $fields[9];
	next if not defined $tss_dist;
	next if $tss_dist eq "NA";
	$tss_dist = abs($tss_dist);
	
	my$gene_id = $fields[15];
	
	if($tss_dist <= $dist){
		print PROX $chr,"\t",$start,"\t",$end,"\n";
	}
	else{
		print DIST $chr,"\t",$start,"\t",$end,"\n";
	}
}

close PROX;
close DIST;

################################################################################################################

exit;

################################################################################################################
