#!/usr/bin/env perl
use Getopt::Long;
use warnings;
use strict;

#############################################################################################################################

my$sample_id_a;
my$sample_id_b;
my$outfile;

GetOptions("a=s" => \$sample_id_a, "b=s" => \$sample_id_b, "o=s" => \$outfile);

if(not defined $sample_id_a){
	$sample_id_a = "Sample_A";
}

if(not defined $sample_id_b){
	$sample_id_b = "Sample_B";
}

#############################################################################################################################

open(OUT, ">", $outfile);
print OUT "chr\tstart\tend\tgene_ID\t$sample_id_a\t$sample_id_b\n";

while(my$line = <STDIN>){
	chomp $line;
	next if $. == 1;

	my@fields = split("\t", $line);

	my$chr = $fields[1];
	my$start = $fields[2];
	my$end = $fields[3];
	my$gene_id = $fields[15];
	my$count_sample_a = $fields[19];
	my$count_sample_b = $fields[20];

	if($count_sample_a eq ""){
		$count_sample_a = 0;
	}

	if(not defined $count_sample_b){
		$count_sample_b = 0;
	}

	print OUT $chr,"\t",$start,"\t",$end,"\t",$gene_id,"\t",$count_sample_a,"\t",$count_sample_b,"\n";
}

close OUT;

#############################################################################################################################

exit;

#############################################################################################################################
