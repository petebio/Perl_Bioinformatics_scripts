#!/usr/bin/env perl
use Getopt::Long;
use warnings;
use strict;

#########################################################################################################################################

my$bed_file_a;
my$bed_file_b;
my$outfile;

GetOptions("a=s" => \$bed_file_a, "b=s" => \$bed_file_b, "o=s" => \$outfile);

#########################################################################################################################################

my$tmp_dir = $outfile."_temp";
mkdir $tmp_dir;

extract_pos_columns($bed_file_a, "$tmp_dir/BedA.bed");
extract_pos_columns($bed_file_b, "$tmp_dir/BedB.bed");

my$intersect_command = "bedtools intersect -wo -a $tmp_dir/BedA.bed -b $tmp_dir/BedB.bed > $tmp_dir/Intersect.bed";
system $intersect_command;

open(INT, "$tmp_dir/Intersect.bed") || die "Cannot open file $tmp_dir/Intersect.bed : $!.\n";
open(MERGE, ">", "$tmp_dir/ToMerge.bed");

while(my$line = <INT>){
	chomp $line;
	my@fields = split("\t", $line);

	print MERGE $fields[0],"\t",$fields[1],"\t",$fields[2],"\n";
	print MERGE $fields[3],"\t",$fields[4],"\t",$fields[5],"\n";
}

close INT;
close MERGE;

my$merge_command = "bedtools sort -i $tmp_dir/ToMerge.bed | bedtools merge -i - > $outfile";
system $merge_command;

system "rm -r $tmp_dir";

#########################################################################################################################################

exit;

#########################################################################################################################################

#### SUBROUTINES ####

sub extract_pos_columns{
	my($infile, $outfile) = @_;

	open(IN, $infile) || die "Cannot open file $infile : $!.\n";
	open(OUT, ">", $outfile);
	
	while(my$line = <IN>){
		chomp $line;
		my@pos = $line =~ m/^(\S+)\s+(\S+)\s+(\S+)\s*/;
		print OUT join("\t", @pos),"\n";
	}

	close IN;
	close OUT;

	return 1;
}

#########################################################################################################################################
