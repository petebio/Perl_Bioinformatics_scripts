#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

#########################################################################################################################################################################

my$dir_path;
my$out_prefix;

my$help = 0;

GetOptions("dir=s" => \$dir_path, "out=s" => \$out_prefix, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

my$args = 0;
if(not defined $dir_path){
	print "\nMissing argument -dir";
	$args = 1;
}


if(not defined $out_prefix){
	print "\nMissing argument -out";
	$args = 1;
}

if($args == 1){
	print "\n\n";
	pod2usage(-verbose => 1);
}

#########################################################################################################################################################################

opendir(DIR, $dir_path) || die "Cannot read directory $dir_path : $!.\n";
my@tsv_file_list = grep {/.+\.tsv/} readdir DIR;
closedir DIR;

my@sample_id_list;
my%gene_exprs_tracking;
foreach my$file (@tsv_file_list){
	my($sample_id) = $file =~ m/^(.+)\.tsv/;
	push(@sample_id_list, $sample_id);

	open(IN, "$dir_path/$file") || die "Cannot open file $dir_path/$file : $!.\n";

	while(my$line = <IN>){
		next if $. == 1;
		chomp $line;

		my@cols = split("\t", $line);
		my$gene_id = $cols[0];
		my$fpkm = $cols[7];
		my$tpm = $cols[8];

		if(not exists $gene_exprs_tracking{$gene_id}{$sample_id}){
			$gene_exprs_tracking{$gene_id}{$sample_id} = [($fpkm, $tpm)];
		}
		else{
			$gene_exprs_tracking{$gene_id}{$sample_id}[1] += $fpkm;
			$gene_exprs_tracking{$gene_id}{$sample_id}[2] += $tpm;
		}
	}

	close IN;
}

#########################################################################################################################################################################

@sample_id_list = sort {$a cmp $b} @sample_id_list;

open(FPKM, ">", "$out_prefix\_FPKM.tsv");
open(TPM, ">", "$out_prefix\_TPM.tsv");

print FPKM "gene_id\t",join("\t", @sample_id_list),"\n";
print TPM "gene_id\t",join("\t", @sample_id_list),"\n";

foreach my$gene_id (sort {$a cmp $b} keys %gene_exprs_tracking){
	print FPKM $gene_id;
	print TPM $gene_id;

	foreach my$sample_id (@sample_id_list){
		my$fpkm = '0.000000';
		my$tpm = '0.000000';

		if(exists $gene_exprs_tracking{$gene_id}{$sample_id}){
			($fpkm, $tpm) = @{$gene_exprs_tracking{$gene_id}{$sample_id}};
		}

		print FPKM "\t",$fpkm;
		print TPM "\t",$tpm;
	}

	print FPKM "\n";
	print TPM "\n";
}

close FPKM;
close TPM;

#########################################################################################################################################################################

exit;

#########################################################################################################################################################################

__END__

=head1 DESCRIPTION

Read all .tsv files from a directory containing stringie abundace estimates, and merge all gene expression values into a single file. 
Produces two output files, one for FPKM values and the other for TPM.

=head1 ARGUMENTS

=over 8

=item B<-dir>

Input directory.

=item B<-out>

Prefix for the output files.

=item B<-help>

Print help and exit.

=back









