#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

##############################################################################################################################################

my$rna_file;
my$dhs_file;
my$outfile;
my$help = 0;

GetOptions("rna=s" => \$rna_file,
	"dhs=s" => \$dhs_file,
	"os=s" => \$outfile,
	"h|?" => \$help);

pod2usage(-verbose => 99) if $help;

my$args = 1;
if(not defined $rna_file){
	print STDERR "\nMissing argument -rna";
	$args = 0;
}

if(not defined $dhs_file){
	print STDERR "\nMissing argument -dhs";
	$args = 0;
}

if(not defined $outfile){
	print STDERR "\nMissing argument -o";
	$args = 0;
}

if($args == 0){
	print STDERR "\n\n";
	pod2usage(-verbose => 1);
}

##############################################################################################################################################

open(RNA, $rna_file) || die "Cannot open file $rna_file : $!.\n";

my%rna;
while(my$line = <RNA>){
	next if $. == 1;
	chomp $line;

	my@fields = split("\t", $line);
	my$id = $fields[0];
	my$fc = $fields[1];

	$rna{$id} = $fc;
}

close RNA;

##############################################################################################################################################

open(DHS, $dhs_file) || die "Cannot open file $dhs_file : $!.\n";
open(CDT, ">", $outfile);

print CDT "GID\tCOORD\tNAME\tGWEIGHT\tFC\n";
print CDT "AID\t1\n";
print CDT "EWEIGHT\t1\n";

my$n = 1;
while(my$line = <DHS>){
	next if $. == 1;
	chomp $line;
	
	my@fields = split("\t", $line);
	my$coord = $fields[0].":".$fields[1]."-".$fields[2];

	my$gene_ID = $fields[3];

	my$fc = 0;
	if(exists $rna{$gene_ID}){
		$fc = $rna{$gene_ID};
	}

	print CDT $n,"\t",$coord,"\t",$gene_ID,"\t1\t",$fc,"\n";
	$n++;
}

close DHS;
close CDT;

##############################################################################################################################################

exit;

##############################################################################################################################################

__END__

=head1 DESCRIPTION

Read a table of fold-change values from RNA-Seq, and an ordered, annotated BED file of DHSs and their corresponding closest gene. 
Produce a CDT file which contains the log fold-change values for the closest gene for each DHS.
If no fold-change value is available for that gene, set the fold-change value to 0.

=head1 ARGUMENTS

=over 8

=item B<-rna>
RNA-Seq file. This file should have at least two columns, where the first column is the gene ID and the second is the fold-change value.
Assumes headers are present. 

=item B<-dhs>
DHS file. Essentially a BED file with a header. First 3 columns should be position (chromosome, start, end), and the third column is
the closest gene ID.

=item B<-o>
Outfile. A CDT file in the same order as the DHS file. 

=item B<-h>
Print help and exit.

=back
