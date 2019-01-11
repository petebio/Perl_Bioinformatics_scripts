#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use POSIX;
use warnings;
use strict;

########################################################################################################################################################

my$infile;
my$outfile;

my$max_dist = 2;
my$help = 0;

GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "dist=i" => \$max_dist, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

my$args = 0;
if(not defined $infile){
	print STDERR "\nMissing argument -infile";
	$args = 1;
}

if(not defined $outfile){
	print STDERR "\nMissing argument -outfile";
	$args = 1;
}

if($args == 1){
	print STDERR "\n\n";
	pod2usage(-verbose => 1);
}

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

		# A motif is considered to be repeated if they are on opposite strands,
		# and the center of the motifs are <= 2bp (default) from each other (distance can be managed using the -dist flag)

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

__END__

=head1 DESCRIPTION

Remove palindromic motif repeats from a motif bed file produced by Homer's annotatePeaks.pl (if the -mbed option is used).
This avoids motif positions being reported twice in the final output.

=head1 ARGUMENTS

=over 8

=item B<-infile>
Input file. The motif bed file produced by Homer's annotatePeaks.pl function (if the -mbed option is used).

=item B<-outfile>
Output file

=item B<-dist>
Maximum distance between the centre of two motifs to be considered a repeat (default = 2bp)

=back
