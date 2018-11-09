#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use POSIX;
use warnings;
use strict;

###############################################################################################################################################

my$infile;
my$outfile;
my$l;
my$help = 0;

GetOptions("i=s" => \$infile, "o=s" => \$outfile, "l=i" => \$l, "h|?" => \$help);

pod2usage(-verbose => 99) if $help;

my$args = 1;
if(not defined $infile){
	print STDERR "\nMissing argument -i";
	$args = 0;
}

if(not defined $outfile){
	print STDERR "\nMissing argument -o";
	$args = 0;
}

if(not defined $l){
	print STDERR "\nMissing argument -l";
	$args = 0;
}

if($args == 0){
	print STDERR "\n\n";
	pod2usage(-verbose => 1);
}

###############################################################################################################################################

if($l % 2 != 0){
	print STDERR "Length of feature must be even\n";
}

my$e = $l / 2;

###############################################################################################################################################

open(IN, $infile) || die "Cannot open file $infile : $!.\n";
open(OUT, ">", $outfile);

while(my$line = <IN>){
	chomp $line;
	my@fields = split("\t", $line);

	my($start, $end) = sort {$a <=> $b} ($fields[1], $fields[2]);
	my$mid = ceil(($start + $end) / 2);

	my$new_start = $mid - $e;
	my$new_end = $mid + $e;

	if($new_start < 0){
		$new_start = 0;
	}

	if($new_end < 0){
		$new_end = 0;
	}

	$fields[1] = $new_start;
	$fields[2] = $new_end;

	print OUT join("\t", @fields),"\n";	
}

close IN;
close OUT;

###############################################################################################################################################

exit;

###############################################################################################################################################

__END__

=head1 DESCRIPTION

Resize all features in a BED file to a specified length.
This script works by first finding the mid-point of each feature, 
and extending +/-(L/2)bp from this mid-point, where L is the specified feature length. Therefore, the -l argument must be an even number.


=head1 ARGUMENTS

=over 8

=item B<-i>
Input BED file.

=item B<-o>
Output BED file.

=item B<-l>
Length of features. Must be an even number.

=item B<-h>
Print help and exit.

=back
