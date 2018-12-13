#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use POSIX;
use warnings;
use strict;

###############################################################################################################################################

my$bed_file;
my$bam_file;
my$outfile;

my$window = 100;

my$help = 0;

GetOptions("bed=s" => \$bed_file, "bam=s" => \$bam_file, "outfile=s" => \$outfile, "window=i" => \$window, "h|?" => \$help);

pod2usage(-verbose => 99) if $help;

###############################################################################################################################################

open(BED, $bed_file) || die "Cannot open file $bed_file : $!.\n";
open(TEMP, ">", "extended.bed");

while(my$line = <BED>){
	chomp $line;
	my@cols = split("\t", $line);

	my$chr = $cols[0];
	my$pos = ceil(($cols[1] + $cols[2]) / 2);

	my$start = $pos - $window;
	my$end = $pos + $window;

	print TEMP $chr,"\t",$start,"\t",$end,"\n";
}

close BED;
close TEMP;

###############################################################################################################################################

my$command = "dnase_wig_tracks.py extended.bed $bam_file forward.wig reverse.wig";
system($command);

unlink "extended.bed";

my@forward_cuts = get_average_counts("forward.wig");
my@reverse_cuts = get_average_counts("reverse.wig");

unlink "reverse.wig";
unlink "forward.wig";

open(OUT, ">", $outfile);
print OUT "pos\tforward\treverse\n";

for(my$i = 0; $i < scalar @forward_cuts; $i++){
	my$pos = $i + 1;
	print OUT $pos,"\t",$forward_cuts[$i],"\t",$reverse_cuts[$i],"\n";
}

close OUT;

###############################################################################################################################################

exit;

###############################################################################################################################################

sub get_average_counts{
	my$wig_file = shift;

	open(WIG, $wig_file) || die "Cannot open file $wig_file : $!.\n";

	my$site_count = 0;
	my@cuts;

	my$i = 0;
	while(my$line = <WIG>){
		chomp $line;
		next if $line =~ m/^track/;

		if($line =~ m/^fixed/){
			$i = 0;
			$site_count++;
		}
		else{
			if(not defined $cuts[$i]){
				push(@cuts, 0);
			}

			$cuts[$i] += abs($line);
			$i++;
		}
	}

	close WIG;

	my@average_cuts;
	foreach my$cut (@cuts){
		my$average = $cut / $site_count;
		push(@average_cuts, $average);
	}

	return @average_cuts;
}

###############################################################################################################################################

__END__

=head1 DESCRIPTION

Retrieve the average number of DNase cuts around a set of footprints.

=head1 ARGUMENTS

=over 8

=item B<-bed>
BED file of footprint coordinates.

=item B<-bam>
BAM file.

=item B<-outfile>
Ouptut file (.tsv)

=item B<-window>
Window size to plot the profile (default is +/- 100bp).

=item B<-help>
Print help and exit.

=back