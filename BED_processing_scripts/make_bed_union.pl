#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use POSIX;
use warnings;
use strict;

###########################################################################################################################################################################

my$bed_file_a;
my$bed_file_b;
my$outfile;

my$summit_dist = 200;
my$final_length = 400;

my$help = 0;

GetOptions("a=s" => \$bed_file_a, 
	"b=s" => \$bed_file_b, 
	"outfile=s" => \$outfile, 
	"dist=s" => \$summit_dist, 
	"length=s" => \$final_length,
	"help|?" => \$help);

pod2usage(-verbose => 99) if $help;

my$args = 0;
if(not defined $bed_file_a){
	print STDERR "\nMissing argument -a";
	$args = 1;
}

if(not defined $bed_file_b){
	print STDERR "\nMissing argument -b";
	$args = 1;
}

if(not defined $outfile){
	print STDERR "\nMissing argument -outfile";
	$args = 1;
}

#if($bed_file_a eq $bed_file_b){
#	print STDERR "\nError : Input files must be different";
#	$args = 1;
#}

if($summit_dist =~ m/\./){
	print STDERR "\nError : -dist must be an integer";
	$args = 1;
}

if($final_length =~ m/\./){
	print STDERR "\nError: -length must be in integer";
	$args = 1
}

if($args == 1){
	print STDERR "\n\n";
	pod2usage(-verbose => 1);
}

###########################################################################################################################################################################

my$extend = ceil($summit_dist / 2);

my$n_lines_a = prep_bed_file($bed_file_a, $extend);
my$n_lines_b = prep_bed_file($bed_file_b, $extend);

print "\n---------------------------\n\n";
print "Read $n_lines_a lines from $bed_file_a\n";
print "Read $n_lines_b lines from $bed_file_b\n";
print "\n---------------------------\n\n";

###########################################################################################################################################################################

system "bedtools intersect -wo -a $bed_file_a.extend -b $bed_file_b.extend > Intersect.bed";

open(INT, "Intersect.bed");

my@peak_list;
my%peak_tracking;

while(my$line = <INT>){
	chomp $line;
	my@fields = split("\t", $line);

	my$chr = $fields[0];
	my$pos_a = ceil(($fields[1] + $fields[2]) / 2);
	my$id_a = $fields[3];

	my$pos_b = ceil(($fields[5] + $fields[6]) / 2);
	my$id_b = $fields[7];

	next if exists $peak_tracking{$id_a};
	next if exists $peak_tracking{$id_b};

	my$new_pos = ceil(($pos_a + $pos_b) / 2);
	my$peak = $chr.":".$new_pos;

	push(@peak_list, $peak);
	
	$peak_tracking{$id_a} = 1;
	$peak_tracking{$id_b} = 1;
}

close INT;

unlink "Intersect.bed";

print "Merged ",scalar @peak_list," peaks\n";
print "\n---------------------------\n\n";

my$unique_peaks_a = get_unique_peaks($bed_file_a);
my$unique_peaks_b = get_unique_peaks($bed_file_b);

print "Found ",$unique_peaks_a," unique peaks in $bed_file_a\n";
print "Found ",$unique_peaks_b," unique peaks in $bed_file_b\n";
print "\n---------------------------\n\n";

unlink "$bed_file_a.extend";
unlink "$bed_file_b.extend";

###########################################################################################################################################################################

my$final_ext = ceil($final_length / 2);

open(UNSORT, ">", "$outfile.unsorted");

my$peak_count = 0;
foreach my$peak (@peak_list){
	my($chr, $pos) = split(":", $peak);

	my$start = $pos - $final_ext;
	my$end = $pos + $final_ext;

	if($start > 0){
		print UNSORT $chr,"\t",$start,"\t",$end,"\n";
		$peak_count++;
	}
}

close UNSORT;

print "Wrote $peak_count peaks to file\n";
print "\n---------------------------\n\n";

system "bedtools sort -i $outfile.unsorted > $outfile";
unlink "$outfile.unsorted";

###########################################################################################################################################################################

exit;

###########################################################################################################################################################################

#### SUBROUTINES ####

sub prep_bed_file{
	my($file, $extend) = @_;

	open(BED, $file) || die "Cannot open file $file : $!.\n";
	open(TEMP, ">", "$file.extend");

	my$i = 0;
	while(my$line = <BED>){
		next unless $line =~ m/^chr/;
		chomp $line;
		my@fields = split("\t", $line);

		my$chr = $fields[0];
		my$pos = ceil(($fields[1] + $fields[2]) / 2);

		my$start = $pos - $extend;
		my$end = $pos + $extend;

		my$id = $file."_".$i;

		if($start > 0){
			print TEMP $chr,"\t",$start,"\t",$end,"\t",$id,"\n";
			$i++;
		}
	}

	close TEMP;
	return $i;
}

sub get_unique_peaks{
	my$file = shift;

	open(IN, "$file.extend") || die "Cannot open file $file.extend : $!.\n";

	my$i = 0;
	while(my$line = <IN>){
		chomp $line;
		my($chr, $start, $end, $id) = split("\t", $line);

		if(not exists $peak_tracking{$id}){
			my$pos = ceil(($start + $end) / 2);
			my$peak = $chr.":".$pos;
			push(@peak_list, $peak);
			$i++;
		}
	}

	return $i;
}

###########################################################################################################################################################################

__END__

=head1 DESCRIPTION

Create a BED union from two BED files.

=head1 ARGUMENTS

=over 8

=item B<-a>
BED file A

=item B<-b>
BED file B

=item B<-dist>
Maximum distance between summits to merge. Default is 200bp

=item B<-length>
Length to set the merged peaks to. Default is 400bp.

=item B<-outfile>
Output file. 

=item B<-help>
Print help and exit

=back
