#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

##############################################################################################################################################################

my$chr_to_remove;
my$help = 0;

GetOptions("chr=s" => \$chr_to_remove, "h|?" => \$help);

pod2usage(-verbose => 99) if $help;

##############################################################################################################################################################

while(my$line = <STDIN>){
	chomp $line;
	
	if($line =~ m/^\@SQ/){
		my($chr_id) = $line =~ m/SN:(chr.+)\s/;
		next if $chr_id eq $chr_to_remove;
		print $line,"\n";
	}
	elsif($line !~ m/^\@/){
		my@fields = split("\t", $line);

		my$chr_a = $fields[2];
		my$chr_b = $fields[6];

		next if $chr_a eq $chr_to_remove;
		next if $chr_b eq $chr_to_remove;

		print $line,"\n";
	}
	else{
		print $line,"\n";
	}
}

##############################################################################################################################################################

exit;

##############################################################################################################################################################

__END__

=head1 DESCRIPTION

Remove a specified chromosome from a SAM file.

=head1 USAGE

cat myAlignment.sam | perl ./filter_chr_SAM.pl -chr chrM

=head1 ARGUMENTS

=over 8

=item B<-chr>

Name of chromosome to remove

=item B<-h>

Print help and exit

=back


























