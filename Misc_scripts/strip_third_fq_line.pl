#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

##########################################################################################################################
#                                                                                                                        #
# Stip characters following the "+" sign on the third line for each fastQ entry.                                         #
# This information is redundant, and is useless anyway... just takes up space.                                           #
#                                                                                                                        #
##########################################################################################################################

my$infile;
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "help|?" => \$help);

# If -h (or --help) is specified, print usage and exit.
pod2usage(-verbose => 99) if $help;

# Check command line arguments are provided.
my$args = 1;
if(not defined $infile){
	print STDERR "\nMissing argument -i";
	$args = 0;
}

if(not defined $outfile){
	print STDERR "\nMissing argument -o";
	$args = 0;
}

# If any command line argument is missing, $args is set to 0.
# In this case, print usage and exit.
if($args == 0){
	print "\n\n"; # Some formatting.
	pod2usage(-verbose => 1);
}

##########################################################################################################################

open(IN, $infile) || die "Cannot open $infile: $!.\n";
open(OUT, ">", $outfile);

my$i = 0; # Counter to keep track of which line of the fastQ entry we are on. 
while(my$line = <IN>){
	$i++;
	if($i == 3){
		print OUT "+\n";
	}
	else{
		if($i == 4){
			$i = 0; # Reset counter. 
		}
	
		print OUT $line;
	}
}

close IN;
close OUT;

##########################################################################################################################

exit;

##########################################################################################################################

__END__

=head1 DESCRIPTION

Strip characters that follow the "+" sign in the third line of each fastQ entry.

=head1 AUTHOR

Peter Keane

=head1 CONTACT

peterakeane@gmail.com

=head1 USAGE

=over 8

=item B<-i>
Input file (fastQ).

=item B<-o>
Output file (fastQ).

=item B<-h>
Print this helpful help message.

=back