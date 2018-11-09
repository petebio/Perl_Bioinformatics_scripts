#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

#################################################################################################################################################################### 

# Name: gtf_to_gff3.pl
# Description: Convert a GTF file (as produced by genePredToGtf from UCSC) to the GFF3 format.
# Author: Peter Keane (peterakeane@gmail.com).
# Last modified: 29th October 2016.
# Version: 1.0

# Example of GTF input format (tab delimited):

# chr1	refGene	exon	975199	976269	.	-	.	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "1"; exon_id "NM_001291366.1"; gene_name "PERM1";
# chr1	refGene	CDS	976175	976269	.	-	2	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "1"; exon_id "NM_001291366.1"; gene_name "PERM1";
# chr1	refGene	exon	976499	976624	.	-	.	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "2"; exon_id "NM_001291366.2"; gene_name "PERM1";
# chr1	refGene	CDS	976499	976624	.	-	2	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "2"; exon_id "NM_001291366.2"; gene_name "PERM1";
# chr1	refGene	exon	978881	981047	.	-	.	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "3"; exon_id "NM_001291366.3"; gene_name "PERM1";
# chr1	refGene	CDS	978881	981029	.	-	0	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "3"; exon_id "NM_001291366.3"; gene_name "PERM1";
# chr1	refGene	exon	982065	982117	.	-	.	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "4"; exon_id "NM_001291366.4"; gene_name "PERM1";
# chr1	refGene	start_codon	981027	981029	.	-	0	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "1"; exon_id "NM_001291366.1"; gene_name "PERM1";
# chr1	refGene	stop_codon	976172	976174	.	-	0	gene_id "PERM1"; transcript_id "NM_001291366"; exon_number "1"; exon_id "NM_001291366.1"; gene_name "PERM1";

# Example of GFF3 output format (tab delimited):

# chr1	refGene	mRNA	975199	982117	.	-	.	ID=NM_001291366;Name=NM_001291366;Parent=PERM1
# chr1	refGene	five_prime_UTR	982065	982117	.	-	.	ID=five_prime_UTR:NM_001291366.4;Parent=NM_001291366
# chr1	refGene	five_prime_UTR	981027	981047	.	-	.	ID=five_prime_UTR:NM_001291366.3;Parent=NM_001291366
# chr1	refGene	start_codon	981027	981029	.	-	0	ID=start_codon:NM_001291366.3;Parent=NM_001291366
# chr1	refGene	exon	975199	976269	.	-	.	ID=exon:NM_001291366.1;Parent=NM_001291366
# chr1	refGene	exon	976499	976624	.	-	.	ID=exon:NM_001291366.2;Parent=NM_001291366
# chr1	refGene	exon	978881	981047	.	-	.	ID=exon:NM_001291366.3;Parent=NM_001291366
# chr1	refGene	exon	982065	982117	.	-	.	ID=exon:NM_001291366.4;Parent=NM_001291366
# chr1	refGene	CDS	976175	976269	.	-	2	ID=CDS:NM_001291366.1;Parent=NM_001291366
# chr1	refGene	CDS	976499	976624	.	-	2	ID=CDS:NM_001291366.2;Parent=NM_001291366
# chr1	refGene	CDS	978881	981029	.	-	0	ID=CDS:NM_001291366.3;Parent=NM_001291366
# chr1	refGene	stop_codon	976172	976174	.	-	0	ID=stop_codon:NM_001291366.1;Parent=NM_001291366
# chr1	refGene	three_prime_UTR	975199	976172	.	-	.	ID=three_prime_UTR:NM_001291366.1;Parent=NM_001291366

#################################################################################################################################################################### 

my$infile;
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$infile, "outfile=s" => \$outfile, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check command line options.
# If any are missing, set $args to 0 and print usage information.  
my$args = 1;
if(not defined $infile){
	print STDERR "\nMissing argument: -i";
	$args = 0;
}

if(not defined $outfile){
	print STDERR "\nMissing argument: -o";
	$args = 0;
}

if($args == 0){
	print "\n\n"; # Some formatting.
	pod2usage(-verbose => 1);
}

#################################################################################################################################################################### 

# Read GTF file and extract information requires for the GFF3 format.

open(GTF, $infile) || die "Cannot open file $infile: $!.\n";
print "Reading GTF file.\n";

# Declare some variables.
my%gene_tracking;
my%transcript_tracking;
my%exon_tracking;
my%cds_tracking;
my%start_codon_tracking;
my%stop_codon_tracking;
my%gene_transcript_map;
my%start_exon;
my%end_exon;

while(my$line = <GTF>){
	chomp $line;
	my@fields = split("\t", $line);
	
	my$chr = $fields[0];
	my$source = $fields[1];
	my$feature_type = $fields[2];
	my$start = $fields[3];
	my$end = $fields[4];
	my$strand = $fields[6];
	my$phase = $fields[7];
	my$attributes = $fields[8];
	
	# Skip entries from non-standard chromosomes.
	# Should fix this later, but this is OK for now. 
	next if $chr !~ m/^chr(:?\d+|M|X|Y)$/;
	
	# Ensure start/end values are in correct order. 
	if($start > $end){
		my$temp = $start;
		$start = $end;
		$end = $start;
	}
	
	# Extract necessary info from attributes
	my($gene_id, $transcript_id, $exon_id) = $attributes =~ m/^gene_id\s+"([^"]+)";\s+transcript_id\s+"([^"]+)";\s+exon_number\s+"(\d+)";/;
	
	# Find the start and end positions for the gene.
	# If no entry for the current gene exists in %gene_tracking, initialise it.
	# In subsequent lines of the GTF, check if the current exon start position < gene start position,
	# and similarly if exon end position > gene end position.
	# If so, replace gene start/end position with the relavant value of current exon.
	if(not exists $gene_tracking{$gene_id}){
		$gene_tracking{$gene_id} = [($start, $end, $chr, $strand, $source)];
	}
	else{
		if($start < $gene_tracking{$gene_id}[0]){
			$gene_tracking{$gene_id}[0] = $start;
		}
		
		if($end > $gene_tracking{$gene_id}[1]){
			$gene_tracking{$gene_id}[1] = $end;
		}
	}
	
	# Do the same for the transcript.
	if(not exists $transcript_tracking{$transcript_id}){
		$transcript_tracking{$transcript_id} = [($start, $end, $strand)];
	}
	else{
		if($start < $transcript_tracking{$transcript_id}[0]){
			$transcript_tracking{$transcript_id}[0] = $start;
		}
		
		if($end > $transcript_tracking{$transcript_id}[1]){
			$transcript_tracking{$transcript_id}[1] = $end;
		}
	}
	
	# Record start/end/phase information for exon + CDS, as well as start and stop codon positions. 
	if($feature_type eq "exon"){
		$exon_tracking{$transcript_id}{$exon_id} = [($start, $end)];
	}
	elsif($feature_type eq "CDS"){
		$cds_tracking{$transcript_id}{$exon_id} = [($start, $end, $phase)];
		
		if(not exists $start_exon{$transcript_id}){
			$start_exon{$transcript_id} = $exon_id;
			$end_exon{$transcript_id} = $exon_id;
		}
		else{
			if($strand eq "+"){
				if($exon_id < $start_exon{$transcript_id}){
					$start_exon{$transcript_id} = $exon_id;
				}
	
				if($exon_id > $end_exon{$transcript_id}){
					$end_exon{$transcript_id} = $exon_id;
				}
			}
			else{
				if($exon_id > $start_exon{$transcript_id}){
					$start_exon{$transcript_id} = $exon_id;
				}

				if($exon_id < $end_exon{$transcript_id}){
					$end_exon{$transcript_id} = $exon_id;
				}
			}
		}
	}
	elsif($feature_type eq "start_codon"){
		$start_codon_tracking{$transcript_id} = [($start, $end, $phase)];
	}
	elsif($feature_type eq "stop_codon"){
		$stop_codon_tracking{$transcript_id} = [($start, $end, $phase)];
	}
	
	# Record gene to transcript mapping.
	$gene_transcript_map{$gene_id}{$transcript_id} = 1;
}

close GTF;

print "Found ",scalar keys %gene_tracking," genes with ",scalar keys %transcript_tracking," transcripts.\n";

#################################################################################################################################################################### 

# Find the position of the 5' UTR for each transcript.
# i.e. exonic regions that are upstream of the transcripts start codon.
# This may include entire exons, or only regions of exons. 

print "Searching for 5' UTRs.\n";

my%five_utr_tracking;
foreach my$isoform (keys %start_codon_tracking){
	my$strand = $transcript_tracking{$isoform}[2];
	my($codon_start, $codon_end, $phase) = @{$start_codon_tracking{$isoform}};

	if($strand eq "+"){
		foreach my$exon (sort {$a <=> $b} keys %{$exon_tracking{$isoform}}){
			my($upstream, $downstream, $exon_phase) = @{$exon_tracking{$isoform}{$exon}};
			
			if($downstream < $codon_end){
				push(@{$five_utr_tracking{$isoform}}, [($exon, $upstream, $downstream)]);	
			}
			elsif($upstream < $codon_start && $downstream > $codon_end){
				push(@{$five_utr_tracking{$isoform}}, [($exon, $upstream, $codon_end)]);
			}
		}
	}
	elsif($strand eq "-"){
		foreach my$exon (sort {$b <=> $a} keys %{$exon_tracking{$isoform}}){
			my($downstream, $upstream, $exon_phase) = @{$exon_tracking{$isoform}{$exon}};

			if($downstream > $codon_start){
				push(@{$five_utr_tracking{$isoform}}, [($exon, $downstream, $upstream)]);
			}
			elsif($upstream > $codon_end && $downstream < $codon_start){
				push(@{$five_utr_tracking{$isoform}}, [($exon, $codon_start, $upstream)]);
			}
		}
	}
}

#################################################################################################################################################################### 

# Find the position of the 3' UTR for each transcript.
# i.e. exonic regions that are downstream of the transcripts stop codon.
# This may include entire exons, or only regions of exons. 

print "Searching for 3' UTRs.\n";

my%three_utr_tracking;
foreach my$isoform (keys %stop_codon_tracking){
	my$strand = $transcript_tracking{$isoform}[2];
	my($codon_start, $codon_end, $codon_phase) = @{$stop_codon_tracking{$isoform}};

	if($strand eq "+"){
		foreach my$exon (sort {$a <=> $b} keys %{$exon_tracking{$isoform}}){
			my($upstream, $downstream, $exon_phase) = @{$exon_tracking{$isoform}{$exon}};

			if($upstream > $codon_end){
				push(@{$three_utr_tracking{$isoform}}, [($exon, $upstream, $downstream)]);
			}
			elsif($upstream < $codon_start && $downstream > $codon_end){
				push(@{$three_utr_tracking{$isoform}}, [($exon, $codon_start, $downstream)]);
			}
		}
	} 
	elsif($strand eq "-"){
		foreach my$exon (sort {$b <=> $a} keys %{$exon_tracking{$isoform}}){
			my($downstream, $upstream, $exon_phase) = @{$exon_tracking{$isoform}{$exon}};

			if($upstream < $codon_start){
				push(@{$three_utr_tracking{$isoform}}, [($exon, $downstream, $upstream)]);
			}
			elsif($upstream > $codon_end && $downstream < $codon_start){
				push(@{$three_utr_tracking{$isoform}}, [($exon, $downstream, $codon_start)]);
			}
		}
	}
}

#################################################################################################################################################################### 

# Sort gene start positions and chromosomes.

my%chr_mapping;
foreach my$gene (keys %gene_tracking){
	my$chr = $gene_tracking{$gene}[2];
	my$start = $gene_tracking{$gene}[0];
	$chr_mapping{$chr}{$start} = $gene;
}

#################################################################################################################################################################### 

# Print gene/transcript annotations in GFF3 format. Write to $outfile.

open(GFF, ">", $outfile);
print "Converting to GFF3 format.\n";

print GFF "##gff-version 3\n";

foreach my$cur_chr (sort keys %chr_mapping){
	foreach my$pos (sort {$a <=> $b} keys %{$chr_mapping{$cur_chr}}){
		my$gene = $chr_mapping{$cur_chr}{$pos};
		my($gene_start, $gene_end, $chr, $strand, $source) = @{$gene_tracking{$gene}};
		print GFF $chr,"\t",$source,"\tgene\t",$gene_start,"\t",$gene_end,"\t.\t",$strand,"\t.\tID=",$gene,";Name=",$gene,"\n";


		foreach my$isoform (keys %{$gene_transcript_map{$gene}}){
			my($isoform_start, $isoform_end, $isoform_strand) = @{$transcript_tracking{$isoform}};

			# For protein-coding transcripts the type is mRNA.
			# For non-coding transcripts the type is transcript.
			my$type = "mRNA";
			if(not exists $cds_tracking{$isoform}){
				$type = "transcript";
			}

			print GFF $chr,"\t",$source,"\t",$type,"\t",$isoform_start,"\t",$isoform_end,"\t.\t";
			print GFF $strand,"\t.\tID=",$isoform,";Name=",$isoform,";Parent=",$gene,"\n";

			if(exists $five_utr_tracking{$isoform}){
				foreach my$utr (@{$five_utr_tracking{$isoform}}){
					my$exon = $utr -> [0];
					my$utr_exon = $isoform.".".$exon;
				
					my$utr_start = $utr -> [1];
					my$utr_end = $utr -> [2];

					print GFF $chr,"\t",$source,"\tfive_prime_UTR\t",$utr_start,"\t",$utr_end,"\t.\t";
					print GFF $strand,"\t.\tID=five_prime_UTR:",$utr_exon,";Parent=",$isoform,"\n";
				}

				my($start_codon_start, $start_codon_end, $start_codon_phase) = @{$start_codon_tracking{$isoform}};

				print GFF $chr,"\t",$source,"\tstart_codon\t",$start_codon_start,"\t",$start_codon_end,"\t.\t";
				print GFF $strand,"\t",$start_codon_phase,"\tID=start_codon:";
				print GFF $isoform,".",$start_exon{$isoform},";Parent=",$isoform,"\n";
			}

			foreach my$exon (sort {$a <=> $b} keys %{$exon_tracking{$isoform}}){
				my($exon_start, $exon_stop) = @{$exon_tracking{$isoform}{$exon}};
				my$exon_id = $isoform.".".$exon;

				print GFF $chr,"\t",$source,"\texon\t",$exon_start,"\t",$exon_stop,"\t.\t";
				print GFF $strand,"\t.\tID=exon:",$exon_id,";Parent=",$isoform,"\n";
			}

			if(exists $cds_tracking{$isoform}){
				foreach my$exon (sort {$a <=> $b} keys %{$cds_tracking{$isoform}}){
					my($exon_start, $exon_end, $exon_phase) = @{$cds_tracking{$isoform}{$exon}};
					my$exon_id = $isoform.".".$exon;

					print GFF $chr,"\t",$source,"\tCDS\t",$exon_start,"\t",$exon_end,"\t.\t";
					print GFF $strand,"\t",$exon_phase,"\tID=CDS:",$exon_id,";Parent=",$isoform,"\n";
				}
			}

			if(exists $three_utr_tracking{$isoform}){
				my($stop_codon_start, $stop_codon_end, $stop_codon_phase) = @{$stop_codon_tracking{$isoform}};

				print GFF $chr,"\t",$source,"\tstop_codon\t",$stop_codon_start,"\t",$stop_codon_end,"\t.\t";
				print GFF $strand,"\t",$stop_codon_phase,"\tID=stop_codon:";
				print GFF $isoform,".",$end_exon{$isoform},";Parent=",$isoform,"\n";

				foreach my$utr (@{$three_utr_tracking{$isoform}}){
					my$exon = $utr -> [0];
					my$utr_exon = $isoform.".".$exon;

					my$utr_start = $utr -> [1];
					my$utr_end = $utr -> [2];

					print GFF $chr,"\t",$source,"\tthree_prime_UTR\t",$utr_start,"\t",$utr_end,"\t.\t";
					print GFF $strand,"\t.\tID=three_prime_UTR:",$utr_exon,";Parent=",$isoform,"\n";
				}
			}
		}
	}
}

close GFF;
print "Done!!\n";

#################################################################################################################################################################### 

exit;

#################################################################################################################################################################### 

__END__

=head1 DESCRIPTION

Convert a GTF file (as produced by genePredToGtf from UCSC) to the GFF3 format.

=head1 AUTHOR

Peter Keane

=head1 CONTACT

peterakeane@gmail.com

=head1 USAGE

=over 8

=item B<-i>
Input file (GTF format).

=item B<-o>
Output file.

=item B<-h>
Print this helpful help message.

=back
