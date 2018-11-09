# Calculate splicing entropy from a matrix of FPKM values. 

# Input is a tab delimited file, with isoform IDs as row names, and gene IDs as the first column.
# Example:

# isoform_id	gene_id	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
# NM_001130164	Oxr1	1.30145	1.60133	1.97587	1.93863	2.24592	0.918745	3.718	3.48884
# NM_001130163	Oxr1	0.0037831	2.67303	3.52705	3.0456	0.414477	4.89627	3.49248	4.33392
# NM_001130165	Oxr1	4.53812	2.85613	2.66847	3.05996	3.49403	0.875206	4.46073	3.45206
# NM_130885	Oxr1	0.462999	0.00103016	0.00771421	0.309547	0.858262	0.000160108	0.330945	0.000159017
# NM_001130166	Oxr1	0.353866	0.682494	0.973115	0.7274	0.381587	0.554147	0.674204	1.47208

# Author: Peter Keane
# Contact: peterakeane@gmail.com

#################################################################################################################################

# Get command line arguments.
args<- commandArgs(TRUE)
infile<- args[1]
outfile<- args[2]

# Read expression matrix from file. 
# First column of file should correspond to isoform IDs.
exprs<- read.delim(infile, header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)

# Remove genes with only 1 isoform.
gene_ids<- exprs$gene_id
gene_counts<- table(gene_ids) # Count number of occurence of each gene in gene ID column.
gene_ids<- names(gene_counts)[gene_counts > 1]

# Define splicing entropy function.
splice_entropy<- function(gene = NULL, gene_exprs = NULL){
  # Get isoform expression values for given gene.
  iso_exprs<- gene_exprs[gene_exprs$gene_id == gene,]
  
  # Remove column with gene IDs.
  iso_exprs<- iso_exprs[,-1]
  
  # Add some small constant to expression values.
  # This avoids getting a possible log of 0 later on.
  iso_exprs<- iso_exprs + 0.00001
  
  # Calculate p, the proportion that each isoform 
  # contribute to the overall expression of its gene.
  p<- apply(iso_exprs, 2, function(x){x/sum(x)})
  
  # Calculate splicing entropy.
  h<- -1 * colSums(p * log2(p))
  
  return(h)
}

# Apply splicing entropy function.
entropy<- t(sapply(gene_ids, splice_entropy, gene_exprs = exprs))

# Write results to file.
write.table(entropy, file = outfile, sep = '\t', quote = FALSE)

#################################################################################################################################