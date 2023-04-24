library(RegressHaplo)
setwd("/Volumes/Temp/Lukas/LCMV/Run_0355/RegressHaplo")
out_dir <- "output2/"
bam_file <- "/Volumes/Temp/Lukas/LCMV/Run_0355/Single_Chrom/BSF_0355_S01_L.bam"
haps <- readDNAStringSet("output/final_haplo.fasta")
haps
# determine the positions on the reference that were considered variable 
get_variable_positions.pipeline(out_dir)
# get haplotype information
info <- get_fasta.pipeline(out_dir)
# info is a list containing the elements haplotypes and freq
# freq gives the frequency of the reconstructed haplotypes
info$freq
# info$haplotype is a character vector containing the haplotypes.
class(info$haplotypes)
length(info$haplotypes)
# we can use Biostrings to see the haplotypes since the sequences are rather long
DNAStringSet(info$haplotypes)
df <- get_solutions_summary.pipeline(out_dir)
# let's look at every 100th solution
df[seq(from=100,to=700,by=100),]
plot(df$K, df$fit, xlab="K", ylab="fit")
parse_solutions.RegressHaploSolutions("output/solutions.csv")

bam_to_variant_calls.pipeline(bam_file, out_dir, start_pos = 1,end_pos = 500, sig = 0.05)

variant_calls_to_read_table.pipeline(bam_file, out_dir, sig=0.05)
read_table_to_loci.pipeline(out_dir,max_num_haplotypes = 500,min_cover = 100) #outputs loci.csv.
loci_to_haplotypes.pipeline(out_dir,max_num_haplotypes = 500) #outputs h.csv.
haplotypes_to_parameters.pipeline(out_dir) #outputs y.csv, P.csv
parameters_to_solutions.pipeline(out_dir, num_trials=700,rho_vals = c(0.1, 1, 5, 10, 20)) #outputs solutions.csv
#parameters_to_solutions.pipeline(out_dir, num_trials=700,rho_vals = c(0.01,1.0,2.5)) #outputs solutions.csv
solutions_to_haplotypes.pipeline(out_dir) #final_haplo.csv
haplotypes_to_fasta.pipeline(bam_file, out_dir) #final_haplo.fasta



