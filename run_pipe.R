library(RegressHaplo)
args = commandArgs(trailingOnly=TRUE)
bam_file <- args[1]
out_dir <- args[2]
start_pos <- args[3]
end_pos <- args[4]
dir.create(out_dir)
full_pipeline(bam_file, out_dir, start_pos=start_pos, end_pos=end_pos)
#full_pipeline(bam_file, out_dir, start_pos=1, end_pos=500, rho = 0.184942158)

