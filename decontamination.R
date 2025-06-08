#check for remotes
if(!require(remotes))
  install.packages("remotes")
#install mmgenome2 using remotes
remotes::install_github("kasperskytte/mmgenome2")

setwd("C:\Users\sipe3\OneDrive\Desktop\Speciale\Bioinformatics files")
library("mmgenome2")

#Run everything in this section with one genome at a time
cov <- read.table("EG12_bbmap_coverage_statistics.tsv", sep='\t',
                  head=TRUE,comment.char="")[,c("X.ID", "Avg_fold")]

ass <- Biostrings::readDNAStringSet("EG12_scaffolds.fasta",
                                    format = "fasta")


genome <- mmload(
  assembly = ass,
  coverage = cov)

mmplot(genome,
       x = "gc",
       y = "cov_Avg_fold",
       y_scale = "log10",
       locator = TRUE)

#Removal of Decontamination
selection <- data.frame(gc = c(53.415, 74.38, 69.699, 53.008, 47.106),
                        cov_Avg_fold = c(1476.234, 337.18, 48.702, 48.702, 261.393))


decon_genome <- mmextract(genome, selection)

mmexport(decon_genome,
         assembly=ass,
         file = "EG12_high_GC_decontaminated_genome_assembly.fasta")
