
# Created on 02/10/2020

# compute a Fisher's exact test from 3 lists: setA, setB, and setU=universe

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript fisher_test_from_lists.R <setA> <setB> <setU> <out_prefix>\n\n")
	cat("each set should be specified as a file, with one element per line\n")
	cat("U is the universe and should contain both set A and set B\n\n")
	q()
}

SET_A <- args[1]
SET_B <- args[2]
SET_U <- args[3]
OUT_PREFIX <- args[4]

set_A <- as.character(drop(as.matrix(read.table(SET_A))))
set_B <- as.character(drop(as.matrix(read.table(SET_B))))
set_U <- as.character(drop(as.matrix(read.table(SET_U))))

if (length(set_A %in% set_U) < length(set_A))
	write("Set A is not fully contained in the universe set U", stderr()) 

if (length(set_B %in% set_U) < length(set_B))
	write("Set B is not fully contained in the universe set U", stderr()) 

INTERS_FILE <- paste(OUT_PREFIX, "intersection.txt", sep = "_")
FISHER_FILE <- paste(OUT_PREFIX, "fisher.txt", sep = "_")

intersection <- set_A[set_A %in% set_B]
write.table(intersection, file = INTERS_FILE, quote = FALSE, col.names = FALSE, row.names = FALSE)

# run the Fisher tests

a_and_b <- length(intersection)
a_not_b <- length(set_A) - length(intersection)
b_not_a <- length(set_B) - length(intersection)
not_a_not_b <- length(set_U) - a_and_b - a_not_b - b_not_a
F <- matrix(c(a_and_b, b_not_a, a_not_b, not_a_not_b), nrow=2, ncol=2)

sink(file = FISHER_FILE)
write.table(F)
fisher.test(x = F, alternative = "greater")
sink()


q()


