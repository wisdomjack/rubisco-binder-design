library(ggmsa)
library(ggplot2)
library(Biostrings)
library(patchwork)

setwd("/Users/jackarcher/Desktop/plos_one_comp_bio/consensus/naive")

seqs <- readAAStringSet("rbcs_aligned.fasta")

# --- Alignment panel ---
p_aln <- ggmsa(seqs, start = 1, end = 134,
               color = "Chemistry_AA",
               seq_name = TRUE) +
  facet_msa(field = 60) +
  geom_msaBar()

# --- Entropy panel ---
seqs_mat <- do.call(rbind, strsplit(as.character(seqs), ""))

shannon_entropy <- apply(seqs_mat, 2, function(col) {
  col <- col[col != "-"]
  if (length(col) == 0) return(0)
  freq <- table(col) / length(col)
  -sum(freq * log2(freq))
})

entropy_df <- data.frame(
  position = 1:length(shannon_entropy),
  entropy  = shannon_entropy
)

p_entropy <- ggplot(entropy_df, aes(x = position, y = entropy)) +
  geom_col(fill = "steelblue", width = 0.8) +
  labs(x = "Position", y = "Shannon entropy (bits)") +
  theme_classic()

# --- Combine and save ---
p_combined <- wrap_plots(p_aln, p_entropy, ncol = 1, heights = c(4, 1))

ggsave("rubisco_alignment.pdf", p_combined, width = 14, height = 10)