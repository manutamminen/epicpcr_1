
library(tidyverse)
library(ape)


bact_abund_tips <-
  # read.tree("../../data/final/16S.tre") %>%
  read.tree(snakemake@input[[1]]) %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


bact_tre <-
  # read.tree("../../data/final/16S.tre") %>%
  read.tree(snakemake@input[[1]]) %>%
  keep.tip(bact_abund_tips) %>%
  root(outgroup = "JQ837894.1.1415",
       resolve.root = TRUE)


euk_abund_tips <-
  # read.tree("../../data/final/18S.tre") %>%
  read.tree(snakemake@input[[2]]) %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


euk_tre <-
  # read.tree("../../data/final/18S.tre") %>%
  read.tree(snakemake@input[[2]]) %>%
    keep.tip(euk_abund_tips) %>%
    root(outgroup = "Human_18S_rRNA_gene",
         resolve.root = TRUE)


is_tip <- function(tree) tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- function(tree) tree$edge[is_tip(tree), 2]
tip_coord <- function(tree) {
    tre_len <- length(tree$tip.label)
    tibble(Tip = tree$tip.label[ordered_tips(tree)]) %>%
        mutate(Ix = seq(0, 1, 1 / (tre_len - 1)))
}


bact_tip_coords <-
  tip_coord(bact_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  select(-Abundance)


euk_tip_coords <-
  tip_coord(euk_tre) %>%
  separate(Tip, c("Taxonomy", "Abundance"), sep="____") %>%
  select(-Abundance)


bact_tip_percs <-
  # read_tsv("../../tables/16S_abunds.txt") %>%
  read_tsv(snakemake@input[[3]]) %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(bact_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))


euk_tip_percs <-
  # read_tsv("../../tables/18S_abunds.txt") %>%
  read_tsv(snakemake@input[[4]]) %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(euk_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))


plot_hist <- function(sample, tip_percs) {
  plot(NULL, xlim = c(0, 5e3), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
  tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(0, .y), c(.x, .x))))
}


samples <-
  c("Sample1", "Sample2", "Sample3", "Sample4",
  "Fluosample1", "Fluosample2", "Fluosample3", "Fluosample4")



pdf(snakemake@output[[1]])

lmat <- matrix(1:10, ncol = 10)
layout(lmat, widths = c(2, 1, 1, 1, 1, 1, 1, 1, 1, 5), heights = 1)
par(mar=c(1, 1, 1, 1))

plot(ladderize(bact_tre), cex = 0.4,
     align.tip.label = TRUE, show.tip.label = FALSE)

walk(samples, ~plot_hist(., tip_percs=bact_tip_percs))

plot(NULL, xlim = c(0, 5e3), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
bact_tip_coords %>%
  with(walk2(Ix, Taxonomy, ~ text(0, .x, .y, cex=0.3)))

dev.off()



# pdf("../../figures/euk_abunds.pdf")
pdf(snakemake@output[[2]])

lmat <- matrix(1:10, ncol = 10)
layout(lmat, widths = c(2, 1, 1, 1, 1, 1, 1, 1, 1, 5), heights = 1)
par(mar=c(1, 1, 1, 1))

plot(ladderize(euk_tre), cex = 0.4,
     align.tip.label = TRUE, show.tip.label = FALSE)

walk(samples, ~plot_hist(., tip_percs=euk_tip_percs))

plot(NULL, xlim = c(0, 5e3), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
euk_tip_coords %>%
  with(walk2(Ix, Taxonomy, ~ text(0, .x, .y, cex=0.6)))

dev.off()



