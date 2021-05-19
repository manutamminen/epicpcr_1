
library(tidyverse)
library(ape)


bact_abund_tips <-
  read.tree(snakemake@input[[1]]) %>%
  # read.tree("../../data/final/16S.tre") %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


bact_tre <-
    read.tree(snakemake@input[[1]]) %>%
    # read.tree("../../data/final/16S.tre") %>%
    keep.tip(bact_abund_tips) %>%
    root(outgroup = "JQ837894.1.1415",
         resolve.root = TRUE)


euk_abund_tips <-
  read.tree(snakemake@input[[2]]) %>%
  # read.tree("../../data/final/18S.tre") %>%
  .$tip.label %>%
  tibble %>%
  separate(".", c("Taxonomy", "Abundance"), sep="____", remove=FALSE) %>%
  mutate(Abundance = as.numeric(Abundance),
         Abundance = ifelse(is.na(Abundance), 51, Abundance)) %>%
  filter(Abundance > 50) %>%
  pull(".")


euk_tre <-
    read.tree(snakemake@input[[2]]) %>%
    # read.tree("../../data/final/18S.tre") %>%
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


connections <-
    read_tsv(snakemake@input[[3]]) %>%
    # read_tsv("../../tables/nonmock_euk_bact_connections.txt") %>%
    mutate(Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, ":", "_"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\(", "__"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\),", "__"),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, "\\)", ""),
           Eukaryotic_taxonomy = str_replace_all(Eukaryotic_taxonomy, " ", "_"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, ":", "_"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\(", "__"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\),", "__"),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, "\\)", ""),
           Bacterial_taxonomy = str_replace_all(Bacterial_taxonomy, " ", "_")) %>%
    inner_join(bact_tip_coords, by = c("Bacterial_taxonomy" = "Taxonomy")) %>%
    inner_join(euk_tip_coords, by = c("Eukaryotic_taxonomy" = "Taxonomy"))


bact_tip_percs <-
  read_tsv(snakemake@input[[4]]) %>%
  # read_tsv("../../tables/16S_abunds.txt") %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(bact_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))


euk_tip_percs <-
  read_tsv(snakemake@input[[5]]) %>%
  # read_tsv("../../tables/18S_abunds.txt") %>%
  mutate(Taxonomy = str_replace_all(Taxonomy, ":", "_"),
         Taxonomy = str_replace_all(Taxonomy, "\\(", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\),", "__"),
         Taxonomy = str_replace_all(Taxonomy, "\\)", ""),
         Taxonomy = str_replace_all(Taxonomy, " ", "_")) %>%
  inner_join(euk_tip_coords, by="Taxonomy") %>%
  group_by(Sample) %>%
  mutate(Perc = Count / sum(Count))


# samples <-
#   c("Fluosample1", "Fluosample2", "Fluosample3", "Fluosample4",
#     "Sample1", "Sample2", "Sample3", "Sample4")

mineral_samples <-
  c("Sample1", "Sample2", "Sample3", "Sample4")

fluor_samples <-
  c("Fluosample1", "Fluosample2", "Fluosample3", "Fluosample4")

plot_bact_hist <- function(sample) {
  plot(NULL, xlim = c(0, 1e4), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
  bact_tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(0, .y), c(.x, .x))))
}

plot_euk_hist <- function(sample) {
  plot(NULL, xlim = c(-1e4, 0), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)
  euk_tip_percs %>%
    filter(Sample == sample) %>%
    with(walk2(Ix, Count, ~ lines(c(1, (1 - .y)), c(.x, .x))))
}

# pdf(snakemake@output[[1]])
png(snakemake@output[[1]], units="in", width=5, height=5, res=300)

lmat <- matrix(1:11, ncol = 11)
layout(lmat, widths = c(2, rep(1, 4), 3, rep(1, 4), 2), heights = 1)
par(mar=c(1, 1, 1, 1))

plot(ladderize(bact_tre), cex = 0.4,
     align.tip.label = TRUE, show.tip.label = FALSE)

walk(mineral_samples, plot_bact_hist)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

connections %>%
  filter(Sample == "Sample1") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("green", 0.05))))

connections %>%
  filter(Sample == "Sample2") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("red", 0.05))))

connections %>%
  filter(Sample == "Sample3") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("black", 0.05))))

connections %>%
  filter(Sample == "Sample4") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("orange", 0.05))))

walk(rev(mineral_samples), plot_euk_hist)

plot(ladderize(euk_tre), align.tip.label = TRUE,
     direction = "leftwards", show.tip.label = FALSE)

dev.off()


# pdf(snakemake@output[[2]])
png(snakemake@output[[2]], units="in", width=5, height=5, res=300)

lmat <- matrix(1:11, ncol = 11)
layout(lmat, widths = c(2, rep(1, 4), 3, rep(1, 4), 2), heights = 1)
par(mar=c(1, 1, 1, 1))

plot(ladderize(bact_tre), cex = 0.4,
     align.tip.label = TRUE, show.tip.label = FALSE)

walk(fluor_samples, plot_bact_hist)

plot(NULL, xlim = c(0, 1), ylim = c(0, 1), type='n', axes=FALSE, ann=FALSE)

connections %>%
  filter(Sample == "Fluosample1") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("green", 0.05))))

connections %>%
  filter(Sample == "Fluosample2") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("red", 0.05))))

connections %>%
  filter(Sample == "Fluosample3") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("black", 0.05))))

connections %>%
  filter(Sample == "Fluosample4") %>%
  with(pwalk(list(Ix.x, Ix.y),
             ~ lines(c(0, 1), c(.x, .y), col = alpha("orange", 0.05))))

walk(rev(fluor_samples), plot_euk_hist)

plot(ladderize(euk_tre), align.tip.label = TRUE,
     direction = "leftwards", show.tip.label = FALSE)

dev.off()

