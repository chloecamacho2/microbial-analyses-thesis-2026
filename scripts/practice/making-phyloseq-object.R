#installing and loading packages & libraries
library("phyloseq")
BiocManager::install("phyloseq")
library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
theme_set(theme_bw())
library("dplyr")
library("tidyverse")

#loading in example data
data(GlobalPatterns)
data(esophagus)
data(enterotype)
data(soilrep)

example(enterotype, ask=FALSE)
?phyloseq

# Create a pretend OTU table that you read from a file, called otumat
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat
#assigning sample and OTU names
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat
#making taxonomy table
taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

#looking at structure of objects 
class(otumat)
class(taxmat)

#combining matrices into a phyloseq object
library(phyloseq)
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX

physeq = phyloseq(OTU,TAX)
physeq
plot_bar(physeq, fill = "Family")

#creating random data to add to dataset
sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

#create a random phylogenetic tree with the ape package, and add it to your dataset.
install.packages("ape")
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

#combining phylo tree with phyloseq
#1
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1
#2
physeq2 = phyloseq(OTU, TAX, sampledata, random_tree)
physeq2
#double check
identical(physeq1, physeq2)

#making tree plots
plot_tree(physeq1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

plot_tree(physeq1, color="Depth", shape="Location", label.tips="taxa_names", ladderize="right", plot.margin=0.3)

#heat maps
plot_heatmap(physeq1)
plot_heatmap(physeq1, taxa.label="Phylum")
