# Bulk RNA-seq Expression Mini-Analysis (Healthy vs TB)


1.	Inspect shape, row/col names of counts and str(samples).
2.	Compute per-gene mean and per-gene sd (standard deviation) across samples (use apply()), and report the top 2 most variable genes.
3.	Compute per-sample totals (use colSums); which sample has the highest library size?
4.	Create a gene–by–condition mean table for each gene: mean in Healthy vs TB (use tapply() per gene).
5.	List genes upregulated in TB (TB mean > Healthy mean).
6.	Among upregulated genes, which has the largest fold change proxy (TB mean / Healthy mean)? (handle division by zero carefully.)
7.	 Label samples "HighLoad" if their total counts > overall median total; else "Normal". Add this as a new column to samples.
8.	Use table (condition, HighLoad) to summarize; briefly interpret.
9.	Create a row-standardized version of counts (subtract row mean), explain what this helps reveal (in report).
10.	Build a list  called rna_seq  containing counts and samples. Extract the TB mean for gene G3 using nested indexing.
11.	Identify the gene with the highest single value in counts.
12.	Compute per-condition total counts by summing samples in each group. Which group has higher total expression overall?
13.	Use logical indexing to zero out any count values < 5 (simulate low-count filtering). How many entries changed?
14.	Recompute per-gene mean after filtering; did the top mean gene change?
15.	Create a factor of sample order by condition (Healthy first, then TB). Reorder columns of counts accordingly and show the new column order.
16.	Using apply(), compute each gene’s max–min range; which gene is most dynamic?
17.	Create a small barplot of per-sample totals .
18.	Rename gene G2 to G2A in counts (row name operation).
19.	Remove the lowest-mean gene from counts and justify (in report) a scenario where that might be useful
20.	Write a 1–2 paragraph conclusion summarizing which genes differ most by condition and which samples look high-load.
