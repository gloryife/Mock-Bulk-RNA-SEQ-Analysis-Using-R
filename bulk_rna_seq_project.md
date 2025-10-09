Gene Expression Profiling to Identify Biomarkers for Tuberculosis
Infection Status using Bulk RNA-Seq: A BioiformHer Mini Project
================
Glory Jayeoba
10th October, 2025

# Background

Tuberculosis infection (TB) remains a major global health concern due to
its high transmissibility and persistent burden across populations. It
is also one of the leading cause of death with an estimated 1.25 million
death worldwide.

RNA-sequencing (RNA-seq) is a powerful technique that can provide
valuable insight into organisms’ biological mechanism and molecular
pathways. Understanding the transcriptomic differences between a
diseased and healthy state in organisms can provide insight to molecular
signatures of diseases that can enhance diagnosis, treatment and disease
monitoring.

In this project, we will be investigating gene expression differences of
a small panel of immune-related genes in TB and Healthy individuals to
understand the immune mechanisms underlying TB infection.

Identifying key transcriptomic differences between TB and Healthy
individuals provides a basis for potential biomarkers for disease
diagnosis and severity. Additionally, distinguishing high-load TB cases
from moderate ones at the transcriptomic level can also provide valuable
insight for disease monitoring and individualized treatment strategies.

# Aim

The aim of this project is to investigate host gene expression changes
in tuberculosis patients compared to healthy individuals and identify
transcriptomic signatures of disease and bacterial load.

# Methodology

In this project, a small RNA-seq dataset was simulated to reflect the
expression of immune-related genes in TB and Healthy individuals. A
total of 6 gene were identified for this study, 4 of which are
immune-related genes while the reminding 2 are housekeeping genes which
were included in the analysis as control.

***This dataset serves as a training resource for learning and
practicing RNA-seq data analysis in R. This project is also in
collaboration with the BioinformHer Initiative***

# Analysis

## Setting up the notebook output

This is an R notebook that will be knitted into a github document. The
setting below therefore ensures that the codes, figures and writings are
displayed in clean and consistent manner in the final knitted document.
The codes suppresses warnings, messages and sets all figures to a
consistent size.

``` r
##setting the final output format
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 8,
  fig.height = 5
)
```

## Loading Dataset (Setting up the R enviroment)

To begin analysis, it is important that the R environment is well set up
and ready for analysis. This may include loading necessary packages,
uploading the required dataset etc.

The data used for analysis in this project include a simulated raw
expression counts of a mock RNA-seq experiment and a metadata.

The datset will be uploaded as a matrix while the metadata will be
uploaded as a dataframe for convenient use and manipulation.

``` r
##creating the dataset as a matrix and assigning it to an object called data.

data <- matrix(c(100, 180, 50, 80, 60, 140, 230, 110, 110, 170, 62, 136, 320, 65, 195,
                 280, 65, 143, 270, 85, 140, 220, 63, 139, 110, 190, 55, 75, 58, 135),
               ncol = 5, 
               nrow = 6,
               dimnames = list(c("STAT1", "IFIT3", "ISG15", "GBP5", "ALAS2", "RPL13A"),
                              paste0("Sample", "_", 1:5)))


##creating a metadata for the dataset
metadata <- data.frame(Sample = paste0("Sample", "_", 1:5),
                       Health_status = c("Healthy", "TB", "TB", "TB", "Healthy"))
```

## 1. Inspecting the dataset

One good practice in R is inspecting your dataset before beginning
analysis. This include inspecting the dimensions of your dataset (the
shape), the structure, completeness (missing values), data type, data
class etc.

In this project, we will be inspecting the shape, the structure, class
and type of data used in this project.

``` r
##inspecting the dimensions, column and row names, type and class of data and structure of data
dim(data); dim(metadata)
```

    ## [1] 6 5

    ## [1] 5 2

``` r
nrow(data); nrow(metadata)
```

    ## [1] 6

    ## [1] 5

``` r
ncol(data); ncol(metadata)
```

    ## [1] 5

    ## [1] 2

``` r
colnames(data); colnames(metadata)
```

    ## [1] "Sample_1" "Sample_2" "Sample_3" "Sample_4" "Sample_5"

    ## [1] "Sample"        "Health_status"

``` r
rownames(data); rownames(metadata)
```

    ## [1] "STAT1"  "IFIT3"  "ISG15"  "GBP5"   "ALAS2"  "RPL13A"

    ## [1] "1" "2" "3" "4" "5"

``` r
typeof(data); typeof(metadata)
```

    ## [1] "double"

    ## [1] "list"

``` r
class(data); class(metadata)
```

    ## [1] "matrix" "array"

    ## [1] "data.frame"

``` r
str(data); str(metadata)
```

    ##  num [1:6, 1:5] 100 180 50 80 60 140 230 110 110 170 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:6] "STAT1" "IFIT3" "ISG15" "GBP5" ...
    ##   ..$ : chr [1:5] "Sample_1" "Sample_2" "Sample_3" "Sample_4" ...

    ## 'data.frame':    5 obs. of  2 variables:
    ##  $ Sample       : chr  "Sample_1" "Sample_2" "Sample_3" "Sample_4" ...
    ##  $ Health_status: chr  "Healthy" "TB" "TB" "TB" ...

The dataset is a matrix with 6 rows and 5 columns. Each row represents a
gene and each column represents a sample (either TB or healthy). Lastly,
the type of data contained in the matrix is said to be a double.

The metadata is a dataframe with a 5 rows and 2 columns. Each row
represents a sample and each each column represents information about
the sample and health status

## 2. Determining the top 2 most variable genes

A quick way to begin analysis and interpretation of your data is to
identify hit genes i.e which showed the most differential expression
between the 2 health status. To do this, we will need to first compute
the mean and standard deviation for each gene across all samples. A high
standard deviation value correspond to a high variance in gene
expression.

For this project, we will be identifying the top 2 most variable gene.

``` r
##computing per-gene mean
gene_mean <- apply(data, 1, mean)    # 1 applies function to all rows of the data

##computing per-gene sd
gene_sd <- apply(data, 1, sd)

##sorting the per-gene mean in decreasing order by the per-gene variance
variable_gene <- sort_by(gene_mean, gene_sd, decreasing = TRUE)

##determining the top 2 most variable gene
names(variable_gene[1:2])
```

    ## [1] "STAT1" "GBP5"

We can observe that the STAT1 gene and GBP5 gene were the 2 most
variable gene. These genes had higher expression counts (high mean) in
TB samples compared with healthy samples.

# 3. Determining the library size of each sample

Another basic analysis is to determine the library size for each sample.
Library size also known as sequencing depth refers to the total number
of reads (expression counts) obtained per sample. To compute the library
size for each sample, you only need to calculate the total expression
count for each sample.

``` r
##computing per-sample sum
sample_sum <- colSums(data)

##determining the sample with the highest sample size
sample_sum[which.max(sample_sum)]
```

    ## Sample_3 
    ##     1068

Sample_3 is shown to be the sample with the highest library size with a
total expression count of 1068.

## 4. Creating a gene-by-condition mean table

When analyzing gene expression data, evaluating the spread of expression
across all samples can obscure meaningful patterns due to sample-level
variability. A good approach to facilitate better visualization and
interpretation is to compare the expression levels for each gene by
biological state (e.g. diseased vs Healthy status).

For this project, we will computing the mean of expression count per
gene by different health status in a table. To do this, we will use
tapply() to compute the mean for each gene and then combine them into a
table

``` r
##computing health condition as factor
factor <- as.factor(metadata$Health_status)

##using tapply to compute the mean expression count of each gene by factor
STAT1 <- tapply(data["STAT1", ], factor, mean)
IFIT3 <- tapply(data["IFIT3", ], factor, mean)
ISG15 <- tapply(data["ISG15", ], factor, mean)
GBP5 <- tapply(data[ "GBP5", ], factor, mean)
ALAS2 <- tapply(data["ALAS2", ], factor, mean)
RPL13A <- tapply(data["RPL13A", ], factor, mean)

##combining all values for each gene to a table
gene_by_con <- rbind(STAT1, IFIT3, ISG15, GBP5, ALAS2, RPL13A)

##rounding up the mean to one decimal place
gene_by_con <- round(gene_by_con, 1)

gene_by_con
```

    ##        Healthy    TB
    ## STAT1    105.0 273.3
    ## IFIT3    185.0  86.7
    ## ISG15     52.5 148.3
    ## GBP5      77.5 223.3
    ## ALAS2     59.0  63.3
    ## RPL13A   137.5 139.3

## 5. Identifying upregulated genes in TB

With the counts group by condition, it will be easier to observe
differences in the expression by condition.

``` r
##listing up-regulated genes in TB
up_gene <- rownames(gene_by_con)[ gene_by_con[, "TB"] > 
                                    gene_by_con[, "Healthy"] ]

up_gene
```

    ## [1] "STAT1"  "ISG15"  "GBP5"   "ALAS2"  "RPL13A"

The result listed 5 out of the 6 genes to be up-regulated genes. All 5
genes listed had higher expression counts in healthy individuals vs TB
individuals. The sixth gene IFIT3 had higher expression counts in
healthy individuals vs TB individuals.

### 6. Determining which gene has the largest fold change proxy

Another good way to identify the up-regulated genes is to determine the
fold change proxy for each gene.

- Gene with fold change proxy with values equal to 1 shows no difference
  between the 2 health status.

- Gene with fold change with values less than 1 shows down-regulated
  genes

- Genes with fold change with values greater than 1 shows up-regulated
  genes

``` r
##dividing the TB mean by Healthy mean 
fold_change <- gene_by_con[, "TB"] / gene_by_con[, "Healthy"]


##top-regulated genes
top_fold_change <- names(fold_change[fold_change < 1 | fold_change > 1.1])

top_fold_change
```

    ## [1] "STAT1" "IFIT3" "ISG15" "GBP5"

The genes identified to be up-regulated in this project are “GBP5”
“ISG15” “STAT1” and “IFIT3”. Identifying such gene can be helpful as
they can serve as potential biomarker for TB infection.

## 7. Determing TB load levels in samples

``` r
##recall that total count was stored in an object called sample_sum
sample_sum
```

    ## Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 
    ##      610      818     1068      917      623

``` r
##computing the median of overall total median per sample
sample_median <- median(sample_sum)

##Identifying sample that is high load and adding it to a new column in metadata
metadata$TB_Load <- ifelse(sample_sum > sample_median, "Highload", "Normal")

##extracting samples that are highload
metadata$Sample[which(metadata$TB_Load == "Highload")]
```

    ## [1] "Sample_3" "Sample_4"

## 8. Interpretation

Highload implies large library size which corresponds to high gene
expression activity. Sample_3 and Sample_4 identified as samples with
the highest load were all TB samples while all healthy samples had
normal expression levels. Since the gene studied are all immune-related
genes, it can be implied that the increased expression of these genes
may have been triggered by the presence of TB infection.

## 9. Creating a row-standardized version of counts

Because we are measuring the expression of several genes in RNA-seq, it
is highly likely for one gene to have very high expression counts and
another gene to have very low expression count. Using raw expression
counts may therefore not be ideal for proper interpreatation because
gene with high expression count may be very distracting and obscure true
pattern. One way to combat this is by normalizing your dataset.

Creating a row standardized version of your counts (subtracting the mean
from each expression count per gene) is a good way to normalized the
expression counts. Applying row standardization to the expression counts
allows you to focus on relative up- or down-regulation across conditions
rather than absolute expression levels.

``` r
##recall that computed mean per gene was stored in a object called gene_mean
gene_mean
```

    ##  STAT1  IFIT3  ISG15   GBP5  ALAS2 RPL13A 
    ##  206.0  126.0  110.0  165.0   61.6  138.6

``` r
##Creating a row-standardized version of counts
data_standardized <- sweep(data, 1, gene_mean, "-")

data_standardized
```

    ##        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
    ## STAT1    -106.0     24.0    114.0     64.0    -96.0
    ## IFIT3      54.0    -16.0    -61.0    -41.0     64.0
    ## ISG15     -60.0      0.0     85.0     30.0    -55.0
    ## GBP5      -85.0      5.0    115.0     55.0    -90.0
    ## ALAS2      -1.6      0.4      3.4      1.4     -3.6
    ## RPL13A      1.4     -2.6      4.4      0.4     -3.6

## Creating a list with dataset

``` r
##creating a list called rna_seq
rna_seq <- list(count = data,
                sample = metadata)


##extracting only TB samples from metadata
TB_samples <- rna_seq$sample$Sample[rna_seq$sample$Health_status == "TB" ]

##extracting the expression counts of only TB samples for gene 3
gene3_TB_counts <- rna_seq$count[3, TB_samples]

##computing the mean for the expression counts
gene3_TB_mean <- round(mean(gene3_TB_counts), 1)

gene3_TB_mean
```

    ## [1] 148.3

## 11. Identifying gene with the highest single value in counts

``` r
##checking which gene has the highest value in count
max_gene <- rownames(which(data == max(data), arr.ind = TRUE))

max_gene
```

    ## [1] "STAT1"

## 12. Computing total expression per group (condition)

``` r
##extracting healthy from metadata (recall TB samples have been extracted previously)

Healthy_samples <- rna_seq$sample$Sample[rna_seq$sample$Health_status == "Healthy" ]

##finding the sum of all expression counts for both group

All_TB_counts <- sum(rna_seq$count[, TB_samples])

All_healthy_counts <- sum(rna_seq$count[, Healthy_samples])

All_healthy_counts ; All_TB_counts
```

    ## [1] 1233

    ## [1] 2803

## 13. Simulating low-count filtering using logical indexing

``` r
##filtering out counts less than 5 and storing in an object low_count
low_count <- data[data < 5]

## creating a filtered dataset where all reads < 5 are changed to 0
data_filtered <- data
data_filtered[low_count] <- 0

##checking how many entries changed
which(data_filtered == 0)
```

    ## integer(0)

## 14 Recomputing per-gene mean after filtering

``` r
##computing per-gene mean for using filtered data
new_gene_mean <- apply(data_filtered, 1, mean)


##checking if the new gene mean is the same as the previous mean
gene_mean == new_gene_mean
```

    ##  STAT1  IFIT3  ISG15   GBP5  ALAS2 RPL13A 
    ##   TRUE   TRUE   TRUE   TRUE   TRUE   TRUE

## 15. Reordering the columns of counts according to condition

``` r
##recall health condition as factor was stored in an object called factor
factor
```

    ## [1] Healthy TB      TB      TB      Healthy
    ## Levels: Healthy TB

``` r
##reordering samples in metadata by condition
order <- metadata$Sample[order(factor)]

##reordering data by new order
new_data <- data[, order]

new_data
```

    ##        Sample_1 Sample_5 Sample_2 Sample_3 Sample_4
    ## STAT1       100      110      230      320      270
    ## IFIT3       180      190      110       65       85
    ## ISG15        50       55      110      195      140
    ## GBP5         80       75      170      280      220
    ## ALAS2        60       58       62       65       63
    ## RPL13A      140      135      136      143      139

## 16. Computing maximum and minimum range of counts for each gene

``` r
##computing range of counts for each gene
range <- apply(data, 1, function(x) {
  max(x) - min(x)
})


##checking which gene is the most dynamic
max(range) ; names(which.max(range))
```

    ## [1] 220

    ## [1] "STAT1"

## 17. Creating a barplot using samples totals

``` r
##saving barplot as pdf
pdf("barplot_counts.pdf",
    width = 7,
    height =5)

##recall that total counts has been stored as sample_sum
barplot(sample_sum, 
        ylab = "count total", 
        xlab = "Samples",
        col = "lightblue",
        ylim = c(0, 1100))
abline(h = 0, col = "black")
dev.off()
```

    ## png 
    ##   2

## 18. Renaming Gene 2 to Gene 2A

``` r
##renaming gene 2 to Gene2A
rownames(data_filtered)[2] <- "Gene2A"
```

## 19. Removing the lowest-mean gene from count

``` r
##Identifying which gene has the lowest mean
which.min(gene_mean)
```

    ## ALAS2 
    ##     5

``` r
##removing the counts for the gene with the lowest mean
data_filtered <- data_filtered[-5, ]
 

data_filtered
```

    ##        Sample_1 Sample_2 Sample_3 Sample_4 Sample_5
    ## STAT1       100      230      320      270      110
    ## Gene2A      180      110       65       85      190
    ## ISG15        50      110      195      140       55
    ## GBP5         80      170      280      220       75
    ## RPL13A      140      136      143      139      135

## 20. Conclusion

This mock study aimed at investigating the differential gene expression
between Healthy and Tb individuals using a small panel of immune-related
genes. The dataset included expression profile six genes: “STAT1”
“ISG15” “GBP5” “ALAS2” “RPL13A” in 5 samples (3 of which are TB samples
while the other 2 were healthy samples).

Several immune-related genes including STAT1, ISG15, and GBP5 showed
clear higher average (mean) in TB samples than in Healthy controls,
consistent with an activated immune response in TB. STAT1 and GBP5 were
shown to have the largest expression variability across samples and the
greatest difference between conditions. While STAT1 has the highest
overall mean (206.0) and the largest dynamic range (220), GBP5 however
had the highest fold-change proxy (≈2.88) between TB and Healthy. Also,
IFIT3 was shown to have the lowest fold-change proxy (≈0.47) with TB
samples having signifcant low expression counts compared with Healthy
samples.

Considering expression count, Sample_3 has the highest library size
(1068) and together with Sample_4 was labelled as HighLoad samples
(total \> median). It is important to note that HighLoad were only
observed among TB cases in this dataset, which may reflect a strong
correlation between the expression of these genes with TB infection.
