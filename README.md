# Functional disease-associated SNPs prediction (FDSP)

# Introduction

Genome-wide association studies (GWASs) is an effective strategy to identify susceptibility loci for human complex diseases. However, missing heritability is still a big problem. Here we developed a pipeline named functional disease-associated SNPs prediction (FDSP), to identify novel susceptibility loci for complex diseases based on the interpretation of the functional features for known disease-associated variants with machine learning.

# License

This software is distributed under the terms of GPL 2.0

# Source

https://github.com/xjtugenetics/FDSP

# Contact

You can contact yangtielin@mail.xjtu.edu.cn when you have any questions, suggestions, comments, etc. Please describe in details, and attach your command line and log messages if possible.

# Requirements

awk, R (>= 3.3.2), PLINK, bedtools

R packages: caret, methods, randomForest, doMC, parallel, rrcovHD, C50, kernlab, plyr, ROCR, mlbench, RSNNS

# Preparation

### Setup folder framework

There are three files in the /path/prep directory. Users MUST set them or use the default files before starting calculation.

```
cd /path/prep
vi files.txt
ls
annotation.txt  files.txt  parameters.txt
```

### Locate whole-genome SNPs files with position and MAF

We provide series of files contating whole-genome SNPs information, including position, names and MAF provided by 1000 Genomes Project building on hg19(Human Genome version 19). You can also use files provided by youself with the following format:

```
(chr,position,rs ID, MAF)
1,15444772,rs12744957,0.3728
```

# Acquisition of labeled SNPs

The labeled SNPs were consisted of positive and negative SNPs. For the labeled positive SNPs, we firstly obtained index SNPs from the public SNP-trait association databases: GWAS Catalog (https://www.ebi.ac.uk/gwas/) with the threshold of P value < 5 ×10-8. We secondly obtained SNPs that were in strong LD (r2 ≥ 0.8) with each index SNP using the 1000 Genomes phase III data. 

For the negative SNPs, we created four sets of SNPs with different distances to positive SNPs. SNPs with similar MAF to the positive SNP were remained (MAF difference < 0.01). The full set using 20 negatives per positive. All negative SNPs were filtered to remove overlap with positive SNPs.

Before starting labeled SNPs selection, please edit preparation files first.

### Edit /path/prep/files.txt

This file contains the list of files paths necessary in analyses. 

Leave the format as is 

```
File_pur	File_path
lead_SNPs	../exp/lead_SNPs.txt 
bfile	NA
bfile_list	../exp/ldlist.txt
```

Lead_SNPs.txt is location of a file including index SNPs associated with complex diseases found by previous GWASs. 

We suggest you use PLINK genotype files of 1000 Genomes Project for calculating LD. Because the data could be very large so that computation speed may quitely slow, you can devide dataset into several files according to the positions of SNPs on chromosomes. You can locate a PLINK file or a directory including several PLINK files foe calculating LD.

If you have PLINK .bed/bim/fam files including all samples and all chromosomes information, please add the files name on bfile line (without postfix).

If you have multiple PLINK files devided by chromosomes, please create a file containing a list of PLINK .bed/bim/fam files. Then add the directory on bfile_list.

### Calculate LD

```
cd /path/pre-prosessing/scripts
sh calculate_LD.sh
```
PLEASE WAIT until this job finishes. LD calculating results store in data folder.

### Select positive set

```
cd /path/pre-prosessing/scripts
sh select_positive.sh
```
You now have the positive set in bed format store in data folder.

### Edit /path/prep/parameters.txt

This file contains the list of parameters paths necessaried in analyses. 

Leave the format as is

```
Para	Value
r2_ld	0.1
po_threshold	0.8
na_distant	40000,200000,1000000,5000000
```
Users can use the default settings or update the values.

R2_ld is the threshold of calculating LD at which SNPs are excluded from selecting negative set. i.e. if a SNP LD r-square with SNPs in the positive set is greater than eg 0.1, it is excluded from negative set.

Po_threshold limits the LD r-square value of SNPs included in positive set. SNPs with greater than eg 0.8 value of LD r-square with lead SNPs consist of positive set with lead SNPs. It is set as 0.8 in general.

Na_distant is the distance between SNPs in positive set and SNPs in selected negative set. i.e all SNPs will be devided into four sets SNPs with different distances to positive SNPs eg the maximum distance in each group was 40 kb, 200 kb, 1000 kb, and 5000 kb.

### Select negative set

```
cd /path/pre-prosessing/scripts
sh select_negative.sh
```
PLEASE WAIT until this job finishes. It may takes several hours. You will get the negative set in bed format.

# Feature annotation

In this step, we provided a script of functional annotation. It can annotate all SNPs based on the epigenomic data (or any type of data in bed format).

### Edit /path/prep/annotation.txt

This file consist of a set of genomic elements which be used as features in establishing model. The genomic features could be genomic coordinates in bed format. For example, we have add the ChIP-seq peak data in bed format for histone marks.

```
Element	File
E003-H3K23ac	/path/E003-H3K23ac.bed
E001-H3K27me3	/path/E001-H3K27me3.bed
```

### Annotate SNPs with features

```
cd /path/pre-prosessing/scripts
sh bedtools_annotation.sh
```
PLEASE WAIT until this job finishes.

Now you get the results of annotation steps. You can train models with these results for predicting SNPs.

# Model generation, evaluation and optimization

### Running MP package in R

As an example, we are going to predict candidate risk SNPs of diabetes. We have annotated some SNPs with epigenomic elements for training model. We have also provided a example file including some function unknown SNPs to test the model.

Make sure you have installed MP package before running below commands.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
library('MP')
SNPanno<-read.csv(file.path(system.file('extdata', 'SNP_ori.csv', 
			package="MP")), header=T)
			
rownames(SNPanno)<-SNPanno[,1]
SNPanno <- SNPanno[, -1]
save(SNPanno, file='example.SNPanno.Rda')
```
Make sure the first column of input files is "SNP", as well as the last column is "Class". "SNP" includes the name of known SNPs. "Class" 
means the status of SNPs. Tag "1" means risk SNPs (positive set) and "0" means non-risk SNPs (negative set).

### Filter features

Filtering high correlation features. Features lower than threshold remains for next step.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
data(SNPanno)
 SNPdatafilter <- filter_features(SNPanno)
save(SNPdatafilter, file="example.SNPdatafilter.Rda")
```

### Create dataset

Create train dataset and dataset after filtering high correlation features.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
data(SNPdatafilter)
dataset<-create_dataset(SNPdatafilter,numbercv=5)
no_cv <- 1
test_data <- dataset[[no_cv]]
train_data <- do.call(rbind, dataset[setdiff(1:length(dataset),no_cv)])
save(dataset, file="example.dataset.Rda")
```

### Train model

This step will select a best-performed model with appropriate number of features. You can choose type of model and number of cross-validation you want. This package supports to select the best performance model from four types of model: CSimca(CSimca), radial basis function kernel(svmRadial), C5.0(C5.0) and random forest(rf). 


```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
data(SNPdatafilter)
model <-  model_train(SNPdatafilter,method="cv", model="all", cores = 10,start=10, end=60, sep=10)
model_best <- model$model
feature_importance <- model$feature_importance
save(model, file="example.model.Rda")
```

You can just train a model you want with changing "method" option(shown in brackets). All avalible models show in https://topepo.github.io/caret/available-models.html

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
data(SNPdatafilter)
model_rf <-  model_train(SNPdatafilter,method="cv", model="rf", cores = 10,start=10, end=60, sep=10)
model_best <- model$model
feature_importance <- model$feature_importance
save( model_rf, file='example.model_rf.Rda')
```

### Predict SNPs

Predict new candidate SNPs accroding to their epigenomic elements and trained model.
Make sure the first column of customer files is "SNP".

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
data(model)
test <- read.csv(file.path(system.file('extdata','test.csv', package='MP')),header=T)
rownames(test) <- test[,1]
test <- test[,-1]
predict_test <- SNP_predict(model,test)
save(predict_test, file='example.predict_test.Rda')
```

# Save results

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
write.table(predict_test, file='example.results.txt', row=F, col=T, quote=F, sep="\t")
```

