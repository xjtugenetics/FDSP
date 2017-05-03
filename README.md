# Functional disease-associated SNPs prediction (FDSP)

# Introduction

We developed functional disease-associated SNPs prediction (FDSP), a toolkit for predicting disease risk SNPs. Our tools includes two parts: pre-processing pipeline and FDSP package in R.

Several scripts consist of pre-processing pipeline. These scripts perform calculating linkage disequilibrium(LD) and minor allele frequency(MAF), creating positive dataset and negative dataset according to lead SNPs of diseases, and annotating SNPs with epigenetics elements.

FDSP package includes three sections: training set and test set creation, training prediction model with machine learning algorithm, and disease-associated SNPs prediction.

FDSP is open source and can be run on different platforms.

# Author

FDSP is written by Chen Yi-xiao, Yao shi and Dong shan-shan from Genetics and Bioinformatics Lab, Xi An Jiao Tong Univ, Sch Life Sci & Technol, Minist Educ, Key Lab Biomed Informat Engn, Xian 710049, Peoples R China.


# License

This software is distributed under the terms of GPL 3.0

# Source

https://github.com/xjtugenetics/SNPprediction

# Contact

You can contact yangtielin@mail.xjtu.edu.cn when you have any questions, suggestions, comments, etc. Please describe with more details, and attach your command line and log messages if possible.

# Requirements

awk, R (>= 3.3.2), PLINK, bedtools

R packages: caret, methods, randomForest, doParallel, parallel

# Preparation

## Setup folder framework

In /path/prep, there are three for calculating preparation. You MUST set them before starting calculation.

```
cd /path/prep
vi files.txt
ls
annotation.txt  files.txt  parameters.txt
```
There are three files in this directory. These files include some files or options involving in scripts. You can use this procedure with defualt setting or editing and updating by youself.

## Locate whole-genome SNPs files with position and MAF

We provide series of files contating whole-genome SNPs information, including position, names and MAF provided by 1000 Genomes Project building on hg19(Human Genome version 19). You can also use files provided by youself. If you prepare files by youself, leave the format as example files and put the files into /path/extdata.

```
(chr,position,rs ID, MAF)
1,15444772,rs12744957,0.3728
```

# Acquisition of labeled SNPs

The labeled SNPs were consisted of positive and negative SNPs. For the labeled positive SNPs, we firstly obtained index SNPs from the public SNP-trait association databases: GWAS Catalog (https://www.ebi.ac.uk/gwas/) with the threshold of P value < 5 ×10-8. We secondly obtained SNPs that were in strong LD (r2 ≥ 0.8) with each index SNP using the 1000 Genomes phase III data. 

For the negative SNPs, we created four sets of SNPs with different distances to positive SNPs. SNPs with similar MAF to the positive SNP were remained (MAF difference < 0.01). The full set using 20 negatives per positive. All negative SNPs were filtered to remove overlap with positive SNPs.

Before starting labeled SNPs selection, please edit preparation files first.

## Edit /path/prep/files.txt

This file contains the list of files paths necessaried in analyses. 

Leave the format as is

```
File_pur	File_path
lead_SNPs	../exp/lead_SNPs.txt 
bfile	NA
bfile_list	../exp/ldlist.txt
```

Lead_SNPs is location of a file including lead SNPs associated with complex diseases found by previous GWASs. 

We suggest you use PLINK genotype files of 1000 Genomes Project for calculating LD. Because the data could be very large so that computation speed may quitely slow, you can devide dataset into several files according to the positions of SNPs on chromosomes. You can locate a PLINK file or a directory including several PLINK files foe calculating LD.

If you have PLINK .bed/bim/fam files including all samples and all chromosomes information, please add the files name on bfile line (without postfix).

If you have multiple PLINK files devided by chromosomes, please create a file containing a list of PLINK .bed/bim/fam files. Then add the directory on bfile_list.

## Calculate LD

```
cd /path/pre-prosessing/scripts
sh calculate_LD.sh
```
PLEASE WAIT until this job finishes. LD calculating results store in data folder.

## Select positive set

```
cd /path/pre-prosessing/scripts
sh select_positive.sh
```
You now have the positive set in bed format store in data folder.

## Edit /path/prep/parameters.txt

This file contains the list of parameters paths necessaried in analyses. 

Leave the format as is

```
Para	Value
r2_ld	0.1
po_threshold	0.98
na_distant	40000,200000,1000000,5000000
```
All values should be updated by you.

R2_ld is the threshold of calculating LD at which SNPs are excluded from selecting negative set. I.e. if a SNP LD r-square with SNPs in the positive set is greater than eg 0.1, it is excluded from negative set.

Po_threshold limits the LD r-square value of SNPs included in positive set. SNPs with greater than eg 0.8 value of LD r-square with lead SNPs consist of positive set with lead SNPs. It is set as 0.8 in general.

Na_distant is the distance between SNPs in positive set and SNPs in selected negative set. I.e all SNPs will be devided into four sets SNPs with different distances to positive SNPs eg the maximum distance in each group was 40 kb, 200 kb, 1000 kb, and 5000 kb.

## Select negative set

```
cd /path/pre-prosessing/scripts
sh select_negative.sh
```
PLEASE WAIT until this job finishes. It may takes several hours. You will get the negative set in bed format.

# Feature annotation

In this step, we provided a script of functional annotation. It can annotate all SNPs based on the epigenomic data (or any type of data in bed format).

## Edit /path/prep/annotation.txt

This file consist of a set of genomic elements which be used as features in establishing model. The genomic features could be genomic coordinates in bed format. For example, we have add the ChIP-seq peak data in bed format for histone marks.

```
Element	File
E003-H3K23ac	/path/E003-H3K23ac.bed
E001-H3K27me3	/path/E001-H3K27me3.bed
```

## Annotate SNPs with features

```
cd /path/pre-prosessing/scripts
sh bedtools_annotation.sh
```
PLEASE WAIT until this job finishes.

Now you get the results of annotation steps. You can train models with these results for predicting SNPs.

# Model generation, evaluation and optimization

## Running FDSP package in R

As an example, we are going to predict candidate risk SNPs of diabetes. We have annotated some SNPs with epigenomic elements for training model. We have also provided a example file including some function unknown SNPs to test the model.

Make sure you have installed FDSP package before running below commands.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
library('FDSP')
SNPanno<-read.csv(file.path(system.file('extdata', 'SNP_ori.csv', 
			package="FDSP")), header=T)
save(SNPanno, file='example.SNPanno.Rda')
```
Make sure the first column of input files is "SNP", as well as the last column is "status". "SNP" includes the name of known SNPs. "Status" 
means the status of SNPs. Tag "1" means risk SNPs (positive set) and "0" means non-risk SNPs (negative set).

## Filter features

Filtering high correlation features. Features lower than threshold remains for next step.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
load(file.path(system.file('extdata', 'example.SNPanno.Rda', package='FDSP')))
 SNPdatafilter <- filter_features(SNPanno)
save(SNPdatafilter, file="example.SNPdatafilter.Rda")
```

## Create dataset

Create train dataset and dataset after filtering high correlation features.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
load(file.path(system.file('extdata', 'example.SNPdatafilter.Rda', package='FDSP')))
dataset<-create_dataset(SNPdatafilter,numbercv=5)
train_data<-dataset[[1]]
test_data<-dataset[[2]]
save(dataset, file="example.dataset.Rda")
```

## Train model

Train a model with training set. Choose type of model and number of cross-validation you want. This package supports five types of model: CSimca(CSimca), KNN(knn), radial basis function kernel(svmRadial), C5.0(C5.0) and random forest(rf). You can choose type of model you want with changing "method" option(shown in brackets).

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
load(file.path(system.file('extdata', 'example.dataset.Rda', package='FDSP')))
load(file.path(system.file('extdata', 'example.SNPdatafilter.Rda', package='FDSP')))
train_data<-dataset[[1]]
model <- model_train(SNPdatafilter, train_data, method="rf", numbercv = 5)
save(model, file="example.model.Rda")
```

## Model evaluation

Get prediction results, confusion matrix, F1 score and feature importance.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
load(file.path(system.file('extdata', 'example.model.Rda', package='FDSP')))
load(file.path(system.file('extdata', 'example.SNPdatafilter.Rda', package='FDSP')))
test_data<-dataset[[2]]
evaluate_data <- model_evaluate(model, test_data)
prediction_results <- evaluate_data[[1]]
confusion_matrix <- evaluate_data[[2]]
F1_score <- evaluate_data[[3]]
feature_importance <- evaluate_data[[4]]
save(evaluate_data,  file='example.evaluate_data.Rda')
```

## Select model

Select a best-performed model with appropriate number of features.

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
load(file.path(system.file('extdata', 'example.SNPdatafilter.Rda', package='FDSP')))
load(file.path(system.file('extdata', 'example.evaluate_data.Rda', package='FDSP')))
model_best_set<-select_features(SNPdatafilter,evaluate_data,from=10,to=100,sep=10)
model_best <- model_best[[1]]
feature_best <- model_best[[4]]
save( model_best_set, file='example.model_best_set.Rda')
```

# Predict SNPs

Predict new candidate SNPs accroding to their epigenomic elements and trained model.
Make sure the first column of customer files is "SNP".

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
test <- read.csv(file.path(system.file('extdata','test.csv',header=T)))
rownames(test) <- test[,1]
feature_best <- test[, feature_best ]
predict_test <- SNP_predict(model_best,feature_best)
predict_test
```

# Save results

```{r warning=FALSE, message=FALSE, tidy=TRUE, eval=FALSE}
save(predict_test, file='example.predict_test.Rda')
write.table(as.data.frame( predict_test), file='example.results.txt', row=T, col=F, quote=F, sep="\t")
```

