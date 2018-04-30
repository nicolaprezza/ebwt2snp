# eBWTclust

### Overview

This suite is composed (for now) by two main modules: **eBWTclust** partitions the eBWT of a set of reads in clusters corresponding to the same nucleotide in the reference genome (paper describing this theory coming soon!), while **clust2snp** analyzes the clusters produced by eBWTclust and detects SNPs.

Together, these two tools can be used to discover SNPs between two sets of reads (fasta) *without* aligning them to  a reference genome (alignment-free) with just a scan of the extended Burrows Wheeler Transform (e-BWT) of the sets of reads. The output is a fasta file (in KisSNP2 format) where sequences are the contexts surrounding the identified SNPs. 

 **eBWTclust** and **clust2snp** require the Enhanced Generalized Suffix Array (EGSA) of the sets of reads (https://github.com/felipelouza/egsa) to be built beforehand. 

### Install

~~~~
#download bw-snp and EGSA
git clone https://github.com/nicolaprezza/eBWTclust
git clone https://github.com/felipelouza/egsa

#build bw-snp
cd eBWTclust
mkdir build
cd build
cmake ..
make

#build egsa
cd egsa
egsa_dir=`pwd`
make compile BWT=1
~~~~

### Run

Enter the folder with the two fasta files _reads1.fasta_  and _reads2.fasta_. It is assumed that the variable _egsa_dir_ stores the installation path of _EGSA_ (see installation above).

~~~~
#Step 1: optional, but considerably increases sensitivity of the tool. Insert in reads1.fasta also the reverse-complement of the reads. Repeat with reads2.fasta

#Count reads in first fasta and concatenate fasta files
nreads1=`grep ">" reads1.fasta | wc -l`
cat reads1.fasta reads2.fasta > ALL.fasta

#Build the EGSA of the sets of reads
make -C ${egsa_dir} run DIR=`pwd` INPUT=ALL.fasta K=0 CHECK=0

#Build cluster file (do this in the same folder containing all other files)
eBWTclust -i ALL.fasta

#Call SNPs (do this in the same folder containing all other files)
clust2snp -i ALL.fasta -n ${nreads1}

#File ALL.snp.fasta now contains identified SNPs
~~~~
