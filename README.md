# eBWTclust

### Overview

This suite is composed (for now) by two main modules: **eBWTclust** partitions the eBWT of a set of reads in clusters corresponding to the same nucleotide in the reference genome, while **clust2snp** analyzes the clusters produced by eBWTclust and detects SNPs (only haploid samples supported for now).

The paper describing the theory behind the tool (eBWT positional clustering) can be found here: https://arxiv.org/abs/1805.01876

Together, these two tools can be used to discover SNPs between two sets of reads (fasta) *without* aligning them to  a reference genome (alignment-free, reference-free) with just a scan of the extended Burrows Wheeler Transform (e-BWT) of the sets of reads and the LCP and gSA arrays. The output is a fasta file (in KisSNP2 format) where sequences are the contexts surrounding the identified SNPs. Note that most of the SNPs will be found (and reported) twice: one time for the forward strand, and one for the reverse strand.

 **eBWTclust** and **clust2snp** require the Enhanced Generalized Suffix Array (EGSA) of the sets of reads (https://github.com/felipelouza/egsa) to be built beforehand. 
 
 The additional tool **validate** can be used to compute precision and recall (sensitivity) of the output (in KisSNP2 format), given the fasta reference file of the first individual and the .vcf file describing the ground truth (i.e. the real variants between the two samples). All details about the validation step can be found in the paper. 

### Install

~~~~
#download eBWTclust and EGSA
git clone https://github.com/nicolaprezza/eBWTclust
git clone https://github.com/felipelouza/egsa

#build eBWTclust
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

Enter the folder with the two fasta files _reads1.fasta_  and _reads2.fasta_ (i.e. the reads of the two samples). It is assumed that the variable _egsa_dir_ stores the installation path of _EGSA_ (see installation above).

~~~~
#Step 1: optional, but considerably increases sensitivity of the tool. Insert in reads1.fasta also the reverse-complement of the reads. Repeat with reads2.fasta

#Count reads in the first sample and concatenate fasta files
nreads1=`grep ">" reads1.fasta | wc -l`
cat reads1.fasta reads2.fasta > ALL.fasta

#Build the EGSA of the sets of reads
make -C ${egsa_dir} run DIR=`pwd` INPUT=ALL.fasta K=0 CHECK=0

#Build cluster file (do this in the same folder containing all other files)
eBWTclust -i ALL.fasta

#Call SNPs (do this in the same folder containing all other files)
clust2snp -i ALL.fasta -n ${nreads1}

#File ALL.snp.fasta now contains identified SNPs

#Compute sensitivity and precision of the output. Be aware that the argument of -f must match the vcf specified with -v!
validate -v true_variants_between_two_samples.vcf -c ALL.snp.fasta -f ref_of_first_sample.fasta
~~~~
