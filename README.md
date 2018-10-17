# eBWTclust

### Overview

The *eBWTclust* suite can be used to discover SNPs/indels between two sets of reads (fasta) *without* aligning them to  a reference genome (alignment-free, reference-free) with just a scan of the extended Burrows Wheeler Transform (e-BWT) of the sets of reads and the LCP and gSA arrays. The output is a fasta file (in KisSNP2 format) where sequences are the contexts surrounding the identified SNPs.  This suite is composed by the following modules: 

- **eBWTclust** partitions the eBWT of a set of reads in clusters corresponding to the same nucleotide in the reference genome. Output: a ".clusters" file.
- **clust2snp** analyzes the clusters produced by eBWTclust and detects SNPs and indels. Output: a ".snp" file (this is actually a fasta file containing pairs of reads testifying the variations).
- **snp2fastq** converts the ".snp" file produced by clust2snp into a ".fastq" file (with fake base qualities) ready to be aligned.
- **sam2vcf** converts the ".sam" file produced by aligning the above ".fastq" into a ".vcf" file containing the variations. 
- **compareVCF** compares the above ".vcf" file against a ground truth ".vcf" file, and computes sensitivity/specificity of the eBWTclust procedure. 

Note: *eBWTclust* and *clust2snp* require the Enhanced Generalized Suffix Array (EGSA) of the sets of reads (https://github.com/felipelouza/egsa or https://github.com/giovannarosone/BCR_LCP_GSA) to be built beforehand. 

The suite finds its main use in applications where no reference genome is known (alignment-free, reference-free variation discovery). For this, use the following pipeline:

**[EGSA|BCR] -> eBWTclust -> clust2snp** 

Note that the above pipeline finds many SNPs/indels twice: one time on the forward strand and one on the reverse complement strand. If a reference is available, one can extend the above pipeline to produce a vcf file. For this, one can use the extended pipeline 

**[EGSA|BCR] -> eBWTclust -> clust2snp -> snp2fastq -> bwa-mem -> sam2vcf**

The latter pipeline keeps only one copy of SNPs/indels that have been detected twice. To conclude, one can validate the generated VCF against a ground-truth VCF generated usng a standard pipeline (for example, BWA + samtools + bcftools) by using the tool **compareVCF**. The script **pipeline.sh** automates the whole pipeline.

The paper describing the theory behind the tool (eBWT positional clustering) has been published in:

---

*Nicola Prezza, Nadia Pisanti, Marinella Sciortino and Giovanna Rosone: Detecting mutations by eBWT. WABI 2018. Leibniz International Proceedings in Informatics, LIPIcs , 2018, 113, art. no. 3, Schloss Dagstuhl--Leibniz-Zentrum für Informatik.*

---

A pre-print version can be found here: https://arxiv.org/abs/1805.01876. 


### Install

~~~~
#download eBWTclust, EGSA, and BCR
git clone https://github.com/nicolaprezza/eBWTclust
git clone https://github.com/felipelouza/egsa
git clone https://github.com/giovannarosone/BCR\_LCP\_GSA

#build eBWTclust
cd eBWTclust
mkdir build
cd build
cmake ..
make

#build egsa
cd egsa
make compile BWT=1
~~~~

### Run

Enter the folder with the two fasta files _reads1.fasta_  and _reads2.fasta_ (i.e. the reads of the two samples). We assume that executables 'egsa', 'eBWTclust', and 'clust2snp' are global. 

~~~~
#Step 1: optional, but considerably increases sensitivity of the tool. Insert in reads1.fasta also the reverse-complement of the reads. Repeat with reads2.fasta

#Count reads in the first sample and concatenate fasta files
nreads1=`grep ">" reads1.fasta | wc -l`
cat reads1.fasta reads2.fasta > ALL.fasta

#Build the EGSA of the sets of reads
egsa ALL.fasta 0

#Build cluster file (do this in the same folder containing all other files)
eBWTclust -i ALL.fasta

#Call SNPs (do this in the same folder containing all other files)
clust2snp -i ALL.fasta -n ${nreads1}

#File ALL.snp.fasta now contains identified SNPs/indels

~~~~
