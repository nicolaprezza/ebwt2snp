# Usage: pipeline.sh <reads1.fastq> <reads2.fastq> <reference.fasta> <c>
#
# behaviour (the steps are skipped if the file they will produce already exists)
# 1.  Converts input reads to fasta -> reads1.fasta reads2.fasta 
# 2.  Adds the reverse complements to the reads and concatenates the two read files -> reads1.reads2.frc.fasta
# 3.  Builds EGSA -> reads1.reads2.frc.fasta.gesa
# 4.  Run eBWTclust -> reads1.reads2.frc.fasta.clusters 
# 5.  Run clust2snp with parameter -c <c> -> reads1.reads2.frc.<c>.snp
# 6.  Run snp2fastq -> reads1.reads2.frc.<c>.snp.fastq
# 7.  Builds BWA MEM index of reference.fasta -> reference.fasta.{amb,ann,bwt,fai,pac,sa} files
# 8.  Create reference of reads1.fasta using BWA MEM + bcftools + vcfconsensus -> reads1.reference.fasta
# 9.  Builds BWA MEM index of reads1.reference.fasta -> reads1.reference.fasta.{amb,ann,bwt,fai,pac,sa} files
# 10.  Aligns reads1.reads2.frc.<c>.snp.fastq on reads1.reference.fasta -> reads1.reads2.frc.<c>.snp.sam
# 11. Generates VCF (eBWTclust calls) using sam2vcf -> reads1.reads2.frc.<c>.snp.sam.vcf
# 12. Generates VCF (bcftools calls) using BWA MEM + bcftools -> reads1.reads2.bcftools.vcf
# 13. Generate report containing running times of eBWTclust pipeline, BWA+bcftools pipeline, and precision/recall of eBWTclust pipeline (using BWA+bcftools pipeline as ground truth)

# Requires the following executables to be globally visible (in addition to the executables of eBWTclust): 
# - fastq2fasta.sh (https://github.com/nicolaprezza/bioinfo-tools)
# - seqtk
# - egsa (https://github.com/felipelouza/egsa)
# - gsufsort (https://github.com/felipelouza/gsufsort)
# - bwa
# - samtools
# - bcftools
# - bgzip
# - vcftools

#Byte size of EGSA fields
GSAtext=4
GSAsuff=1
LCP=1

WD=$(dirname $(readlink -f "$1"))

#extract full path of working directory without trailing /
echo "Working directory: "${WD}

#extract file name and extension of reads1
READS1=`basename $1`
EXT1=${READS1##*.}
READS1="${READS1%.*}"

echo "input 1: "${READS1}"."${EXT1}

#extract file name and extension of reads2
READS2=`basename $2`
EXT2=${READS2##*.}
READS2="${READS2%.*}"

#extract file name of reference
REF=`basename $3`

#parameter -c in clust2snp (max reads in consensus. Ideally big, but bigger slows down computation)
C=$4

#parameter -m in clust2snp (we keep this fixed: at least 5 reads per individual in cluster)
M=5

TIME_EGSA=${WD}/egsa.time
TIME_EBWTCLUST=${WD}/eBWTclust.time
TIME_CLUST2SNP=${WD}/clust2snp_${C}.time
TIME_BWAMEM=${WD}/bwamem.time
TIME_BWAIDX=${WD}/bwaindex.time
TIME_BCFTOOLS=${WD}/bcftools.time
REPORT_OUT=${WD}/report_${C}

samtool_version=`samtools 2>&1 >/dev/null | grep Version | cut -d' ' -f 2 | cut -d. -f 1`

echo "input 2: "${READS2}"."${EXT2}

# 1.  Converts input reads to fasta, if not already in fasta -> reads1.fasta reads2.fasta 

if [ ! -f ${WD}/${READS1}.fasta ]; then
	echo ${WD}/${READS1}.fasta" not found. Converting fastq to fasta ..."
	fastq2fasta.sh ${WD}/${READS1}.fastq > ${WD}/${READS1}.fasta
fi

if [ ! -f ${WD}/${READS2}.fasta ]; then
	echo ${WD}/${READS2}.fasta" not found. Converting fastq to fasta ..."
	fastq2fasta.sh ${WD}/${READS2}.fastq > ${WD}/${READS2}.fasta
fi

#number of reads in individual 1
N=`cat ${WD}/${READS1}.fasta | grep '^>' | wc -l`
N=$((N*2))

# 2.  Adds the reverse complements to the reads and concatenates the two read files -> reads1.reads2.frc.fasta

if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta ]; then
	echo "Adding reverse complement to the reads and building main fasta "${WD}/${READS1}.${READS2}.frc.fasta" ..."
	seqtk seq -r ${WD}/${READS1}.fasta > ${WD}/${READS1}.rc.fasta
	seqtk seq -r ${WD}/${READS2}.fasta > ${WD}/${READS2}.rc.fasta
	cat ${WD}/${READS1}.fasta ${WD}/${READS1}.rc.fasta ${WD}/${READS2}.fasta ${WD}/${READS2}.rc.fasta > ${WD}/${READS1}.${READS2}.frc.fasta
	rm ${WD}/${READS1}.rc.fasta ${WD}/${READS2}.rc.fasta
fi

# 3.  If EGSA/BCR files do not exist (both), then builds EGSA -> reads1.reads2.frc.fasta.gesa
if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.gesa ]; then
	if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.out ]; then
		echo "building EGSA ..."
		#/usr/bin/time -v egsa -vvv -c ${WD}/${READS1}.${READS2}.frc.fasta 0 > ${TIME_EGSA} 2>&1
		/usr/bin/time -v gsufsort ${WD}/${READS1}.${READS2}.frc.fasta --gesa ${GSAtext} ${GSAsuff} ${LCP} > ${TIME_EGSA} 2>&1
		mv ${WD}/${READS1}.${READS2}.frc.fasta.${GSAtext}.${GSAsuff}.${LCP}.1.gesa ${WD}/${READS1}.${READS2}.frc.fasta.gesa
		
		rm -rf ${WD}/tmp
		rm -rf ${WD}/partition
	fi
fi

# 4.  Run eBWTclust -> reads1.reads2.frc.fasta.clusters 

if [ ! -f ${WD}/${READS1}.${READS2}.frc.fasta.clusters ]; then
	echo "running eBWTclust ..."
	/usr/bin/time -v eBWTclust -i ${WD}/${READS1}.${READS2}.frc.fasta -x ${LCP} -y ${GSAtext} -z ${GSAsuff} > ${TIME_EBWTCLUST} 2>&1
fi

# 5.  Run clust2snp with parameters m, c -> reads1.reads2.frc.<c>.snp

if [ ! -f ${WD}/${READS1}.${READS2}.frc.${C}.snp ]; then
	echo "running clust2snp ..."
	/usr/bin/time -v clust2snp -i ${WD}/${READS1}.${READS2}.frc.fasta -n $N -m $M -c $C -x ${LCP} -y ${GSAtext} -z ${GSAsuff} > ${TIME_CLUST2SNP} 2>&1
	mv ${WD}/${READS1}.${READS2}.frc.snp ${WD}/${READS1}.${READS2}.frc.${C}.snp
fi

# 6.  Run snp2fastq -> reads1.reads2.frc.<c>.snp.fastq

if [ ! -f ${WD}/${READS1}.${READS2}.frc.${C}.snp.fastq ]; then
	echo "converting .snp to .fastq ..."
	snp2fastq ${WD}/${READS1}.${READS2}.frc.${C}.snp
fi


# 7.  Builds BWA MEM index of reference.fasta -> reference.fasta.{amb,ann,bwt,fai,pac,sa} files

if [ ! -f ${WD}/${REF}.amb ]; then
	echo "Indexing "${WD}/${REF}" (BWA) ..."
	bwa index ${WD}/${REF}
fi

# 8.  Create reference of reads1.fasta using BWA MEM + bcftools + vcfconsensus

if [ "$samtool_version" -eq "1" ]; then
	if [ ! -f ${WD}/${READS1}.reference.fasta ]; then
		echo "Computing reference of "${READS1}". Storing it to file "${WD}/${READS1}.reference.fasta" ..."
		bwa mem ${WD}/${REF} ${WD}/${READS1}.fastq > ${WD}/alignment.tmp.sam
		samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam
		samtools sort ${WD}/alignment.tmp.bam > ${WD}/alignment.tmp.sorted.bam
		samtools index ${WD}/alignment.tmp.sorted.bam
		bcftools mpileup -f ${WD}/${REF} ${WD}/alignment.tmp.sorted.bam | bcftools call -mv -o ${WD}/calls.tmp.vcf
		bgzip ${WD}/calls.tmp.vcf
		vcf2fasta.sh ${WD}/calls.tmp.vcf.gz ${WD}/${REF} > ${WD}/${READS1}.reference.fasta
		rm *.tmp*
	fi
else
	if [ ! -f ${WD}/${READS1}.reference.fasta ]; then
		echo "Computing reference of "${READS1}". Storing it to file "${WD}/${READS1}.reference.fasta" ..."
		bwa mem ${WD}/${REF} ${WD}/${READS1}.fastq > ${WD}/alignment.tmp.sam
		samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam
		samtools sort ${WD}/alignment.tmp.bam ${WD}/alignment.tmp.sorted
		samtools index ${WD}/alignment.tmp.sorted.bam

		#OLD SAMTOOLS/BCFTOOLS
		samtools mpileup -uD -f ${WD}/${REF} ${WD}/alignment.tmp.sorted.bam > ${WD}/mpileup.tmp
		bcftools view -vc ${WD}/mpileup.tmp > ${WD}/calls.tmp.vcf
		
		bgzip ${WD}/calls.tmp.vcf
		vcf2fasta.sh ${WD}/calls.tmp.vcf.gz ${WD}/${REF} > ${WD}/${READS1}.reference.fasta
		rm ${WD}/*.tmp*
	fi
fi

# 9.  Builds BWA MEM index of reads1.reference.fasta -> reads1.reference.fasta.{amb,ann,bwt,fai,pac,sa} files

if [ ! -f ${WD}/${READS1}.reference.fasta.amb ]; then
	echo "Indexing "${WD}/${READS1}.reference.fasta" (BWA) ..."
	/usr/bin/time -v bwa index ${WD}/${READS1}.reference.fasta > ${TIME_BWAIDX} 2>&1
fi

# 10.  Aligns reads1.reads2.frc.<c>.snp.fastq on reads1.reference.fasta -> reads1.reads2.frc.<c>.snp.sam

if [ ! -f ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam ]; then
	echo "Aligning "${WD}/${READS1}.${READS2}.frc.${C}.snp.fastq" on "${WD}/${READS1}.reference.fasta" ..."
	bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS1}.${READS2}.frc.${C}.snp.fastq > ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam 
fi

# 11. Generates VCF (eBWTclust calls) using sam2vcf -> reads1.reads2.frc.<c>.snp.sam.vcf

if [ ! -f ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam.vcf ]; then
	echo "Generating eBWTclust's VCF in "${WD}/${READS1}.${READS2}.frc.${C}.snp.vcf" ..."
	sam2vcf -s ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam
fi

# 12. Generates VCF (bcftools calls) using BWA MEM + bcftools -> reads1.reads2.bcftools.vcf

if [ "$samtool_version" -eq "1" ]; then
	if [ ! -f ${WD}/${READS1}.${READS2}.bcftools.vcf ]; then
		echo "Calling SNPs using BWA + bcftools. Storing SNPs to file "${WD}/${READS1}.${READS2}.bcftools.vcf" ..."
		/usr/bin/time -v bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS2}.fastq > ${WD}/alignment.tmp.sam 2> ${TIME_BWAMEM}
		/usr/bin/time -v samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools sort ${WD}/alignment.tmp.bam > ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools index ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v bcftools mpileup -f ${WD}/${READS1}.reference.fasta ${WD}/alignment.tmp.sorted.bam | bcftools call -mv -o ${WD}/${READS1}.${READS2}.bcftools.vcf 2>> ${TIME_BCFTOOLS}
		rm *.tmp*
	fi
else
	if [ ! -f ${WD}/${READS1}.${READS2}.bcftools.vcf ]; then
		echo "Calling SNPs using BWA + bcftools. Storing SNPs to file "${WD}/${READS1}.${READS2}.bcftools.vcf" ..."
		/usr/bin/time -v bwa mem ${WD}/${READS1}.reference.fasta ${WD}/${READS2}.fastq > ${WD}/alignment.tmp.sam 2> ${TIME_BWAMEM}
		/usr/bin/time -v samtools view -b -S ${WD}/alignment.tmp.sam > ${WD}/alignment.tmp.bam 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools sort ${WD}/alignment.tmp.bam ${WD}/alignment.tmp.sorted 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v samtools index ${WD}/alignment.tmp.sorted.bam 2>> ${TIME_BCFTOOLS}

		#OLD SAMTOOLS/BCFTOOLS
		/usr/bin/time -v samtools mpileup -uD -f ${WD}/${READS1}.reference.fasta ${WD}/alignment.tmp.sorted.bam > ${WD}/mpileup.tmp 2>> ${TIME_BCFTOOLS}
		/usr/bin/time -v bcftools view -vc ${WD}/mpileup.tmp > ${WD}/${READS1}.${READS2}.bcftools.vcf

		rm ${WD}/*.tmp*
	fi
fi


# 13. Generate report containing precision/recall of eBWTclust pipeline (using BWA+bcftools pipeline as ground truth)

if [ ! -f ${WD}/${READS1}.${READS2}.report.${C}.tsv ]; then
	compareVCF -1 ${WD}/${READS1}.${READS2}.frc.${C}.snp.sam.vcf -2 ${WD}/${READS1}.${READS2}.bcftools.vcf > ${REPORT_OUT}
fi






