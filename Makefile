# Annotate and visualize the Sitka spruce (Picea sitchensis) mitochondrion genome
# Written by Shaun Jackman @sjackman

# Target genome, Picea sitchensis
name=psitchensis

# Reference genome, Picea glauca
ref=LKAM00000000

# Number of threads
t=64

# gzip compression program
gzip=pigz -p$t

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

all:

install-deps:
	brew install pigz
	brew tap homebrew/science
	brew install bcftools bwa edirect fastqc samtools seqtk

.PHONY: all clean install-deps
.DELETE_ON_ERROR:
.SECONDARY:

# Fetch data from NCBI

$(name).fa $(ref).fa: %.fa:
	efetch -db nuccore -id $* -format fasta | seqtk seq | sed 's/^>/>$* /' >$@

# seqtk

# Interleave paired-end reads
%.fq.gz: %_1.fq.gz %_2.fq.gz
	seqtk mergepe $^ | $(gzip) >$@

# FASTQC

%.fastqc: %.fq.gz
	mkdir -p $@
	fastqc -t $t -o $@ $<

# BWA

# Index the target genome.
%.fa.bwt: %.fa
	bwa index $<

# Align sequences to the target genome.
$(ref).%.sam: %.fa $(ref).fa.bwt
	bwa mem -t$t -xintractg $(ref).fa $< >$@

# Align paired-end reads to the target genome.
$(ref).%.bam: %.fq.gz $(ref).fa.bwt
	bwa mem -t$t -p $(ref).fa $< | samtools view -h -F4 | samtools sort -@$t -o $@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.bam: %.sam
	samtools sort -@$t -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index $<

# Select properly paired reads.
%.proper.bam: %.bam
	samtools view -Obam -f2 -o $@ $<

# Calculate depth of coverage.
%.depth.tsv: %.bam
	(printf "Seq\tPos\tDepth\n"; samtools depth -a $<) >$@

# bcftools

# Call variants of reads aligned to a reference.
$(ref).%.vcf.gz: $(ref).%.bam $(ref).fa
	samtools mpileup -u -f $(ref).fa $< | bcftools call --threads=$t -c -v --ploidy=1 -Oz >$@

# Filter variants to select locations that differ from the reference.
%.filter.vcf.gz: %.vcf.gz
	bcftools filter -Oz -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && DP4[2] > 0 && DP4[3] > 0' $< >$@

# Index a VCF file.
%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

# Call the consensus FASTA sequence of a VCF file
$(ref).%.vcf.fa: $(ref).%.vcf.gz $(ref).%.vcf.gz.csi $(ref).fa
	bcftools consensus -f $(ref).fa $< | seqtk seq >$@
