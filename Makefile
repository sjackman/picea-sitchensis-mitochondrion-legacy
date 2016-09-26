# Assemble and annotate the Sitka spruce (Picea sitchensis) mitochondrion genome
# Written by Shaun Jackman @sjackman

# Target genome, Picea sitchensis
name=psitchensis

# Reference genome
ref=organelles
ref_fa=$(ref).fa
ref_gff=$(ref).gff

# Picea sitchensis plastid
psitchensiscp=KU215903

# Picea glauca plastid
pglaucacp=KT634228

# Picea glauca mitochondrion
pglaucamt=LKAM00000000

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

# LongRanger

# Index the target genome.
refdata-%/fasta/genome.fa.bwt: %.fa
	longranger mkref $<

# Lanes associated with this sample.
lane=4

# Align reads to the target genome.
$(ref)_$(name)_longranger_align/outs/possorted_bam.bam: fastq_path/read-RA_si-TCAAGGCC_lane-004-chunk-006.fastq.gz refdata-$(ref)/fasta/genome.fa.bwt
	longranger align --id=$(ref)_$(name)_longranger_align --reference=refdata-$(ref) --fastqs=$(<D) --lanes=$(lane)

# Symlink the longranger align bam file.
$(ref).$(name).longranger.wgs.align.bam: $(ref)_$(name)_longranger_align/outs/possorted_bam.bam
	ln -sf $< $@

# Align reads to the target genome, call variants, phase variants, and create a Loupe file.
$(ref)_$(name)_longranger_wgs/outs/phased_possorted_bam.bam: fastq_path/read-RA_si-TCAAGGCC_lane-004-chunk-006.fastq.gz refdata-$(ref)/fasta/genome.fa.bwt
	longranger wgs --id=$(ref)_$(name)_longranger_wgs --sex=female --reference=refdata-$(ref) --fastqs=$(<D) --lanes=$(lane)

# Symlink the longranger wgs bam file.
$(ref).$(name).longranger.wgs.bam: $(ref)_$(name)_longranger_wgs/outs/phased_possorted_bam.bam
	ln -sf $< $@

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

# Convert BAM to FASTQ.
%.bam.fq.gz: %.bam
	samtools collate -Ou $< $@ | samtools fastq /dev/stdin | $(gzip) >$@

# Extract the alignment score.
%.bam.as.tsv: %.bam
	(echo AS; samtools view -F4 $< | gsed -r 's/.*AS:i:([0-9]*).*/\1/') >$@

# Extract the alignment score, barcode and molecule identifier.
%.bam.bx.tsv: %.bam
	samtools view -F4 $< | gawk -F'\t' ' \
		BEGIN { print "Flags\tRname\tMapq\tAS\tBX\tMI" } \
		{ as = bx = mi = "NA" } \
		match($$0, "AS:.:([^\t]*)", x) { as = x[1] } \
		match($$0, "BX:Z:([^\t]*)", x) { bx = x[1] } \
		match($$0, "MI:i:([^\t]*)", x) { mi = x[1] } \
		{ print $$2 "\t" $$3 "\t" $$5 "\t" as "\t" bx "\t" mi }' >$@

# Extract barcodes with at least 4 good aligned reads per barcode.
%.bam.bx.atleast4.txt: %.bam.bx.tsv
	awk 'NR > 1 && $$2 >= 4 {print $$1}' $< >$@

# Extract those reads from a set of barcodes from a BAM file.
%.bam.bx.atleast4.bam: %.bam.bx.atleast4.txt %.bam
	samtools view -h $*.bam | grep -Ff $< -e '@' | samtools view -b -o $@

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

# BFC

# Correct errors in reads.
bfc/%.fq.gz: %.fq.gz
	mkdir -p $(@D)
	bfc -t$t -s1G $< | tr '\t' ' ' | $(gzip) >$@

# ABySS
abyssbin190=/gsc/btl/linuxbrew/Cellar/abyss/1.9.0-k96/bin
abyssbin201=/gsc/btl/linuxbrew/Cellar/abyss/2.0.1-k96/bin
k=96
kc=3
B=100G

# The total genome size of P. sitchensis plastid and P. glauca mitochondrion
G=6055308

# Assemble reads with ABySS 1.9.0.
abyss/1.9.0/k$k/%-scaffolds.fa: pglauca.%.longranger.align.bam.bx.atleast4.bam.fq.gz
	mkdir -p $(@D)
	time $(abyssbin190)/abyss-pe -C $(@D) name=$* j=$t k=$k v=-v in=`realpath $<`

# Assemble reads with ABySS 2.0.1 with abyss-dbg.
abyss/2.0.1/k$k/%-scaffolds.fa: pglauca.%.longranger.align.bam.bx.atleast4.bam.fq.gz
	mkdir -p $(@D)
	time $(abyssbin201)/abyss-pe -C $(@D) name=$* j=$t k=$k v=-v in=`realpath $<`

# Assemble reads with ABySS 2.0.1 abyss-bloom-dbg.
abyss/2.0.1/k$k/kc$(kc)/%-scaffolds.fa: pglauca.%.longranger.align.bam.bx.atleast4.bam.fq.gz
	mkdir -p $(@D)
	time $(abyssbin201)/abyss-pe -C $(@D) name=$* j=$t k=$k kc=$(kc) B=$(B) v=-v in=`realpath $<`

# Calculate assembly contiguity statistics with abyss-fac.
%.stats.tsv: %.fa
	$(abyssbin201)/abyss-fac -t500 -G$G $< >$@

# QUAST

# Analayze the assembled genomes using QUAST.
%.quast/$(ref)/transposed_report.tsv: %.fa
	quast.py -t$t -o $(@D) -R $(ref_fa) -G $(ref_gff) -sLe --fragmented $<

# Aggregate assembly statistics.

abyss-fac.tsv: \
		KU215903.stats.tsv \
		pglauca.stats.tsv \
		abyss/2.0.1/k96/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc10/psitchensis-scaffolds.stats.tsv
	mlr --tsvlite cat $^ >$@

quast.tsv: \
		abyss/2.0.1/k96/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc4/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc5/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc10/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv
	mlr --tsvlite cat $^ >$@

# RMarkdown

# Render RMarkdown to HTML.
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<")'
