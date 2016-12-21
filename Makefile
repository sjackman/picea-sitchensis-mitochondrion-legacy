# Assemble and annotate the Sitka spruce (Picea sitchensis) mitochondrion genome
# Written by Shaun Jackman @sjackman

# Target genome, Picea sitchensis mitochondrion
name=psitchensis
draft=psitchensismt

# Reference genome
ref=organelles
ref_fa=$(ref).fa
ref_gff=$(ref).gff

# Picea sitchensis plastid
psitchensiscp=KU215903

# Picea glauca plastid
pglaucacp=KT634228

# Picea glauca mitochondrion
pglaucamt=LKAM01

# Picea glauca nuclear genome
pglaucanuc=ALWZ04

# Number of threads
t=64

# gzip compression program
gzip=pigz -p$t

# Report run time and memory usage
export SHELL=zsh -opipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

all: assembly-stats.html

install-deps:
	brew install pigz
	brew tap homebrew/science
	brew install arcs bcftools bwa edirect fastqc miller samtools seqtk

.PHONY: all clean install-deps
.DELETE_ON_ERROR:
.SECONDARY:

# Entrez Direct

# Fetch data from NCBI.
$(name).fa $(ref).fa: %.fa:
	efetch -db nuccore -id $* -format fasta | seqtk seq | sed 's/^>/>$* /' >$@

# Download the Picea glauca mitochondrion FASTA.
LKAM.fa:
	curl ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/AM/LKAM01/LKAM01.1.fsa_nt.gz | gunzip -c | seqtk seq >$@

# Download the Picea glauca mitochondrion GBFF.
LKAM.gbff:
	curl ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/LK/AM/LKAM01/LKAM01.1.gbff.gz | gunzip -c >$@

# Download the Picea glauca nuclear FASTA.
ALWZ.fa.gz:
	curl ftp://ftp.ncbi.nih.gov/genomes/genbank/plant/Picea_glauca/representative/GCA_000411955.5_PG29_v4.1/GCA_000411955.5_PG29_v4.1_genomic.fna.gz >$@

# Split a FASTA file into 2 GB chunks of sequences.
%.900.seq %.901.seq %.902.seq %.903.seq %.904.seq %.905.seq %.906.seq %.907.seq %.908.seq %.909.seq %.910.seq %.911.seq %.912.seq %.913.seq %.914.seq %.915.seq %.916.seq %.917.seq %.918.seq %.919.seq %.920.seq %.921.seq %.922.seq %.923.seq %.924.seq %.925.seq %.926.seq %.927.seq %.928.seq %.929.seq %.930.seq %.931.seq %.932.seq %.933.seq %.934.seq %.935.seq %.936.seq %.937.seq %.938.seq %.939.seq %.940.seq %.941.seq %.942.seq %.943.seq %.944.seq %.945.seq: %.fa.gz
	gunzip -c $< | seqtk seq | sed '/^>/d' | split -a3 -C536870911 --numeric-suffixes=900 --additional-suffix=.seq - $*.

# Construct a concatemer FASTA file of a chunk of sequences.
ALWZ.%.concat.fa: ALWZ.%.seq
	echo '>$*' | cat - $< | seqtk seq >$@

# Construct a concatemer FASTA file of the Picea glauca nuclear genome.
ALWZ.concat.fa: ALWZ.900.concat.fa ALWZ.901.concat.fa ALWZ.902.concat.fa ALWZ.903.concat.fa ALWZ.904.concat.fa ALWZ.905.concat.fa ALWZ.906.concat.fa ALWZ.907.concat.fa ALWZ.908.concat.fa ALWZ.909.concat.fa ALWZ.910.concat.fa ALWZ.911.concat.fa ALWZ.912.concat.fa ALWZ.913.concat.fa ALWZ.914.concat.fa ALWZ.915.concat.fa ALWZ.916.concat.fa ALWZ.917.concat.fa ALWZ.918.concat.fa ALWZ.919.concat.fa ALWZ.920.concat.fa ALWZ.921.concat.fa ALWZ.922.concat.fa ALWZ.923.concat.fa ALWZ.924.concat.fa ALWZ.925.concat.fa ALWZ.926.concat.fa ALWZ.927.concat.fa ALWZ.928.concat.fa ALWZ.929.concat.fa ALWZ.930.concat.fa ALWZ.931.concat.fa ALWZ.932.concat.fa ALWZ.933.concat.fa ALWZ.934.concat.fa ALWZ.935.concat.fa ALWZ.936.concat.fa ALWZ.937.concat.fa ALWZ.938.concat.fa ALWZ.939.concat.fa ALWZ.940.concat.fa ALWZ.941.concat.fa ALWZ.942.concat.fa ALWZ.943.concat.fa ALWZ.944.concat.fa ALWZ.945.concat.fa
	cat $^ >$@

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

# Align paired-end reads to the target genome and sort.
$(ref).%.bam: %.fq.gz $(ref).fa.bwt
	bwa mem -t$t -p $(ref).fa $< | samtools view -h -F4 | samtools sort -@$t -o $@

# Align paired-end reads to the draft genome and do not sort.
$(abyss_scaffolds).%.sortn.bam: %.fq.gz $(abyss_scaffolds).fa.bwt
	bwa mem -t$t -p $(abyss_scaffolds).fa $< | samtools view -@$t -h -F4 -o $@

# LongRanger

# Index the target genome.
refdata-%/fasta/genome.fa.bwt: %.fa
	longranger mkref $<

# Lanes associated with this sample.
lane=4

# Extract barcodes using longranger basic.
$(name)_longranger_basic/outs/barcoded.fastq.gz: fastq_path/read-RA_si-TCAAGGCC_lane-004-chunk-006.fastq.gz
	longranger basic --id=$(name)_longranger_basic --fastqs=$(<D) --lanes=$(lane)

# Symlink the longranger basic FASTQ file.
$(name).longranger.basic.fq.gz: $(name)_longranger_basic/outs/barcoded.fastq.gz
	ln -sf $< $@

# Align reads to the target genome.
$(ref)_$(name)_longranger_align/outs/possorted_bam.bam: fastq_path/read-RA_si-TCAAGGCC_lane-004-chunk-006.fastq.gz refdata-$(ref)/fasta/genome.fa.bwt
	longranger align --id=$(ref)_$(name)_longranger_align --reference=refdata-$(ref) --fastqs=$(<D) --lanes=$(lane)

# Symlink the longranger align bam file.
$(ref).$(name).longranger.align.bam: $(ref)_$(name)_longranger_align/outs/possorted_bam.bam
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

# Extract those reads from a set of barcodes from a FASTQ file.
%.bx.atleast4.fq.gz: pglauca.%.longranger.align.bam.bx.atleast4.txt %.bx.fq.gz
	gunzip -c $*.bx.fq.gz | paste - - - - - - - - \
		| grep -Ff $< | tr '\t' '\n' | $(gzip) >$@

# Extract those reads from a set of barcodes from a BAM file.
%.bam.bx.atleast4.bam: %.bam.bx.atleast4.txt %.bam
	samtools view -h $*.bam | grep -Ff $< -e '@' | samtools view -b -o $@

# Calculate the genome coverage of a sorted BAM file
%.bam.coverage.tsv: %.bam
	samtools view -u -F260 $< | samtools depth - \
		| mlr --tsvlite --implicit-csv-header stats1 -a sum,count -f 3 \
		| mlr --tsvlite label Aligned,Covered \
		| mlr --tsvlite put '$$GenomeSize = $(GwithoutN); $$GenomeCoverage = $$Covered / $$GenomeSize' >$@

# Calculate the number of mismatches in a SAM file
%.sam.nm.tsv: %.sam
	samtools view -F260 $< | sed '/NM:i:/!d;s/^.*NM:i://;s/[[:space:]].*//' \
		| mlr --tsvlite --implicit-csv-header stats1 -a sum -f 1 \
		| mlr --tsvlite label NM >$@

# Summarize correctness and completeness
%.metrics.tsv: %.bam.coverage.tsv %.sam.nm.tsv
	paste $^ | mlr --tsvlite put '$$Identity = 1 - $$NM / $$Aligned; $$QV = -10 * log10(1 - $$Identity); $$File = "$*.sam"' >$@

# htsbox

# Convert BAM format to pairwise-alignment-format (PAF) using htsbox
%.paf: %.bam
	htsbox samview -p $< >$@

# Convert PAF format to TSV format
%.paf.tsv: %.paf
	(printf "qname\tqlength\tqstart\tqend\tstrand\ttname\ttlength\ttstart\ttend\tdivergence\tmapq\tattributes\n"; \
		awk 'NF == 12' $<) >$@

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
k=75
kc=4
B=100G
abyss_scaffolds=abyss/2.0.1/k$k/kc$(kc)/psitchensis-scaffolds

# The total genome size of P. sitchensis plastid and P. glauca mitochondrion
GwithN=6118703
GwithoutN=6055308
G=$(GwithoutN)

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

# Break scaffolds at gaps to produce scaftigs.
%-scaftigs.fa: %-scaffolds.fa
	seqtk seq $< | tr _ '~' | $(abyssbin201)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity statistics with abyss-fac.
%.stats.tsv: %.fa
	$(abyssbin201)/abyss-fac -t500 -G$G $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; abyss-samtobreak -l500 $<) >$@

# Convert samtobreak.txt to TSV.
%.samtobreak.tsv: %.samtobreak.txt
	abyss-samtobreak-to-tsv $< >$@

# ARCS

# Add the barcode to the read ID, and skip reads without barcodes.
%.bx.fq.gz: psitchensis.longranger.basic.fq.gz
	gunzip -c $< | gawk ' \
		{ bx = "NA" } \
		match($$0, "BX:Z:([ACGT]*)-1", x) { bx = x[1] } \
		bx == "NA" { getline; getline; getline; next } \
		{ print $$1 "_" bx " " $$2; getline; print; getline; print; getline; print }' \
		| $(gzip) >$@

# Create a graph of linked contigs using ARCS.
c=1
e=50000
r=0.220000
%.c$c_e$e_r$r.arcs.gv: %.sortn.bam $(abyss_scaffolds).fa
	bin/arcs -s98 -c$c -l0 -z500 -m4-20000 -d0 -e$e -r$r -v \
		-f $(abyss_scaffolds).fa -a <(echo $<) -b $*.c$c_e$e_r$r.arcs
	mv $*.c$c_e$e_r$r.arcs_original.gv $@

# Convert the ARCS graph to LINKS TSV format.
%.arcs.tsv: %.arcs.gv $(abyss_scaffolds).fa
	bin/arcs-makeTSVfile $< $@ $(abyss_scaffolds).fa

# Scaffold the assembly using the ARCS graph and LINKS.
a=0.999999
l=1
%.arcs.a$a_l$l.links.scaffolds.fa: %.arcs.tsv $(abyss_scaffolds).fa
	cp $< $*.arcs.a$a_l$l.links.tigpair_checkpoint.tsv
	LINKS -k20 -l$l -t2 -a$a -x1 -s /dev/null -f $(abyss_scaffolds).fa -b $*.arcs.a$a_l$l.links

# Rename the scaffolds.
%/arcs/psitchensis-scaffolds.fa: %/psitchensis-scaffolds.psitchensis.bx.atleast4.c$c_e$e_r$r.arcs.a$a_l$l.links.scaffolds.fa
	mkdir -p $(@D)
	gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' $< >$@

# Select scaffolds at least 5 kbp.
%.5kbp.fa: %.fa
	seqtk seq -L5000 $< >$@

# Spearmint

# Optimize the assembly parameters using Spearmint.
%/stamp: %.json %.py
	-if [ -r $(@D)/mongodb/pid ]; then kill `<$(@D)/mongodb/pid`; fi
	if [ -e $(@D) ]; then mv $(@D) $(@D)-`date '+%Y-%m-%dT%H:%M'`; fi
	mkdir -p $(@D)/mongodb
	mongod --fork --dbpath $(@D)/mongodb \
		--logpath $(@D)/mongodb/log --pidfilepath $(PWD)/$(@D)/mongodb/pid
	ln -sfn $(@D) output
	spearmint --config=$< . 2>&1 | tee $*/log
	touch $@

# Scrape the results of running Spearmint from its log files.
%.dkvp: %/log %/*.out
	grep -h NG50= $(<D)/*.out >$@

# Convert DKVP to TSV using Miller
%.tsv: %.dkvp
	mlr --tsvlite --idkvp put '$$Index=NR' <$< >$@

# QUAST

# Analayze the assembled genomes using QUAST.
%.quast/$(ref)/transposed_report.tsv: %.fa
	quast.py -t$t -o $(@D) -R $(ref_fa) -G $(ref_gff) -sLe --fragmented $<

# KAT

# Count k-mers in an assembly.
%.fa.jf: %.fa
	kat_jellyfish count -t$t -m40 -s50000000 -o $@ $<

# Count k-mers in reads that also occur in an assembly.
%.fq.gz.jf: $(name).fq.gz %.fa
	gunzip -c $< | kat_jellyfish count -t$t -m40 -s50000000 -o $@ --if=$*.fa /dev/stdin

# Analyze the joint k-mer spectrum of reads and an assembly.
%.kat.comp: %.fq.gz.jf %.fa.jf
	kat comp -t$t -o $@ $^

# Aggregate assembly statistics.

abyss-fac.tsv: \
		KU215903.stats.tsv \
		pglauca.stats.tsv \
		abyss/2.0.1/k64/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k80/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k83/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k84/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k85/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k88/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k92/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k98/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k100/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k102/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k104/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k108/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k112/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k128/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k83/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k84/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k85/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc2/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k68/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k72/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k74/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k76/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k79/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k80/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k82/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k83/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k84/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k85/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k86/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k87/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k88/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k89/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k112/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k128/kc3/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k74/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k75/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k76/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k78/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k79/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k80/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k81/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k83/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k84/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k85/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k86/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k87/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k88/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k89/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k90/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k91/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc4/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k65/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k66/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k67/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k68/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k69/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k70/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k72/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k73/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k74/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k75/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k76/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k77/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k78/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k79/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k80/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k81/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k82/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k83/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k84/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k85/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k86/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k94/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc5/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k48/kc6/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k56/kc6/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k60/kc6/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/kc6/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k68/kc6/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k96/kc10/psitchensis-scaffolds.stats.tsv \
		abyss/2.0.1/k64/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k80/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k84/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k88/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k92/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k98/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k100/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k102/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k104/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k108/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k112/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k128/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k64/kc2/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k83/kc2/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k84/kc2/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k85/kc2/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/kc2/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k64/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k68/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k72/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k74/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k76/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k79/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k80/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k82/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k83/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k84/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k85/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k86/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k87/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k88/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k89/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k112/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k128/kc3/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k64/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k74/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k75/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k76/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k78/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k79/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k80/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k81/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k83/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k84/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k85/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k86/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k87/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k88/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k89/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k90/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k91/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/kc4/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k64/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k65/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k66/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k67/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k68/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k69/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k70/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k72/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k73/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k74/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k75/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k76/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k77/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k78/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k79/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k80/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k81/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k82/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k83/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k84/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k85/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k86/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k94/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/kc5/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k48/kc6/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k56/kc6/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k60/kc6/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k64/kc6/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k68/kc6/psitchensis-scaftigs.stats.tsv \
		abyss/2.0.1/k96/kc10/psitchensis-scaftigs.stats.tsv
	mlr --tsvlite cat $^ >$@

samtobreak.tsv: \
		abyss/2.0.1/k64/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k80/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k83/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k84/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k85/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k88/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k92/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k96/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k98/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k100/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k102/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k104/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k108/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k112/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k128/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k64/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k68/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k72/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k74/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k76/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k79/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k80/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k82/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k83/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k84/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k85/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k86/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k87/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k88/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k89/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k96/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k112/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k128/kc3/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k64/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k74/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k75/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k76/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k78/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k79/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k80/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k81/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k83/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k84/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k85/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k86/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k87/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k88/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k89/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k90/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k91/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k96/kc4/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k64/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k65/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k66/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k67/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k68/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k69/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k70/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k72/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k73/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k74/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k75/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k76/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k77/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k78/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k79/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k80/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k81/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k82/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k83/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k84/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k85/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k86/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k94/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv \
		abyss/2.0.1/k96/kc5/$(ref).psitchensis-scaftigs.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

quast.tsv: \
		abyss/2.0.1/k96/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k64/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k80/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k84/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k88/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc3/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k84/kc4/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k85/kc4/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc4/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k68/kc5/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv \
		abyss/2.0.1/k96/kc5/psitchensis-scaffolds.quast/$(ref)/transposed_report.tsv
	mlr --tsvlite cat $^ >$@

# RMarkdown

# Render RMarkdown to HTML.
%.html: %.rmd
	Rscript -e 'rmarkdown::render("$<", params=list(input_tsv="$*.tsv"))'

# Dependencies

assembly-stats.html: abyss-fac.tsv samtobreak.tsv

spearmint-arcs.html: spearmint-arcs.tsv

spearmint-arcs-r.html: spearmint-arcs-r.tsv
