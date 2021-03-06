# Assemble and annotate the Sitka spruce (Picea sitchensis) mitochondrion genome
# Written by Shaun Jackman @sjackman

# Target genome, Picea sitchensis mitochondrion
name=psitchensis

# Assembly of Picea sitchensis organelles
#draft=psitchensiscpmt_x

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
	brew install arcs bcftools bedtools bwa edirect fastqc miller samtools seqtk

psitchensiscpmt_1:
	$(MAKE) psitchensiscpmt_1.fa draft=NA

psitchensiscpmt_2:
	$(MAKE) psitchensiscpmt_2.fa \
		draft=abyss/2.0.1/k75/kc4/psitchensis-scaffolds \
		c=1 e=50000 r=0.220000 a=0.999999 l=1

psitchensiscpmt_3:
	$(MAKE) psitchensiscpmt_3.fa draft=psitchensiscpmt_2

psitchensiscpmt_3_fix:
	$(MAKE) psitchensiscpmt_3.path.breakpoints.tigs.rename.fa draft=psitchensiscpmt_3.path

psitchensiscpmt_4:
	$(MAKE) psitchensiscpmt_4.fa draft=psitchensiscpmt_3

psitchensiscpmt_5:
	$(MAKE) psitchensiscpmt_5.fa draft=psitchensiscpmt_4

psitchensiscpmt_6:
	$(MAKE) psitchensiscpmt_6.fa draft=psitchensiscpmt_5

psitchensiscpmt_7:
	$(MAKE) psitchensiscpmt_7.fa draft=psitchensiscpmt_6

psitchensiscpmt_8:
	$(MAKE) psitchensiscpmt_8.fa draft=psitchensiscpmt_7

psitchensiscpmt_9:
	$(MAKE) psitchensiscpmt_9.fa draft=psitchensiscpmt_8

.PHONY: all clean install-deps \
	psitchensiscpmt_1 \
	psitchensiscpmt_2 \
	psitchensiscpmt_3 \
	psitchensiscpmt_3_fix \
	psitchensiscpmt_4 \
	psitchensiscpmt_5 \
	psitchensiscpmt_6 \
	psitchensiscpmt_7 \
	psitchensiscpmt_8 \
	psitchensiscpmt_9
.DELETE_ON_ERROR:
.SECONDARY:

# Symlink the draft genome.

psitchensismt.1.fa: abyss/2.0.1/k75/kc4/arcs/psitchensis-scaffolds.nocp.2kbp.fa
	ln -sf $< $@

psitchensismt_2.fa: abyss/2.0.1/k75/kc4/arcs/psitchensis-scaffolds.nocp.2500bp.fa
	ln -sf $< $@

# Combine the plastome and chondrome into a single reference
psitchensiscpmt_2.fa: KU215903.fa psitchensismt_2.fa
	cat $^ >$@

# Break scaffolds at loci not supported by molecules.
psitchensiscpmt_3.fa: psitchensiscpmt_2.breakpoints.tigs.2500bp.fa
	ln -sf $< $@

# Correct a misassembly in scaffold 3-7 and rename two scaffolds.
psitchensiscpmt_3.path.breakpoints.tigs.rename.bed: psitchensiscpmt_3.path.breakpoints.tigs.bed
	sed 's/3-7-1/4011/;s/3-7-2/3-7/' $< >$@

# Generate a FASTA file of the corrected misassembly in scaffold 3-7.
psitchensiscpmt_3.path.breakpoints.tigs.rename.fa: psitchensiscpmt_3.path.breakpoints.tigs.rename.bed psitchensiscpmt_3.path.fa
	bedtools getfasta -name -fi psitchensiscpmt_3.path.fa -bed $< | sed 's/::/ /;s/^NN*//;s/NN*$$//' >$@

# Symlink the revised draft genome.
psitchensiscpmt_4.fa: psitchensiscpmt_3.path.breakpoints.tigs.rename.fa
	ln -sf $< $@

# Symlink the revised draft genome.
psitchensiscpmt_5.fa: psitchensiscpmt_4.breakpoints.tigs.fa
	ln -sf $< $@

# Symlink the revised draft genome.
psitchensiscpmt_6.fa: psitchensiscpmt_5.path.fa
	ln -sf $< $@

# Symlink the revised draft genome.
psitchensiscpmt_7.fa: psitchensiscpmt_6.breakpoints.tigs.fa
	ln -sf $< $@

# Generate the FASTA file of the scaffolds.
%.path.fa: %.fa %.fa.fai %.path
	MergeContigs -v -k$k -o $@ $^

# Symlink the revised draft genome.
psitchensiscpmt_8.fa: psitchensiscpmt_7.path.fa
	ln -sf $< $@

psitchensiscpmt_9.fa: psitchensiscpmt_8/8003.bed.bx.as100.bam.barcodes.bx.unicycler.l50k.fa psitchensiscpmt_8.fa
	(seqmagick convert --line-wrap=0 --pattern-exclude '^(41|55|4005-2|8003)$$' psitchensiscpmt_8.fa -; \
		seqtk rename $< 900 $<) >$@

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

# Filter a FASTA file by length.
%.l50k.fa: %.fa
	seqtk seq -L50000 $< >$@

# Interleave paired-end reads
%.fq.gz: %_1.fq.gz %_2.fq.gz
	seqtk mergepe $^ | $(gzip) >$@

# Calculate the nucleotide composition of FASTA
%.fa.comp.tsv: %.fa
	(printf "Name\tLength\tNumA\tNumC\tNumG\tNumT\n"; seqtk comp $< | cut -f1-6) | mlr --tsvlite put '$$NumACGT = $$NumA + $$NumC + $$NumG + $$NumT; $$GC = ($$NumC + $$NumG) / $$NumACGT' >$@

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
	bwa mem -t$t -pC $(ref).fa $< | samtools view -h -F4 | samtools sort -@$t -o $@

# Align paired-end reads to the draft genome and do not sort.
$(draft).%.sortn.bam: %.fq.gz $(draft).fa.bwt
	bwa mem -t$t -pC $(draft).fa $< | samtools view -@$t -h -F4 -o $@

# Align paired-end reads to the draft genome.
$(draft)/%.bwa.bam: $(draft)/%.fq.gz $(draft).fa.bwt
	bwa mem -t$t -pC $(draft).fa $< | samtools view -h -F4 | samtools sort -@$t -o $@

# Align sequences to the draft genome.
$(draft)/%.fa.bwa.sam: $(draft)/%.fa $(draft).fa.bwt
	bwa mem -t$t -xintractg $(draft).fa $< >$@

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

# Convert bedpe to TSV.
%.bedpe.tsv: %_longranger_wgs/outs/large_sv_calls.bedpe %_longranger_wgs/outs/large_sv_candidates.bedpe
	(printf "chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\tname\tqual\tstrand1\tstrand2\tfilters\tinfo\n"; \
		sed '/^#/d' $^) >$@

# Convert bedpe.tsv to GraphViz.
%.bedpe.gv: %.bedpe.tsv $(name)cpmt_2.fa.fai
	( echo "graph "$*" {"; \
		mlr --inidx --ifs tab put -q 'print $$1. " [l=" . $$2 . "]" . " [label=\"" . $$1 . "\\n" . $$2 . " bp\"]"' <$(name)cpmt_2.fa.fai; \
		mlr --tsvlite stats1 -a max -g chrom1,chrom2 -f qual \
			then put -q 'print $$chrom1 . " -- ". $$chrom2 . " [weight=" . $$qual_max . "]" . " [label=" . $$qual_max . "]"' $<; \
		echo "}"; \
	) >$@

# Count the number of reads per barcode.
%.fq.gz.count.tsv: %.fq.gz
	gunzip -c $< | sed -n 's/.*BX:Z:/BX=/p' \
	| mlr --otsvlite count-distinct -f BX then rename count,Reads then sort -f BX >$@

# Create a SAM header of barcodes.
%.count.tsv.sam: %.count.tsv
	mlr --tsvlite put -q 'print "@SQ\tSN:" . $$BX . "\tLN:1"' $< >$@

# Create an unaligned BAM file indexed by barcode.
# The flags assume that the reads are interleaved with no missing reads.
%.longranger.basic.bam: %.longranger.basic.fq.gz %.longranger.basic.fq.gz.count.tsv.sam
	( cat $*.longranger.basic.fq.gz.count.tsv.sam; \
	gunzip -c $< | gawk ' \
		{ bx = "NA"; rname = "*" } \
		match($$0, "BX:Z:([ACGT]*-[0-9])", x) { bx = x[1]; rname = bx } \
		{ id = substr($$1, 2); getline; seq = $$0; getline; getline; qual = $$0 } \
		{ which = !which; flags = which ? 77 : 141 } \
		{ printf "%s\t%u\t%s\t1\t0\t*\t*\t0\t0\t%s\t%s\tBX:Z:%s\n", id, flags, rname, seq, qual, bx }' \
	) | samtools sort -@$t -T/var/tmp/$(USER)/$@ -O bam -o $@

# GraphViz

# Extract the largest connected component from a graph using ccomps.
%.comp1.gv: %.gv
	ccomps -zX'#0' -o $@ $< || test -r $@

# Filter scaffolds by length using gvpr.
%.l5k.gv: %.gv
	gvpr -i 'N[l >= 5000]' -o $@ $<

# Filter scaffolds by length using gvpr.
%.l10k.gv: %.gv
	gvpr -i 'N[l >= 10000]' -o $@ $<

# Filter scaffolds by length using gvpr.
%.l20k.gv: %.gv
	gvpr -i 'N[l >= 20000]' -o $@ $<

# Filter edges by number of barcodes using gvpr.
%.m5.gv: %.gv
	gvpr 'E[label >= 5]' -o $@ $<

# Filter edges by their attribute n.
%.n$n.gv: %.gv
	gvpr 'E[n >= $n]' -o $@ $<

# Filter vertices by their degree.
%.deg3.gv: %.gv
	gvpr -i 'N[degree >= 3]' $< >$@

# Filter edges by their attribute q.
q=0.05
%.q$q.gv: %.gv
	gvpr 'E[q < $q]' -o $@ $<

# Filter edges by their boolean attribute best.
%.best.gv: %.gv
	gvpr 'E[best == "T"]' $< | sed 's/best=T,//' >$@

# Render a graph to PNG using dot.
%.gv.dot.png: %.gv
	dot -Tpng -o $@ $<

# Render a graph to PDF using dot.
%.gv.dot.pdf: %.gv
	dot -Tpdf -o $@ $<

# Render a graph to PDF using fdp.
%.gv.fdp.pdf: %.gv
	gvpr -c 'E{weight = n / 100.0}' $< | fdp -Tpdf -o $@

# Render a graph to PNG using dot and gvpack.
%.gv.gvpack.dot.pdf: %.gv
	sed 's/q=0.000000"/"/;s/q=0"/"/' $< | (ccomps -x || true) | dot | gvpack -g | neato -s -n2 -Tpdf -o $@

# Render a graph to PDF using neato.
%.gv.neato.pdf: %.gv
	neato -Tpdf -Goverlap=false -Gsplines=true -Gsep=0.5 -o $@ $<

# Render a graph to PNG using neato.
%.gv.neato.png: %.gv
	neato -Tpng -Goverlap=false -Gsplines=true -Gsep=0.5 -o $@ $<

# Render a graph to SVG using neato.
%.gv.neato.svg: %.gv
	neato -Tsvg -Goverlap=false -Gsplines=true -Gsep=0.5 -o $@ $<

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Sort a SAM file and produce a sorted BAM file.
%.bam: %.sam
	samtools sort -@$t -o $@ $<

# Sort a query-name-sorted BAM file by target.
%.bam: %.sortn.bam
	samtools sort -@$t -o $@ $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index $<

# Select properly paired reads.
%.proper.bam: %.bam
	samtools view -Obam -f2 -o $@ $<

# Convert BAM to FASTQ.
%.bam.fq.gz: %.bam
	samtools collate -Ou $< $@ | samtools fastq /dev/stdin | $(gzip) >$@

# Extract the alignment score.
%.bam.as.tsv: %.bam
	(echo AS; samtools view -F4 $< | gsed -r 's/.*AS:i:([0-9]*).*/\1/') >$@

# Extract the alignment score, barcode and molecule identifier.
%.bam.bx.tsv: %.bam
	samtools view -F4 $< | gawk -F'\t' ' \
		BEGIN { print "Flags\tRname\tPos\tMapq\tAS\tNM\tBX\tMI" } \
		{ as = bx = mi = nm = "NA" } \
		match($$0, "AS:.:([^\t]*)", x) { as = x[1] } \
		match($$0, "NM:.:([^\t]*)", x) { nm = x[1] } \
		match($$0, "BX:Z:([^\t]*)", x) { bx = x[1] } \
		match($$0, "MI:i:([^\t]*)", x) { mi = x[1] } \
		{ print $$2 "\t" $$3 "\t" $$4 "\t" $$5 "\t" as "\t" nm "\t" bx "\t" mi }' >$@

# Count the number of reads per barcode.
%.bam.bx.count.tsv: %.bam.bx.tsv
	mlr --tsvlite filter '($$Flags & 2) != 0 && $$AS >= 40 && $$BX != "NA"' \
		then count-distinct -f BX \
		then rename count,Reads \
		then sort -nr Reads -f BX \
		$< >$@

# Extract barcodes with at least 4 good aligned reads per barcode.
%.bam.bx.atleast4.txt: %.bam.bx.count.tsv
	mlr --tsvlite --headerless-csv-output filter '$$Reads >= 4' then cut -f BX $< >$@

# Extract those reads from a set of barcodes from a FASTQ file.
%.bx.atleast4.fq.gz: pglauca.%.longranger.align.bam.bx.atleast4.txt %.bx.fq.gz
	gunzip -c $*.bx.fq.gz | paste - - - - - - - - \
		| grep -Ff $< | tr '\t' '\n' | $(gzip) >$@

# Extract those reads from a set of barcodes from a FASTQ file.
$(ref).%.bx.bam.atleast4.fq.gz: $(ref).%.bx.bam.bx.atleast4.txt %.bx.fq.gz
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
%.bam.nm.tsv: %.bam
	samtools view -F260 $< | sed '/NM:i:/!d;s/^.*NM:i://;s/[[:space:]].*//' \
		| mlr --tsvlite --implicit-csv-header stats1 -a sum -f 1 \
		| mlr --tsvlite label NM >$@

# Summarize correctness and completeness
%.bam.metrics.tsv: %.bam.coverage.tsv %.bam.nm.tsv
	paste $^ | mlr --tsvlite put '$$Identity = 1 - $$NM / $$Aligned; $$QV = -10 * log10(1 - $$Identity); $$File = "$*.sam"' >$@

# Remove secondary, supplementary and mapq=0 alignments.
%.primary.bam: %.bam
	samtools view -b -F0x900 -q1 -o $@ $<

# Remove alignments with an alignment score below a threshold.
%.as100.bam: %.bam
	samtools view -h -F4 $< | gawk -F'\t' ' \
			/^@/ { print; next } \
			{ as = 0 } \
			match($$0, "AS:.:([^\t]*)", x) { as = x[1] } \
			as >= 100' \
		| samtools view -@$t -b -o $@

# Select alignments with number of mismatches below a threshold.
nm=5
%.nm$(nm).bam: %.bam
	samtools view -h -F4 $< \
		| gawk -F'\t' ' \
			/^@/ { print; next } \
			{ nm = 999999999 } \
			match($$0, "NM:i:([^\t]*)", x) { nm = x[1] } \
			nm < $(nm)' \
		| samtools view -@$t -b -o $@

# Extract the query and target name from a BAM file.
%.bam.names.tsv: %.bam
	(printf "QName\tTName\n"; samtools view $< | cut -f1,3) >$@

# Filter reads by alignment score.
%.as-30.bx.tsv: %.bx.tsv
	mlr --tsvlite filter '$$AS >= -30' $< >$@

# Filter reads by alignment score.
%.as50.bx.tsv: %.bx.tsv
	mlr --tsvlite filter '$$AS >= 50' $< >$@

# Filter reads by alignment score.
%.as100.bx.tsv: %.bx.tsv
	mlr --tsvlite filter '$$AS >= 100' $< >$@

# Filter alignments by number of mismatches.
%.nm$(nm).bx.tsv: %.bx.tsv
	mlr --tsvlite filter '$$NM < $(nm)' $< >$@

# Group reads into molecules and add molecule identifiers.
%.bam.mi.bx.tsv: %.bam.bx.tsv
	./mi.r $< $@

# Create a TSV file of molecule extents.
%.bx.molecule.tsv: %.bx.tsv
	mlr --tsvlite \
		then stats1 -g BX,MI,Rname -a count,min,p50,max -f Pos,Mapq,AS,NM \
		then rename Pos_min,Start,Pos_max,End,Mapq_p50,Mapq_median,AS_p50,AS_median,NM_p50,NM_median,Pos_count,Reads \
		then put '$$Size = $$End - $$Start' \
		then cut -o -f Rname,Start,End,Size,BX,MI,Reads,Mapq_median,AS_median,NM_median \
		then filter '$$Reads >= 4' \
		$< >$@

# Create a BED file of molecule extents.
%.bx.molecule.bed: %.bx.molecule.tsv
	mlr --tsvlite --headerless-csv-output \
		put '$$Start = $$Start - 1; $$End = $$End - 1' \
		then put '$$Name = "Reads=" . $$Reads . ",Size=" . $$Size . ",Mapq=" . $$Mapq_median . ",AS=" . $$AS_median . ",NM=" . $$NM_median . ",BX=" . $$BX . ",MI=" . $$MI' \
		then cut -o -f Rname,Start,End,Name,Reads $< >$@

# Report summary statistics of a Chromium library
%.bx.molecule.summary.html: %.bx.molecule.tsv
	Rscript -e 'rmarkdown::render("summary.rmd", "html_document", "$@", params = list(input_tsv="$<", output_tsv="$*.summary.tsv"))'

# Identify misassemblies in psitchensiscpmt_2

# Calculate molecule depth of coverage.
psitchensiscpmt_2.breakpoints.tsv: %.breakpoints.tsv: %.psitchensis.longranger.wgs.bam.bx.molecule.tsv %.psitchensis.longranger.wgs.bam.as-30.bx.molecule.tsv
	Rscript -e 'rmarkdown::render("molecules.rmd", "html_notebook", "$*.breakpoints.nb.html", params = list(raw_tsv="$<", filtered_tsv="$*.psitchensis.longranger.wgs.bam.as-30.bx.molecule.tsv", output_tsv="$@"))'

# Identify misassemblies

# Calculate molecule depth of coverage.
%.psitchensis.longranger.align.bam.as100.nm5.breakpoints.tsv: %.psitchensis.longranger.align.bam.bx.molecule.tsv %.psitchensis.longranger.align.bam.as100.nm5.bx.molecule.tsv
	Rscript -e 'rmarkdown::render("molecules.rmd", "html_notebook", "$*.psitchensis.longranger.align.bam.as100.nm5.breakpoints.nb.html", params = list(raw_tsv="$<", filtered_tsv="$*.psitchensis.longranger.align.bam.as100.nm5.bx.molecule.tsv", output_tsv="$@"))'

# Exploratory data analysis of molecule size.
%.bx.as100.nm5.bam.mi.bx.molecules.nb.html: %.bx.as100.nm5.bam.mi.bx.molecule.tsv %.bx.bam.mi.bx.molecule.tsv
	Rscript -e 'rmarkdown::render("molecules.rmd", "html_notebook", "$@", params = list(raw_tsv="$*.bx.bam.mi.bx.molecule.tsv", filtered_tsv="$<", output_tsv="$*.bx.as100.nm5.bam.mi.bx.molecule.breakpoints.tsv"))'

# Determine coordinates of subsequences.
%.breakpoints.tigs.bed: %.breakpoints.tsv $(draft).fa.fai
	Rscript -e 'rmarkdown::render("breaktigs.rmd", "html_notebook", "$*.breakpoints.tigs.nb.html", params = list(input_tsv="$<", input_fai="$(draft).fa.fai", output_bed="$@"))'

# Break scaffolds at loci not supported by molecules.
%.breakpoints.tigs.fa: %.breakpoints.tigs.bed $(draft).fa
	bedtools getfasta -name -fi $(draft).fa -bed $< | sed 's/::/ /;s/^NN*//;s/NN*$$//' >$@

# Identify breakpoints

# Select BED records of a given size or larger.
size_threshold=500
%.size$(size_threshold).bed: %.bed
	awk '$$3 - $$2 >= $(size_threshold)' $< >$@

# Count start positions of molecules larger than a threshold size.
starts_size_threshold=2000
%.molecule.starts.tsv: %.molecule.tsv
	mlr --tsvlite \
		then filter '$$Size >= $(starts_size_threshold)' \
		then count-distinct -f Rname,Start \
		then rename Start,Pos,count,Starts \
		then sort -f Rname -n Pos \
		$< >$@

# Join the tables of depth of coverage and number of molecule starts.
%.molecule.size$(size_threshold).depth.starts.tsv: %.molecule.size$(size_threshold).bed.depth.tsv %.molecule.starts.tsv
	mlr --tsvlite join -u -j Rname,Pos -f $^ >$@

# Select positions with low depth of coverage and high numer of molecule starts.
depth_threshold=75
starts_threshold=4
pos_threshold=200
%.depth.starts.breakpoints.tsv: %.depth.starts.tsv
	mlr --tsvlite filter '$$Depth < $(depth_threshold) && $$Starts >= $(starts_threshold) && $$Pos >= $(pos_threshold)' $< >$@

# Identify breakpoints with low depth of coverage and high number of molecule starts.
%.depth.starts.breakpoints.tsv: %.size$(size_threshold).bed.depth.tsv %.starts.tsv
	Rscript -e 'rmarkdown::render("breakpoints.rmd", "html_notebook", "$*.depth.starts.breakpoints.nb.html", params = list(depth_tsv="$<", starts_tsv="$*.starts.tsv", depth_starts_tsv="$*.depth.starts.tsv", breakpoints_tsv="$@"))'

# Break scaffolds at regions of low long read coverage.

%.bedgraph.breakpoints.tsv: $(draft).fa.fai %.bedgraph
	gawk 'BEGIN { print "Rname\tPos" } \
		ARGIND == 1 { len[$$1] = $$2; next } \
		$$4 == 0 && $$2 > 0 && $$3 < len[$$1] { print $$1 "\t" 1+$$2 "\n" $$1 "\t" 1+$$3 }' $^ >$@

# igvtools

# Index a BED file.
%.bed.idx: %.bed
	igvtools index $<

# bedtools

# Convert BED to BAM.
%.bed.bam: %.bed $(draft).fa.fai
	awk '$$2 != $$3' $< | bedtools bedtobam -i - -g $(draft).fa.fai | samtools sort -@$t -Obam -o $@

# Compute statistics on the depth of coverage of a BED file.
%.bed.genomecov.tsv: %.bed $(draft).fa.fai
	(printf "Rname\tDepth\tCount\tRsize\tFraction\n"; awk '$$2 != $$3' $< | bedtools genomecov -g $(draft).fa.fai -i -) >$@

# Compute the depth of coverage of a BAM file.
%.bam.depth.tsv: %.bam
	(printf "Rname\tPos\tDepth\n"; bedtools genomecov -d -ibam $<) >$@

# Compute the depth of coverage of a BED file.
%.bed.depth.tsv: %.bed $(draft).fa.fai
	(printf "Rname\tPos\tDepth\n"; awk '$$2 != $$3' $< | bedtools genomecov -d -g $(draft).fa.fai -i -) >$@

# Calculate depth of coverage statistics.
%.depth.stats.tsv: %.depth.tsv
	mlr --tsvlite stats1 -a count,p25,p50,p75,mean,stddev -f Depth $< >$@

# Calculate depth of coverage statistics from bedtools genomecov.
%.genomecov.stats.tsv: %.genomecov.tsv
	mlr --tsvlite \
		then filter '$$Rname == "genome" && $$Depth > 0' \
		then step -a rsum -f Fraction \
		then put -q '@Depth_count += $$Count; if (is_null(@p25) && $$Fraction_rsum >= 0.25) { @p25 = $$Depth }; if (is_null(@p50) && $$Fraction_rsum >= 0.50) { @p50 = $$Depth }; if (is_null(@p75) && $$Fraction_rsum >= 0.75) { @p75 = $$Depth } end { emitf @Depth_count, @p25, @p50, @p75 }' \
		then rename p25,Depth_p25,p50,Depth_p50,p75,Depth_p75 \
		then put '$$Depth_IQR = $$Depth_p75 - $$Depth_p25' \
		$< >$@

# Calculate depth of coverage statistics per sequence.
%.depth.seqstats.tsv: %.depth.tsv
	mlr --tsvlite stats1 -a count,p25,p50,p75,mean,stddev -f Depth -g Rname $< >$@

# htsbox

# Convert BAM format to pairwise-alignment-format (PAF) using htsbox
%.paf: %.bam
	htsbox samview -p $< >$@

# Convert PAF format to TSV format
%.paf.tsv: %.paf
	(printf "Qname\tQlength\tQstart\tQend\tStrand\tTname\tTlength\tTstart\tTend\tMatches\tLength\tMapq\tAttributes\n"; \
		cat $<) >$@

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

# The total genome size of P. sitchensis plastid and P. glauca mitochondrion
GwithN=6118703
GwithoutN=6055308
G=$(GwithoutN)

# Assemble reads with ABySS 2.0.1 ABYSS.
abyss_k=64
abyss_e=5
abyss_c=5
%.abyss.k$(abyss_k).e$(abyss_e).c$(abyss_c).fa: %.fq.gz
	$(abyssbin201)/ABYSS -v -k$(abyss_k) -e$(abyss_e) -c$(abyss_c) -o $@ $<

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

# Break scaffolds at gaps to produce scaftigs.
%.scaftigs.fa: %.fa
	seqtk seq $< | tr _ '~' | $(abyssbin201)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity statistics with abyss-fac.
%.stats.tsv: %.fa
	abyss-fac -t500 -G$G $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; abyss-samtobreak -l500 $<) >$@

# Convert samtobreak.txt to TSV.
%.samtobreak.tsv: %.samtobreak.txt
	abyss-samtobreak-to-tsv $< >$@

# ARCS

# Add the barcode to the read ID, and skip reads without barcodes.
%.bx.fq.gz: %.longranger.basic.fq.gz
	gunzip -c $< | gawk ' \
		{ bx = "NA" } \
		match($$0, "BX:Z:([ACGT]*)-1", x) { bx = x[1] } \
		bx == "NA" { getline; getline; getline; next } \
		{ print $$1 "_" bx " " $$2; getline; print; getline; print; getline; print }' \
		| $(gzip) >$@

# Create a graph of linked contigs using ARCS.
c=1
e=30000
r=0.200000
%.c$c_e$e_r$r.arcs_original.gv %.c$c_e$e_r$r.arcs.dist.gv %.c$c_e$e_r$r.arcs.dist.tsv: %.sortn.bam $(draft).fa
	bin/arcs -s98 -c$c -l0 -z500 -m4-20000 -d0 -e$e -r$r -v \
		-f $(draft).fa \
		-b $*.c$c_e$e_r$r.arcs \
		-g $*.c$c_e$e_r$r.arcs.dist.gv \
		--tsv=$*.c$c_e$e_r$r.arcs.dist.tsv \
		--barcode-counts=$<.barcode-counts.tsv \
		$<

# Convert the ARCS graph to LINKS TSV format.
%.arcs.links.tsv: %.arcs_original.gv $(draft).fa
	bin/arcs-makeTSVfile $< $@ $(draft).fa

# Scaffold the assembly using the ARCS graph and LINKS.
a=0.999999
l=10
%.arcs.a$a_l$l.links.scaffolds.fa %.arcs.a$a_l$l.links.assembly_correspondence.tsv: %.arcs.links.tsv $(draft).fa
	cp $< $*.arcs.a$a_l$l.links.tigpair_checkpoint.tsv
	LINKS -k20 -l$l -t2 -a$a -x1 -s /dev/null -f $(draft).fa -b $*.arcs.a$a_l$l.links

# Rename the scaffolds.
%/arcs/psitchensis-scaffolds.fa: %/psitchensis-scaffolds.psitchensis.bx.atleast4.c$c_e$e_r$r.arcs.a$a_l$l.links.scaffolds.fa
	mkdir -p $(@D)
	gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' $< >$@

# Convert the LINKS assembly_correspondence.tsv file to GraphViz format.
%.links.path.gv: %.links.assembly_correspondence.tsv
	bin/links-correspondence-to-gv $< >$@

# Combine the ARCS and LINKS graphs
%.arcs.a$a_l$l.links.dist.path.gv: %.arcs.dist.n$l.gv %.arcs.a$a_l$l.links.path.gv
	bin/graph-union-strict $^ >$@

# Add p-values to an ARCS arcs.dist.tsv file
%.arcs.dist.p.tsv: %.arcs.dist.tsv enrichment.rmd
	Rscript -e 'rmarkdown::render("enrichment.rmd", "html_document", "$*.arcs.dist.html", params = list(input_tsv="$<", output_tsv="$@"))'

# Convert an ARCS dist.p.tsv file to GraphViz format
%.dist.p.gv: %.dist.p.tsv $(draft).fa.fai
	( echo 'digraph g {'; \
		abyss-todot $(draft).fa.fai \
			| gvpr -c 'N{label = sprintf("%s\\n%u bp", name, l)}' \
			| sed '1d;$$d'; \
		mlr --tsvlite put -q 'print "\"" . $$U . "\" -> \"" . $$V . "\" [ best=" . $$Best_orientation . " n=" . $$Shared_barcodes . " q=" . $$q . " label=\"n=" . $$Shared_barcodes . "\\nq=" . $$q . "\" ]"' $<; \
		echo '}' ) >$@

# Convert an ARCS GraphViz file to GFA format.
%.gv.gfa: %.gv
	(printf 'H\tVN:Z:1.0\n'; gvpr 'N{print("S\t", name, "\t*\tLN:i:", l)}' $< | gsed -n 's/+\t/\t/p') \
		| abyss-todot -e --gfa - $< | sed 's/*$$/0M/' >$@

# Tigmint

# Create a sequence segment graph using Tigmint-ARCS.
tigmint_e=10000
%.c$c_e$e.arcs_original.gv %.c$c_e$(tigmint_e).tigmint.arcs.dist.gv %.c$c_e$(tigmint_e).tigmint.arcs.dist.tsv: %.sortn.bam $(draft).fa
	bin/tigmint-arcs -s98 -c$c -l0 -z500 -m4-20000 -d0 -e$(tigmint_e) -v \
		-f $(draft).fa \
		-b $*.c$c_e$(tigmint_e).tigmint.arcs \
		-g $*.c$c_e$(tigmint_e).tigmint.arcs.dist.gv \
		--tsv=$*.c$c_e$(tigmint_e).tigmint.arcs.dist.tsv \
		--barcode-counts=$<.barcode-counts.tsv \
		$<

# Colour the segments by the scaffold from which they originate.
%.colours.tsv: %.psitchensis.bx.c$c_e$(tigmint_e).tigmint.arcs.dist.tsv
	Rscript -e 'rmarkdown::render("colours.rmd", "html_document", "$*.colours.html", params = list(input_tsv="$<", output_tsv="$@"))'

# Convert an ARCS dist.p.tsv file to GraphViz format with colours.
%.dist.p.colour.gv: %.dist.p.tsv $(draft).fa.fai $(draft).colours.tsv
	( echo 'strict graph g {'; \
		echo 'node[style = "filled"]'; \
		mlr --tsvlite put -q 'print "\"" . $$Scaffold . ":" . $$Position . "\" [ fillcolor=\"" . $$Colour . "\" ]"' $(draft).colours.tsv; \
		abyss-todot $(draft).fa.fai \
			| gvpr -c 'N{label = sprintf("%s\\n%u bp", name, l)}' \
			| sed '1d;$$d'; \
		mlr --tsvlite put -q 'print "\"" . $$U . "\" -- \"" . $$V . "\" [ best=" . $$Best_orientation . " n=" . $$Shared_barcodes . " q=" . $$q . " label=\"n=" . $$Shared_barcodes . "\\nq=" . $$q . "\" ]"' $< | grep -v KU215903; \
		echo '}' ) >$@

# Add d and e records and remove the label record for abyss-scaffold.
# Remove the attribute C=0.
%.abyss.dist.gv: %.gv
	sed 's/label="[^"]*",$$/d=100, e=100,/;s/C=0,//' $< >$@

# ABySS DistanceEst
de_l=100
de_s=1000
de_n=5

%.l$(de_l).bam %.l$(de_l).hist: %.bam $(draft).fa
	samtools view -h -F0x900 $< | abyss-fixmate -v -l$(de_l) -h $*.l$(de_l).hist | samtools sort -@$t -Obam -o $*.l$(de_l).bam

%.l$(de_l).n$(de_n).dist.gv: %.l$(de_l).bam %.l$(de_l).hist
	samtools view -h $< | DistanceEst -v --dot -j$t -k$k -l$(de_l) -s$(de_s) -n$(de_n) --maxd=1000 -o $@ $*.l$(de_l).hist

# Label the nodes and edges of an ABySS GraphViz file.
%.label.gv: %.gv $(draft).fa.fai
	( echo 'digraph g {'; \
		abyss-todot $(draft).fa.fai \
			| gvpr -c 'N{label = sprintf("%s\\n%u bp", name, l)}' \
			| sed '1d;$$d'; \
		gvpr -c 'E{label = sprintf("d=%d\\nn=%u", d, n)}' $< | sed '1d;$$d'; \
		echo '}' ) | gvpr 'E[1]' >$@

# Combine the paired-end and 10x distance estimates.
%.l$(de_l).n$(de_n).arcs+pe.dist.gv: $(draft).fa.fai $(draft).psitchensis.bx.sortn.l$(de_l).n$(de_n).dist.gv %.gv
	abyss-todot -v -e $^ >$@

# ABySS-Scaffold

# Scaffold using different values of s.
%.dist.gv.abyss-scaffold.n$n.s.txt: $(draft).fa.fai %.dist.gv
	for s in {1..50}000; do \
		abyss-scaffold -k$k -s$$s -n$n -G$G $^ >/dev/null |& tail -n2; \
	done >$@

# Scaffold using different values of n.
%.dist.gv.abyss-scaffold.txt: $(draft).fa.fai %.dist.gv
	for n in {1..20}; do \
		echo n=$$n; \
		abyss-scaffold -k$k -s$s -n$$n -G$G $^ >/dev/null |& tail -n3; \
	done >$@

# Scaffold the assembly using the ARCS graph and abyss-scaffold.
s=500-50000
n=10
%.dist.gv.n$n.abyss-scaffold.path: $(draft).fa.fai %.dist.gv
	abyss-scaffold -v -k$k -s$s -n$n -G$G -o $@ $^

# Generate the FASTA file of the scaffolds.
%.n$n.abyss-scaffold.fa: $(draft).fa $(draft).fa.fai %.n$n.abyss-scaffold.path
	MergeContigs -v -k$k -o $@ $^

# Convert the path file to GraphViz format.
%.path.gv: %.path
	perl -n \
		-e 'chomp;' \
		-e 'sub rc($$) { my $$s = shift; my $$c = chop $$s; $$s . ($$c eq "+" ? "-" : $$c eq "-" ? "+" : $$c) }' \
		-e 'my ($$id, @path) = split;' \
		-e 'print $$id, "+\t", join(" ", @path), "\n";' \
		-e 'for $$i (@path) { $$i = rc($$i) };' \
		-e 'print $$id, "-\t", join(" ", reverse @path), "\n";' \
		$< \
	| { echo 'digraph g {'; \
		gsed -e 's/ [^ ]*N//g' -e 's/ /" -> "/g' -e 's/^[^\t]*/subgraph "cluster_&" { label="&"/' -e 's/\t/ "/' -e 's/$$/" }/'; \
		 echo '}' } >$@

# Combine the ARCS, LINKS and abyss-scaffold graphs
%.arcs.a$a_l$l.links.abyss-scaffold.path.gv: %.arcs.dist.n$l.gv %.arcs.a$a_l$l.links.path.gv %.arcs.n$n.abyss-scaffold.path.gv
	bin/graph-union $^ >$@

# Filter scaffolds

# Select scaffolds at least 2 kbp.
%.2kbp.fa: %.fa
	seqtk seq -L2000 $< >$@

# Select scaffolds at least 2500 bp.
%.2500bp.fa: %.fa
	seqtk seq -L2500 $< >$@

# Select scaffolds at least 5 kbp.
%.5kbp.fa: %.fa
	seqtk seq -L5000 $< >$@

# Identify plastid scaffolds.
%/psitchensis-scaffolds.fa.cp: %/$(ref).psitchensis-scaffolds.primary.bam.names.tsv
	awk 'NR > 1 && $$2 == "$(psitchensiscp)" {print $$1}' $< | sort -n >$@

# Identify non-plastid scaffolds.
%.fa.nocp: %.fa.fai %.fa.cp
	cut -f1 $< | grep -vwFf $*.fa.cp >$@

# Remove plastid scaffolds.
%.nocp.fa: %.fa %.fa.nocp
	seqtk subseq $< $<.nocp >$@

# Create signatures of scaffolds ends

# Extract vertex names
%.gv.ends.name: %.gv
	gvpr 'N { print(name) }' $< | sort >$@

# Identify blunt vertices.
%.gv.blunt.name: %.gv
	gvpr 'N[outdegree == 0] { print(name) }' $< | sort >$@

# Create a BED file of the ends of blunt scaffolds.
%.bed: %.name
	mlr --nidx --fs tab \
		then put '$$Orientation = substr($$1, -1, -1); $$1 = substr($$1, 0, -2)' \
		then rename 2,Orientation \
		then join -u -f $(draft).fa.fai -j 1 \
		then rename 1,Name,2,Size \
		then cut -f Name,Orientation,Size \
		then put 'if ($$Orientation == "-") { $$Start = 1; $$End = $e } else { $$Start = $$Size - $e; $$End = $$Size }' \
		then cut -o -f Name,Start,End $< >$@

# Extact the alignments from the ends of scaffolds.
%.bed.bx.as100.bam: %.bed $(draft).$(name).bx.as100.bam
	samtools view -q1 -F0x904 -L $< -Obam -o $@ $(draft).$(name).bx.as100.bam

# Extact the alignments from the ends of blunt scaffolds.
%.blunt.as100.bam: %.blunt.bed $(draft).$(name).bx.as100.bam
	samtools view -q1 -F0x904 -L $< -Obam -o $@ $(draft).$(name).bx.as100.bam

# Extract the barcodes from a BAM file and store it in a TSV file.
read_length = 150
%.bam.barcodes.tsv: %.bam
	samtools view $< | gsed -r 's/\t(..):.:([^\t]*)/\t\1=\2/g' \
	| mlr --idkvp --ifs tab --otsvlite cat \
		then rename 3,Rname,4,Pos \
		then cut -f Rname,Pos,AS,BX \
		then put '$$Name = $$Rname . ($$Pos < $e + $(read_length) ? "-" : "+")' \
		then count-distinct -o Reads -f Name,BX >$@

# Count the number of barcodes per scaffold end.
%.bam.barcodes.scaffend.tsv: %.bam.barcodes.tsv
	mlr --tsvlite \
		then filter '$$Reads >= 4' \
		then stats1 -a count,sum -g Name -f Reads \
		then rename Reads_count,Barcodes \
		then put '$$Reads_per_barcode = $$Reads_sum / $$Barcodes' \
		$< >$@

# Write the list of barcodes one file per scaffold end.
%.bam.barcodes/all.tsv: %.bam.barcodes.tsv
	mkdir -p $(@D)
	mlr --itsvlite --onidx \
		then filter '$$Reads >= 4' \
		then put 'print >"$(@D)/" . $$Name . ".bx.bed", $$BX . "\t0\t1" ' $< >$@

# Write the list of barcodes one file per scaffold end.
%.bam.barcodes.scaffend/bed.files: %.bam.barcodes.tsv %.bam.barcodes.scaffend.tsv
	mkdir -p $(@D)
	mlr --itsvlite --onidx \
		then join -u -j Name -f $*.bam.barcodes.scaffend.tsv \
		then filter '$$Barcodes < 500 && $$Reads >= 4' \
		then put 'print >"$(@D)/" . $$Name . ".bx.bed", $$BX . "\t0\t1" ' \
		then uniq -f Name \
		then put '$$File = "$(@D)/" . $$Name . ".bx.bed"' \
		then cut -f File \
		$< >$@

# Convert a list of BED files to FASTQ.gz files.
%/fq.gz.files: %/bed.files
	sed 's/.bed$$/.fq.gz/' $< >$@
	make `<$@`

# Extract the reads of a single barcode to a FASTQ file.
$(name).longranger.basic.bx.%.fq.gz: $(name).longranger.basic.bam
	samtools view -@$t -h $< $* | samtools fastq -@$t - | sed '/^@.*\/[12]$$/s/$$/ BX:Z:$*/' | $(gzip) >$@

# Convert a list of barcodes in TSV format to BED format.
%.bam.barcodes.bx.bed: %.bam.barcodes.tsv
	mlr --tsvlite --headerless-csv-output \
		filter '$$Reads >=4' then cut -f BX then put '$$2 = 0; $$3 = 1' $< >$@

# Extract the reads of a set of barcodes to a FASTQ file.
%.bx.fq.gz: %.bx.bed $(name).longranger.basic.bam
	samtools view -@$t -h -L $< $(name).longranger.basic.bam | gsed 's/\tBX:Z:/\tBC:Z:/' | samtools fastq -@$t -t - | gsed 's/BC:Z:/BX:Z:/' | $(gzip) >$@

# Mash

# Construct a sketch of a FASTA file.
%.fa.msh: %.fa
	mash sketch -p $t -k 32 -s 10000 -i $<

# Construct a sketch of a FASTQ.gz file.
%.fq.gz.msh: %.fq.gz
	mash sketch -p $t -k 32 -s 10000 -m 2 $<

# Construct a sketch of a list of FASTQ.gz files.
%.fq.gz.msh: %/fq.gz.files
	mash sketch -p $t -k 32 -s 10000 -m 2 -l $< -o $*.fq.gz

# Compare all pairs of sketches.
%.msh.dist.tsv: %.msh
	mash dist $< $< >$@

# Compare sketches to a reference.
$(draft).%.msh.$(draft).fa.msh.dist.tsv: $(draft).%.msh $(draft).fa.msh
	mash dist $^ >$@

# Convert a Mash distance TSV file to a GraphViz undirected graph
%.msh.dist.tsv.u.gv: %.msh.dist.tsv
	gsed -r 's~[^\t/]*/~~;s~[^\t/]*/~~;s/\.bx\.fq\.gz//g;s~/10000$$~~;/\t10000$$/d' $< \
		| awk '$$5 >= 100' \
		| awk 'BEGIN { print "strict graph g {" } { print "\"" $$1 "\" -- \"" $$2 "\" [label=" $$5 "]" } END { print "}" }' \
		>$@

# Convert a Mash distance TSV file to a GraphViz file
%.msh.dist.tsv.gv: %.msh.dist.tsv
	gsed -r -e 's~[^\t/]*/~~;s~[^\t/]*/~~;s/\.bx\.fq\.gz//g;s~/10000$$~~;/\t10000$$/d' \
			-e '/^[^\t]*-\t/s/-\t/=\t/' \
			-e '/^[^\t]*\+\t/s/\+\t/-\t/' \
			-e '/^[^\t]*=\t/s/=\t/+\t/' $< \
		| awk '$$5 >= 100' \
		| awk 'BEGIN { print "strict digraph g {" } { print "\"" $$1 "\" -> \"" $$2 "\" [label=" $$5 "]" } END { print "}" }' \
		>$@

# Miniasm

# Align sequences using minimap
%.fa.paf: %.fa
	minimap -S $< $< >$@

# Assemble sequences using miniasm
%.miniasm.gfa: %.fa.paf %.fa
	( \
		seqtk seq -C psitchensiscpmt_10.fa | abyss-todot --gfa; \
		/home/sjackman/src/miniasm/miniasm -p sg -12 -e0 -n0 -f $*.fa $< ) >$@

# Unicycler

# Separate the first read from an interleaved FASTQ file.
%.1.fq.gz: %.fq.gz
	seqtk seq -1 $< | $(gzip) >$@

# Separate the second read from an interleaved FASTQ file.
%.2.fq.gz: %.fq.gz
	seqtk seq -2 $< | $(gzip) >$@

# Assemble reads using Unicycler.
%.unicycler/assembly.fasta: %.1.fq.gz %.2.fq.gz
	unicycler -t$t -o $(@D) -1 $*.1.fq.gz -2 $*.2.fq.gz

# Symlink the assembly.
%.unicycler.fa: %.unicycler/assembly.fasta
	seqtk seq $< >$@

# Reassemble targeted assemblies.

psitchensiscpmt_10.fa: \
		psitchensiscpmt_8/8001_0_1000000.bed.bx.as100.bam.barcodes.bx.unicycler.fa \
		psitchensiscpmt_8/8001_1000000_2000000.bed.bx.as100.bam.barcodes.bx.unicycler.polish9.fa \
		psitchensiscpmt_8/8001_2000000_3000000.bed.bx.as100.bam.barcodes.bx.unicycler.fa \
		psitchensiscpmt_8/8001_3000000_4000000.bed.bx.as100.bam.barcodes.bx.unicycler.fa \
		psitchensiscpmt_8/8001_4000000_4413118.bed.bx.as100.bam.barcodes.bx.unicycler.fa \
		psitchensiscpmt_8/8002.bed.bx.as100.bam.barcodes.bx.unicycler.fa \
		psitchensiscpmt_8/8003.bed.bx.as100.bam.barcodes.bx.unicycler.fa
	( cat KU215903.fa; \
		seqmagick convert --name-prefix 8001_0_ $(wordlist 1, 1, $^) -; \
		seqmagick convert --name-prefix 8001_1_ $(wordlist 2, 2, $^) -; \
		seqmagick convert --name-prefix 8001_2_ $(wordlist 3, 3, $^) -; \
		seqmagick convert --name-prefix 8001_3_ $(wordlist 4, 4, $^) -; \
		seqmagick convert --name-prefix 8001_4_ $(wordlist 5, 5, $^) -; \
		seqmagick convert --name-prefix 8002_ $(wordlist 6, 6, $^) -; \
		seqmagick convert --name-prefix 8003_ $(wordlist 7, 7, $^) - ) \
	| seqtk seq -L10000 >$@

# Bandage

# Render a GFA file to PNG using Bandage.
%.gfa.png: %.gfa
	Bandage image $< $@

# Render a GFA file to SVG using Bandage.
%.gfa.svg: %.gfa
	Bandage image $< $@

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

# BLAST

# Align sequences to the nt database using blastn
%.nt.blastn: %.fa
	blastn -num_threads $t -db nt -query $< -out $@

# Align sequences to the nt database using blastn and report TSV
%.nt.blastn.tsv: %.fa
	(printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tsskingdom\tssciname\tstitle\n"; \
		blastn -num_threads $t -db nt -outfmt '7 std staxid sskingdom ssciname stitle' -query $<) >$@

# Align sequences to the nt database using blastn and report TSV of the best alignment
%.nt.blastn.best.tsv: %.fa
	(printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tsskingdom\tssciname\tstitle\n"; \
		blastn -num_threads $t -db nt -outfmt '7 std staxid sskingdom ssciname stitle' -num_alignments 1 -query $<) >$@

# Align sequences to the nt database using blastn and report SAM
%.nt.blastn.sam: %.fa
	blastn -num_threads $t -db nt -outfmt 17 -query $< >$@

# Align sequences to the nt database using blastn and output an organism report
%.nt.blastn.organism: %.fa
	blastn -num_threads $t -db nt -outfmt 18 -query $< -out $@

# Align sequences to the nt database using blastx and report TSV of the best alignment
%.nr.blastx.best.tsv: %.fa
	(printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxid\tsskingdom\tssciname\tstitle\n"; \
		blastx -num_threads $t -db nr -outfmt '7 std staxid sskingdom ssciname stitle' -num_alignments 1 -query $<) >$@

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
