# Manual of SOAPdenovo-Trans

## Introduction

SOAPdenovo-Trans is a de novo transcriptome assembler basing on the SOAPdenovo framework, adapt to alternative splicing and different expression level among transcripts.The assembler provides a more accurate, complete and faster way to construct the full-length transcript sets.

## System Requirement

SOAPdenovo-Trans aims for the transcript assembly. It runs on 64-bit Linux systems. For animal transcriptomes like mouse, about 30-35GB memory would be required.

## Update Log
1.04 | 2014-04-22 15:00:00 +0800 (Tue, 22 Apr 2014)
Fixes a number of 'seqmentation fault' errors on different kinds of data.
(Thanks for Chris Boursnell (twitter: @chrisboursnell) fixing the bugs.)

1.03 | 2013-07-19 12:00:00 +0800 (Fri, 19 Jul 2013)
Add the function: calculate RPKM (Reads per Kilobase of assembled transcripts per Million mapped reads).

## Installation

1. You can download the pre-compiled binary according to your platform, unpack and execute directly.
2. Or download the source code, unpack to ${destination folder} with the method above, and compile by using GNU make with command "sh make.sh" at ${destination folder} and generate the executable files "SOAPdenovo-Trans-31mer" and "SOAPdenovo-Trans-127mer".

## How to use it

### 1. Configuration file
The configuration file in SOAPdenovo-Trans is mostly the same as SOAPdenovo, but there is no "rank" parameter. The configuration file tells the assembler where to find these files and the relevant information. "example.config" demonstrates how to organize the information and make configuration file.

The configuration file has a section for global information, and then multiple library sections. Right now only "max_rd_len" is included in the global information section. Any read longer than max_rd_len will be cut to this length. The library information and the information of sequencing data generated from the library should be organized in the corresponding library section. Each library section starts with tag [LIB] and includes the following items:
<pre>
1) avg_ins
This value indicates the average insert size of this library or the peak value position in the insert size distribution figure.
2) reverse_seq
This option takes value 0 or 1. It tells the assembler if the read sequences need to be complementarily reversed. Illumima GA produces two types of paired-end libraries: a) forward-reverse, generated from fragmented DNA ends with typical insert size less than 500 bp; b) reverse-forward, generated from circularizing libraries with typical insert size greater than 2 Kb. The parameter "reverse_seq" should be set to indicate this: 0, forward-reverse; 1, reverse-forward.
3) asm_flags
This indicator decides in which part(s) the reads are used. It takes value 1(only contig assembly), 2 (only scaffold assembly), 3(both contig and scaffold assembly).
4) rd_len_cutof
The assembler will cut the reads from the current library to this length.
5) map_len
This takes effect in the "map" step and is the mininum alignment length between a read and a contig required for a reliable read location. The minimum length for paired-end reads and mate-pair reads is 32 and 35 respectively. The assembler accepts read file in three kinds of formats: FASTA, FASTQ and BAM. Mate-pair relationship could be indicated in two ways: two sequence files with reads in the same order belonging to a pair, or two adjacent reads in a single file (FASTA only) are belonging to a pair. In the configuration file single end files are indicated by "f=/path/filename" or "q=/path/filename" for fasta or fastq formats separately. Paired reads in two fasta sequence files are indicated by "f1=" and "f2=". While paired reads in two fastq sequences files are indicated by "q1=" and "q2=". Paired reads in a single fasta sequence file is indicated by "p=" item. Reads in bam sequence files is indicated by "b=". All the above items in each library section are optional. The assembler assigns default values for most of them. If you are not sure how to set a parameter, you can remove it from your configuration file.
</pre>

### 2. Get it started
Once the configuration file is available, the simplest way to run the assembler is:
<pre>
./SOAPdenovo-Trans all -s config_file -o output_prefix
</pre>

User can also choose to run the assembly process step by step as:
<pre>
./SOAPdenovo-Trans pregraph -s config_file -K 31 -o outputGraph 
./SOAPdenovo-Trans contig -g outputGraph
./SOAPdenovo-Trans map -s config_file -g outputGraph
./SOAPdenovo-Trans scaff -g outputGraph -F
</pre>

NOTE: SOAPdenovo-Trans has two versions: SOAPdenovo-Trans-31mer and SOAPdenovo-Trans-127mer.

### 3. Options:
<pre>
SOAPdenovo-Trans all -s configFile -o outputGraph [-R -f -S -F] [-K kmer -p n_cpu -d kmerFreqCutoff -e EdgeCovCutoff -M mergeLevel -L minContigLen -t locusMaxOutput -G gapLenDiff] 

-s	<string>		configFile: the config file of reads
-o	<string>		outputGraph: prefix of output graph file name
-g	<string>		inputGraph: prefix of input graph file names
-R	(optional)		output assembly RPKM statistics, [NO]
-f	(optional)		output gap related reads for SRkgf to fill gap, [NO]
-S	(optional)		scaffold structure exists, [NO]
-F	(optional)		fill gaps in scaffolds, [NO]
-K	<int>		kmer (min 13, max 31/127): kmer size, [23]
-p	<int>		n_cpu: number of cpu for use, [8]
-d	<int>		kmerFreqCutoff: kmers with frequency no larger than KmerFreqCutoff will be deleted, [0]
-e	<int>		EdgeCovCutoff: edges with coverage no larger than EdgeCovCutoff will be deleted, [2]
-M	<int>		mergeLevel (min 0, max 3): the strength of merging similar sequences during contiging, [1]
-L	<int>		minContigLen: shortest contig for scaffolding, [100]
-t	<int>		locusMaxOutput: output the number of transcripts no more than locusMaxOutput in one locus, [5]
-G	<int>		gapLenDiff: allowed length difference between estimated and filled gap, [50]
</pre>

### 4. Output files
These files are output as assembly results:
*.contig contig sequence file
*.scafSeq scaffold sequence file
There are some other files that provide useful information for advanced users, which are listed in Appendix B.

### 5. Parameter adjustment
<pre>
-K:	
The kmer value is always depended on data size and its transcript features. At the current stage, SOAPdenovo-Trans has two versions: 
1.	SOAPdenovo-Trans-31mer, which accepts odd Kmer value from 13 to 31,
2.	SOAPdenovo-Trans-127mer, which accepts odd Kmer value from 13 to 127.
Ordinarily, SOAPdenovo-Trans always assembles the RNA-seq data by small kmer (~35-mer) as some of the transcripts are in low expression level.

-M:
The parameter is used to adjust the strength of merging similar sequences (bubble structure in the de Bruijn Graph). The similar sequences may be caused by sequencing error, repeat sequence, paralogs and heterozygosis.

-L:
The parameter is sensitive to the accuracy and coverage rate of the assemble result. The larger value we set, the higher the accuracy and lower coverage result we get. For example in the Mouse data assembly, the accuracy and coverage rate are 56.05% and 99.03% respectively when "-L" is set to 50. The accuracy and coverage rate are 88.20% and 98.22% respectively when the value is 150. The accuracy is increased by 32.15%. However, the coverage rate is decreased by only 0.79%.

-e:
The parameter is always 1~3, default 2. It deletes the low-coverage edge.

-t:
In the software, sub-graph may be a group of transcripts. The parameter is to set the upper limit of the transcripts output from each sub-graph. 5 is the experiential value. If the complexity of transcriptome is high, increase it. If the user needs higher true positive, decrease it.

-F:
The parameter is optimal. It is set to utilize the gap filling step in SOAPdenovo-Trans. The user can use other gap-filling softwares instead.

-R:
The parameter is optional. It is set to calculate RPKM (Reads per Kilobase of assembled transcripts per Million mapped reads). Note that the RPKM is calculated by unique reads base on Kmer mapping. The information could be found from the file "*.RPKM.Stat". 
Furthermore, the function will output two intermediate files: *. readInformation and *.scafStatistics. They record the information about the reads' locations on contigs and scaffolds.

-S:
The parameter is optional. It is set to skip the step of constructing the scaffold and fill gap directly. It is used to try new methods of gap filling for developers. Or, if users have assembled the data with (or without) gap filling, they can get the assemblies rapidly without (or with) gap filling on the premise of keeping all the output files after scaffolding.
</pre>

# APPENDIX A: 	example.config 

<pre>
#maximal read length
max_rd_len=50
[LIB]
#maximal read length in this lib
rd_len_cutof=45
#average insert size
avg_ins=200
#if sequence needs to be reversed 
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
#fastq file for read 1 
q1=/path/**LIBNAMEA**/fastq_read_1.fq
#fastq file for read 2 always follows fastq file for read 1
q2=/path/**LIBNAMEA**/fastq_read_2.fq
#fasta file for read 1 
f1=/path/**LIBNAMEA**/fasta_read_1.fa
#fastq file for read 2 always follows fastq file for read 1
f2=/path/**LIBNAMEA**/fasta_read_2.fa
#fastq file for single reads
q=/path/**LIBNAMEA**/fastq_read_single.fq
#fasta file for single reads
f=/path/**LIBNAMEA**/fasta_read_single.fa
#a single fasta file for paired reads
p=/path/**LIBNAMEA**/pairs_in_one_file.fa
</pre>

# APPENDIXA B:

## 1. Output files from the command "pregraph"
<pre>
a. *.kmerFreq
Each row shows the number of Kmers with a frequency equals the row number. Note that those peaks of frequencies which are the integral multiple of 63 are due to the data structure.
b. *.edge
  Each record gives the information of an edge in the pre-graph: length, Kmers on both ends, average kmer coverage, whether it's reverse-complementarily identical and the sequence.
c. *.preArc
  Connections between edges which are established by the read paths.
d. *.vertex
  Kmers at the ends of edges.
e. *.preGraphBasic
  Some basic information about the pre-graph: number of vertex, K value, number of edges, maximum read length etc.
</pre>

## 2. Output files from the command "contig"
<pre>
a. *.contig
  Contig information: corresponding edge index, length, kmer coverage, whether it's tip and the sequence. Either a contig or its reverse complementary counterpart is included. Each reverse complementary contig index is indicated in the *.ContigIndex file.
b. *.Arc
  Arcs coming out of each edge and their corresponding coverage by reads.
c. *.updated.edge
  Some information for each edge in graph: length, Kmers at both ends, index difference between the reverse-complementary edge and this one.
d. *.ContigIndex
  Each record gives information about each contig in the *.contig: it's edge index, length, the index difference between its reverse-complementary counterpart and itself.
</pre>

## 3. Output files from the command "map"
<pre>
a. *.peGrads
  Information for each clone library: insert-size, read index upper bound, rank and pair number cutoff for a reliable link. This file can be revised manually for scaffolding tuning.
b. *.readOnContig
  Reads' locations on contigs. Here contigs are referred by their edge index. However about half of them are not listed in the *.contig file for their reverse-complementary counterparts are included already.
c. *.readInGap
  This file includes reads that could be located in gaps between contigs. This information will be used to close gaps in scaffolds if "-F" is set.
d. *. readInformation
  Reads' locations on contigs: read id, start position of read, contig id, start position of contig, the align length and orientation.
</pre>

## 4. Output files from the command "scaff"
<pre>
a. *.newContigIndex
  Contigs are sorted according their length before scaffolding. Their new indexes are listed in this file.  This is useful if one wants to corresponds contigs in *.contig with those in *.links.
b. *.links
  Links between contigs which are established by read pairs. New index are used.
c. *.scaf_gap
  Contigs in gaps are found by contig graph outputted by the contiging procedure. Here new index are used.
d. *.scaf
  Contigs for each scaffold: contig index (concordant to index in *.contig), approximate start position on scaffold, orientation, contig length, and its links to others contigs.
e. *.gapSeq
  Gap sequences between contigs.
f. *.scafSeq
  Sequences of each scaffolds.
g. *.contigPosInscaff
  Contigs' positions in each scaffold.
h. *.readOnScaf
  Reads' locations on each scaffold: read id, start position of read, start position of scaffold, orientation and align length in order.
i. *.scafStatistics
  Statistic information of final scaffold and contig.
j. *.RPKM.Stat
  RPKM information: transcript ID, Transcript Length, Unique reads number, RPKM value.
</pre>
