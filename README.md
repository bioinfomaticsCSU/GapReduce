License
=========

Copyright (C) 2014 Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Junwei Luo(luojunwei@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


GapReduce
=================
1) Introduction

	GapReduce is a gap filling tool based on paired reads.
	The input includes the draft sequences of a genome (scaffolds or contigs which have gaps), and one or multiple paired read sets. 

2) Before installing and running
	
	Users should install Bowtie2, BWA and Samtools firstly and add them to your PATH. GapReduce uses Bowtie2 or BWA map read to contigs in the step of scaffolding, and convert ".sam" file to ".bam" file. 
	Users can download BWA from https://github.com/lh3/bwa
	Users can download Bowtie2 from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml 
	Samtools is available from http://samtools.sourceforge.net/index.shtml
	Please type "bwa" to check whether BWA works.
	Please type "bowtie2" to check whether Bowtie2 works.
	Please type "samtools" to check whether samtool works.

3) Installing

	GapReduce is written C++ and therefore will require a machine with GNU C++ pre-installed.
	Create a main directory (eg:GapReduce). Copy all source code to this directory.
	Type "make all".

4) Running

	Run command line:  
	perl GapReduce.pl library.txt draft_sequence.fasta
	<library.txt>:
		Each line represents one read library.
		The first column is the first mate read file (*.fastq);
		Tthe sencond column is the second mate read file (*.fastq);
		The third column is length of read;
		The fourth column is insert size of read library;
		The fifth column is standard deviation of insert size;
		The sixth column represents whether the read library is paired-end reads (1 denotes paired-end reads, 0 denotes mate-paired reads);
		The seventh column denotes which mapping tool will be used, it equals bwa or bowtie2;
	<draft_sequence.fasta>:
		The draft sequences of a genome which includes some gaps.

4) Output:

	The final filling result is named "scaffold_set_fill_gap.fa".

5) Example:

	one line in library.txt:
	frag_1.fastq frag_2.fastq 101 180 20 0 bwa

Evaluation
=================
1) Before evaluation

	Users should install MUMmer and add them to your PATH, which is used for aligning contigs to reference genomes.
	MUMmer is available from http://mummer.sourceforge.net/;
	Please type "nucmer" to check whether MUMmer works.
	Please enter directory ./GapReduce/evaluation

2) Simulating the draft sequences of a genome by contigs produced by a assembler

	This module is a stand-alone module of GapReduce, which is available after installing GapReduce.
	First, it aligns contigs against the references of a genome by MUMmer. The regions of the references covered by contigs are kept. Any position of the references not covered by any contig is replaced by symbol ‘N’, which will be identified as gaps. Then, it outputs the changed references as the draft sequences of the genome, and the real sequences corresponding to gaps. 
	Run command line:
		perl perfect_scaffold.pl contig_file reference_file perfect_scaffold_file gap_file
		<contig_file>:
			The file which contains the contigs produced by a assembler;
		<reference_file>:
			This is the reference genome file;
		<perfect_scaffold_file>:
			The output name of the draft sequences of the genome;
		<gap_file>:
			The output name of the real sequence of gap regions in the draft sequences;

3) Evaluating the performance of gap filling tools

	This module is a stand-alone module of GapReduce, which is avaiable after installing GapReduce.
	It first breaks the draft sequences of a genome at the positions of the gaps, and produces contigs containing no ‘N’. Then, it performs MUMmer by using the contigs and the draft sequences filled by a gap filling tool as queries and targets, respectively. Based on the alignment position of the queries to targets, it extracts the sequence between two adjacent contigs as filled results of corresponding gap. It adopts the global algorithm Needleman-Wunsch to align the extracted gap sequences to the real gap sequences.
	Run command line:
		perl evaluate_gap.pl draft_sequences.fasta gap_file.fasta draft_sequences_filled.fasta
		<draft_sequences.fasta>:
			The original draft sequences of a genome before gap filling;
		<gap_file.fasta>:
			The real sequences of gaps in the reference genome;
		<draft_sequences_filled.fasta>:
			The gap filling result by a gap filling tool;
	The evaluation results is stored in a file named "fscore.fa". 
	In this file, it includes gap_count, gap_length, match_count, mismatch_count precision, recall, and fscore.
		gap_count is the count of gaps in the draft sequences of genomes;
		gap_length is the length of all gaps;
		match_count is the count of positions in gaps filled correctly;
		mismatch_count is the count of positions in gaps filled incorrectly;
		precision is the ratio of match_count to match_count and mismatch_count; precision = match_count/(match_count + mismatch_count);
		recall is the ratio of match_count to gap_length; recall = match_count/gap_length;
		fscore is a composite metric to evaluate the gap filling results; fscore = 2*precision*recall/(precision + recall);
