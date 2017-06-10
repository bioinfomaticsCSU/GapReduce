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
```
	GapReduce is a gap filling tool based on paired reads.
	The input includes the draft sequences of a genome (scaffolds or contigs which have gaps), and one or multiple paired read sets. 
```
2) Before installing and running
```	
	Users should install Bowtie2, BWA and Samtools firstly and add them to your PATH. GapReduce uses Bowtie2 or BWA map read to contigs in the step of scaffolding, and convert ".sam" file to ".bam" file. 
	Users can download BWA from https://github.com/lh3/bwa
	Users can download Bowtie2 from http://bowtie-bio.sourceforge.net/bowtie2/index.shtml 
	Samtools is available from http://samtools.sourceforge.net/index.shtml
	Please type "bwa" to check whether BWA works.
	Please type "bowtie2" to check whether Bowtie2 works.
	Please type "samtools" to check whether samtool works.
```
3) Installing
```
	GapReduce should run on Linux operating sysetm with gcc. We test GapReduce using gcc4.6.3 on Ubuntu.
	Create a main directory (eg:GapReduce). Copy all source code to this directory.
	Please add current path to the enviroment variable LD_LIBRARY_PATH.
	Type "make all".
```
4) Running
```
	Run command line:  
	perl GapReduce.pl draft_sequence.fasta library.txt 
	<draft_sequence.fasta>:
		The draft sequences of a genome which includes some gaps.
	<library.txt>:
		Each line represents one read library.
		The 1-th column: the first mate read file (*.fastq);
		The 2-th column: the second mate read file (*.fastq);
		The 3-th column: the length of read;
		The 4-th column: insert size of read library;
		The 5-th column: standard deviation of insert size;
		The 6-th column: 1 denotes paired-end reads, 0 denotes mate-paired reads;
		The 7-th column: the value of the large k-mer length which should be shorter than read length;
		The 8-th column: the value of the small k-mer length which should be shorter than read length and the integer of the 7-th column;  
		The 9-th column: an integer represents the step from large k-mer length to small k-mer length;
		The 10-th column: an integer represents the largest k-mer frequency threshold for constructing de bruijn graph;
		The 11-th column: denotes which mapping tool will be used, it equals 'bwa' or 'bowtie2';
```
4) Output:
```
	The final filling result is named "scaffold_set_fill_gap.fa".
```
5) Example:
```
	one line in library.txt:
	frag_1.fastq frag_2.fastq 101 180 20 0 61 31 15 2 bwa
```
