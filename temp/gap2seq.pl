use strict;
	
	my $scaffold_name = "genome.scf.fasta";
	my $output_name = "staph0";
	my $left_read = "frag_1.fastq";
	my $right_read = "frag_2.fastq";
	
	#command line for dataset 1 of S.aureus
	my @temp;
	@temp = (./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read");
	`@temp`;
	
	#command line for dataset 2 of S.aureus
	$output_name = "staph1";
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	@temp = ("./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read");
	#`@temp`;
	
	
	#command line for dataset 1 and 2 of S.aureus
	$output_name = "staph01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	my $left_read1 = "short_1.fastq";
	my $right_read1 = "short_2.fastq";
	@temp = ("./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read,$left_read1,$right_read1");
	`@temp`;
	
	
	#command line for dataset 3 of R.sphaeroides
	$genome_name = "rhod";
	$output_name = "rhod0";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	@temp = ("./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read");
	#`@temp`;
	
	
	#command line for dataset 4 of R.sphaeroides
	$genome_name = "rhod";
	$output_name = "rhod1";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	@temp = ("./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read");
	#`@temp`;
	
	
	#command line for dataset 3 and 4 of R.sphaeroides
	$output_name = "rhod01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	my $left_read1 = "shortjump_1.fastq";
	my $right_read1 = "shortjump_2.fastq";
	@temp = ("./Gap2Seq.sh -scaffolds $scaffold_name -filled $output_name -reads $left_read,$right_read,$left_read1,$right_read1");
	`@temp`;
	
	
	
