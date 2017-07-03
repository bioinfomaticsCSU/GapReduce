use strict;
	

	my $scaffold_name = "genome.scf.fasta";
	my $output_name = "$scaffold_set\_staph0";
	my $left_read = "frag_1.fastq";
	my $right_read = "frag_2.fastq";
	my $read_length = 101;
	my $insert_size = 180;
	my $std = 18;
	my $paired_library = 1;
	my $max_kmer_length = 61;
	my $min_kmer_length = 31;
	my $step = 15;
	my $fre = 2;
	my $mapping_tool = "bwa";
	
	
	#command line for dataset 1 of S.aureus
	open Output, ">library0.txt";
	my $line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $step $fre $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library0.txt $output_name");
	system(@temp);

	
	#command line for dataset 2 of S.aureus
	$output_name = "staph1";
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$read_length = 37;
	$insert_size = 3500;
	$std = 350;
	$paired_library = 0;
	$max_kmer_length = 29;
	$min_kmer_length = 13;
	open Output, ">library1.txt";
	my $line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library1.txt $output_name");
	system(@temp);
	
	
	
	#command line for dataset 1 and 2 of S.aureus
	$output_name = "staph01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$read_length = 101;
	$insert_size = 180;
	$std = 18;
	$paired_library = 1;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	$mapping_tool = "bwa";
	open Output, ">library01.txt";
	$line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool\n";
	print Output $line;
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$read_length = 37;
	$insert_size = 3500;
	$std = 350;
	$paired_library = 0;
	$max_kmer_length = 29;
	$min_kmer_length = 13;
	$line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library01.txt $output_name");
	system(@temp);
	
	
	#command line for dataset 3 of R.sphaeroides
	$genome_name = "rhod";
	$output_name = "rhod0";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$read_length = 101;
	$insert_size = 180;
	$std = 18;
	$paired_library = 1;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	open Output, ">library2.txt";
	my $line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library2.txt $output_name");
	system(@temp);
	
	
	#command line for dataset 4 of R.sphaeroides
	$genome_name = "rhod";
	$output_name = "rhod1";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$read_length = 101;
	$insert_size = 3500;
	$std = 350;
	$paired_library = 0;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	open Output, ">library3.txt";
	my $line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library3.txt $output_name");
	system(@temp);
	
	
	#command line for dataset 3 and 4 of R.sphaeroides
	$genome_name = "rhod";
	$output_name = "rhod01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$read_length = 101;
	$insert_size = 180;
	$std = 18;
	$paired_library = 1;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	open Output, ">library23.txt";
	$line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool\n";
	print Output $line;
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$read_length = 101;
	$insert_size = 3500;
	$std = 350;
	$paired_library = 0;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	$line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library23.txt $output_name");
	system(@temp);
	
	
	
	#command line for dataset 5 of Human 14
	$genome_name = "human";
	$output_name = "human0";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$read_length = 101;
	$insert_size = 2500;
	$std = 250;
	$paired_library = 0;
	$max_kmer_length = 61;
	$min_kmer_length = 31;
	open Output, ">library4.txt";
	my $line = "$left_read $right_read $read_length $insert_size $std $paired_library $max_kmer_length $min_kmer_length $mapping_tool";
	print Output $line;
	my @temp;
	@temp = ("perl GapReduce.pl $scaffold_name library4.txt $output_name");
	system(@temp);

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
