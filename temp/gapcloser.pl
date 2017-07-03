use strict;
	
	#command line for dataset 1 of S.aureus
	my $scaffold_name = "genome.scf.fasta";
	my $output_name = "staph0";
	my $left_read = "frag_1.fastq";
	my $right_read = "frag_2.fastq";
	my $insert_size = 180;
	my $mate_library = 0;

	open Output, ">library.txt";
	my $line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = (./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);

    #command line for dataset 2 of S.aureus
	$output_name = "staph1";
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$insert_size = 3500;
	$mate_library = 1;
	open Output, "library.txt";
	my $line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 37 -o $output_name");
	system(@temp);

	#command line for dataset 1 and 2 of S.aureus
	$output_name = "staph01";
	open Output, ">library.txt";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$mate_library = 0;
	$line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$insert_size = 3500;
	$mate_library = 1;
	$line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);

    #command line for dataset 3 of R.sphaeroides
	$output_name = "rhod0";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$mate_library = "0";
	open Output, ">library.txt";
	my $line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);

	
	#command line for dataset 4 of R.sphaeroides
	$output_name = "rhod1";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 3500;
	$mate_library = 1;
	open Output, ">library.txt";
	my $line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);
	
	#command line for dataset 3 and 4 of R.sphaeroides
	$output_name = "rhod01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$mate_library = "0";
	open Output, ">library.txt";
	$line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 3500;
	$mate_library = 1;
	$line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);
	
	#command line for dataset 5 of Human 14
	$output_name = "human0";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 2500;
	$mate_library = 1;
	open Output, ">library.txt";
	my $line = "[LIB]\navg_ins=$insert_size\nreverse_seq=$mate_library\nasm_flags=3\nq1=$left_read\nq2=$right_read\n";
	print Output $line;
	my @temp;
	@temp = ("./GapCloser -b library.txt -a $scaffold_name -l 101 -o $output_name");
	system(@temp);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
