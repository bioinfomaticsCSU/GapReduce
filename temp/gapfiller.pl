use strict;
	
	my $scaffold_name = "genome.scf.fasta";
	my $output_name = "$scaffold_set\_staph0";
	my $left_read = "frag_1.fastq";
	my $right_read = "frag_2.fastq";
	my $insert_size = 180;
	my $std = 0.25;
	my $paired_library = "FR";
	my $mapping_tool = "bwa";
    
	#command line for dataset 1 of S.aureus
	open Output, ">library.txt";
	my $line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	#command line for dataset 2 of S.aureus
	$output_name = "staph1";
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$insert_size = 3500;
	$paired_library = "RF";
	open Output, ">library.txt";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	#command line for dataset 1 and 2 of S.aureus
	$output_name = "staph01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$std = 0.25;
	$paired_library = "FR";
	$mapping_tool = "bwa";
	open Output, ">library.txt";
	my $line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	$left_read = "short_1.fastq";
	$right_read = "short_2.fastq";
	$insert_size = 3500;
	$paired_library = "RF";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	#command line for dataset 3 of R.sphaeroides
	$output_name = "rhod0";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$paired_library = "FR";
	open Output, ">library.txt";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	#command line for dataset 4 of R.sphaeroides
	$output_name = "rhod1";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 3500;
	$paired_library = "RF";
	open Output, ">library.txt";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	#system(@temp);
	unlink("$assembler\_$scaffold_set\_library.txt");
	
	
	#command line for dataset 3 and 4 of R.sphaeroides
	$output_name = "rhod01";
	$left_read = "frag_1.fastq";
	$right_read = "frag_2.fastq";
	$insert_size = 180;
	$paired_library = "FR";
	open Output, ">library.txt";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 3500;
	$paired_library = "RF";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	#command line for dataset 5 of Human 14
	$output_name = "human0";
	$left_read = "shortjump_1.fastq";
	$right_read = "shortjump_2.fastq";
	$insert_size = 2500;
	$paired_library = "RF";
	open Output, ">library.txt";
	$line = "lib $mapping_tool $left_read $right_read $insert_size $std $paired_library\n";
	print Output $line;
	my @temp;
	@temp = ("perl GapFiller.pl -l library.txt -s $scaffold_name -m 30 -o 2 -r 0.7 -n 10 -d 50 -t 10 -g 0 -T 1 -i 3 -b $output_name");
	system(@temp);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
