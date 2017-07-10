use strict;
	if(@ARGV != 3){
		die "\nPlease input correct command line:  
perl GapReduce.pl draft_sequence.fasta library.txt output_address
<draft_sequence.fasta>:
	The draft sequences of a genome.
<library.txt>:
	Each line represents one read library.
	The 1-th column: 
		the first mate read file (*.fastq);
	The 2-th column: 
		the second mate read file (*.fastq);
	The 3-th column: 
		length of read;
	The 4-th column: 
		insert size of read library;
	The 5-th column: 
		standard deviation of insert size;
	The 6-th column:  
		1 denotes paired-end reads, 0 denotes mate-paired reads;
	The 7-th column: 
		the value of the large k-mer length which should be shorter 
		than read length;
	The 8-th column: 
		the value of the small k-mer length which should be shorter than 
		read length and the integer of the 7-th column;  
	The 9-th column: 
		an integer represents the step from large k-mer length to 
		small k-mer length;
	The 10-th column: 
		an integer represents the largest k-mer frequency threshold for constructing
		de bruijn graph;
	The 11-th column: 
		denotes which mapping tool will be used, it equals 
		bwa or bowtie2;
<output_address>:
		The directory address of output files;
";
	}
	
	my $library_number = 0;
	my @library_name;
	my @library_insertsize;
	my @library_readLength;
	my @library_std;

	my @library_orientation;
	my @library_sam;
	my @library_bam;
	
	my @temp;
	
	my $scaffold_set = shift;
	my $library_information = shift;
	my $output_address = shift;
	
	if(-e $output_address){
		unlink glob "$output_address/*";
		unlink glob "$output_address/*";
	}else{
		mkdir($output_address);
	}
	my $min_gap_distance = 0;
	my @large_kmer_length;
	my @kmer_step;
	my @large_fre;
	my @kmer_length;
	my @mapping_tool;
	my $line;
	my $infor_output = "./$output_address/infor.txt";
	my $max_insert_size = 0;
	my $min_insert_size = 0;
	my $end_contig_length = 0;
	my $iterative_count = 1;
	my $real_iterative_count = 0;
	
	open(OUT, ">$infor_output")or die $!;
	open(INPUT,$library_information)or die $!;
	
	print OUT scalar localtime;
	print OUT "\n";
	while($line=<INPUT>){
	    chomp $line;
		my @infor = split /\s+/, $line;
		$library_name[$library_number*2]=$infor[0];
		$library_name[$library_number*2+1]=$infor[1];
		$library_readLength[$library_number] = $infor[2];
		$library_insertsize[$library_number]=$infor[3];
		$library_std[$library_number]=$infor[4];
		$library_orientation[$library_number]=$infor[5];
		$large_kmer_length[$library_number]=$infor[6];
		$kmer_length[$library_number]=$infor[7];
		$kmer_step[$library_number]=$infor[8];
		$large_fre[$library_number]=$infor[9];
		$mapping_tool[$library_number]=$infor[10];
		$library_number++;
	}
	
	my $folder;
	my $contig_set = "./$output_address/contig_set.fa";
	my $end_contig_set = "./$output_address/end_contig_set.fa";
	my $interval_length = 500;
	my $gap_information = "./$output_address/gap_information.fa";
	
	my $scaffold_end_contig_index = "./$output_address/scaffold_end_contig_index.fa";
	my $scaffold_set_fill_gap = "./$output_address/scaffold_set_fill_gap.fa";
	
	
	
	for(my $k=0;$k<$library_number;$k++){
		$real_iterative_count = 0;
		$iterative_count = 1;
		my $previous_all_gap_length = 0;
		while($iterative_count > 0){
			$end_contig_length = $library_insertsize[$k] + 3*$library_std[$k];
			$end_contig_length = 2*$end_contig_length;
			
			if($k != 0 || $real_iterative_count != 0){
				$scaffold_set = $scaffold_set_fill_gap;
			}
	
			print OUT "./splitScaffoldSet $scaffold_set $min_gap_distance $gap_information $contig_set $end_contig_length $interval_length $end_contig_set $scaffold_end_contig_index\n";
			@temp = ("./splitScaffoldSet $scaffold_set $min_gap_distance $gap_information $contig_set $end_contig_length $interval_length $end_contig_set $scaffold_end_contig_index");
			`@temp`;
	
			my $current_iterative_count;
			my $all_gap_length = 0;
			my $gap_length_max_insertsize = 0;
			my $gap_count = 0;
			my $max_gap_distance = 0;
			my @gap_distance;
			my @gap_left_contig_index;
			my @gap_right_contig_index;
			my $t = 0;
			open(INPUT_GAP,$gap_information)or die $!;
			while($line=<INPUT_GAP>){
				if($t==0){
					$gap_count = $line;
				}else{
					chomp $line;
					my @infor = split /\s+/, $line;
					$gap_distance[$t-1] = $infor[0];
					$gap_left_contig_index[$t-1] = $infor[1];
					$gap_right_contig_index[$t-1] = $infor[2];
					$all_gap_length = $all_gap_length + $gap_distance[$t-1];
					if($max_gap_distance < $gap_distance[$t-1]){
						$max_gap_distance = $gap_distance[$t-1];
					}
				}
				$t++;
			}
			if($real_iterative_count == 0){
				$iterative_count = $max_gap_distance/(2*($library_insertsize[$k] - $library_readLength[$k]));
				if($iterative_count > 3){
					$iterative_count = 3;
				}
			}else{
				if($previous_all_gap_length - $all_gap_length < 200){
					$iterative_count = 1;
				}
			}
			$previous_all_gap_length = $all_gap_length;
			
			$library_sam[2*$k] = "./$output_address/library_$k"."_left.sam";
			$library_sam[2*$k+1] = "./$output_address/library_$k"."_right.sam";
			$library_bam[2*$k] = "./$output_address/library_$k"."_left.bam";
			$library_bam[2*$k+1] = "./$output_address/library_$k"."_right.bam";
			if($mapping_tool[$k] eq "bowtie2"){
				@temp = ("bowtie2-build $end_contig_set ./$output_address/contigs");
				system(@temp) == 0 or die "\nThe command 'bowtie2-build' can not be found! Please install the mapping tool BWA\n";;
				@temp = ("bowtie2 -x ./$output_address/contigs $library_name[2*$k] -S $library_sam[2*$k] 2>&-");
				`@temp`;
				@temp = ("bowtie2 -x ./$output_address/contigs $library_name[2*$k+1] -S $library_sam[2*$k+1] 2>&-");
				`@temp`;
			}
			if($mapping_tool[$k] eq "bwa"){
				@temp = ("bwa index $end_contig_set");
				
				system(@temp) == 0 or die "\nThe command 'bwa' can not be found! Please install the mapping tool BWA\n";
				if($library_readLength[$k] > 33){
					@temp = ("bwa mem $end_contig_set $library_name[2*$k] > $library_sam[2*$k] 2>&-");
					`@temp`;
					@temp = ("bwa mem $end_contig_set $library_name[2*$k+1] > $library_sam[2*$k+1] 2>&-");
					`@temp`;
				}else{
					@temp = ("bwa aln $end_contig_set $library_name[2*$k] > ./$output_address/reads.sai 2>&-");
					`@temp`;
					@temp = ("bwa samse $end_contig_set ./$output_address/reads.sai $library_name[2*$k] > $library_sam[2*$k] 2>&-");
					`@temp`;
					
					@temp = ("bwa aln $end_contig_set $library_name[2*$k+1] > ./$output_address/reads.sai 2>&-");
					`@temp`;
					@temp = ("bwa samse $end_contig_set ./$output_address/reads.sai $library_name[2*$k+1] > $library_sam[2*$k+1] 2>&-");
					`@temp`;
				}
				
			}
			@temp = ("samtools view -Sb $library_sam[2*$k] > $library_bam[2*$k]");
			system(@temp) == 0 or die "\nThe command 'samtools' can not be found! Please install the tool Samtools\n";
			@temp = ("samtools view -Sb $library_sam[2*$k+1] > $library_bam[2*$k+1]");
			`@temp`;
			unlink glob $end_contig_set.".*";
			unlink glob $library_sam[2*$k];
			unlink glob $library_sam[2*$k+1];
			
			unlink glob "./$output_address/contigs.*";
			
			$max_insert_size = $library_insertsize[$k] + 3*$library_std[$k];
			$min_insert_size = $library_insertsize[$k] - 3*$library_std[$k];
			$folder = "./$output_address/gap_read_out_put_$k";
			if(-e $folder){
				unlink glob "$folder/*";
				unlink glob "$folder/*";
			}else{
				mkdir($folder);
			}
			print OUT "./fillGapRead $scaffold_set $contig_set $min_gap_distance $library_bam[2*$k] $library_bam[2*$k+1] $max_insert_size $min_insert_size $interval_length $library_orientation[$k] $k $output_address\n";
			@temp = ("./fillGapRead $scaffold_set $contig_set $min_gap_distance $library_bam[2*$k] $library_bam[2*$k+1] $max_insert_size $min_insert_size $interval_length $library_orientation[$k] $k $output_address");
			`@temp`;
			
			unlink $library_bam[2*$k], $library_bam[2*$k+1];
			
			my $gap_region = "./$output_address/fill_gap_region.fa";
			my $fill_gap_region = "./$output_address/fill_gap_region.fa";
			my $fill_parameter;
			
			if(-e $fill_gap_region){
				unlink $fill_gap_region;
			}
			print OUT "gap_cout:$gap_count\n";
			for($t=0;$t<$gap_count;$t++){
				$fill_parameter = "./fillGapRegion";
				my $left_gap_read = "leftReadSet_gapIndex_$t.fa";
				my $right_gap_read = "rightReadSet_gapIndex_$t.fa";
				
				$folder = "./$output_address/gap_read_out_put_$k/";
				$left_gap_read = $folder.$left_gap_read;
				$right_gap_read = $folder.$right_gap_read;
				$fill_parameter = $fill_parameter.' '.$left_gap_read.' '.$right_gap_read.' '.$library_insertsize[$k].' '.$library_std[$k].' '.$library_orientation[$k];
				
				$fill_parameter = $fill_parameter.' '.$large_kmer_length[$k].' '.$kmer_length[$k].' '.$gap_left_contig_index[$t].' '.$gap_right_contig_index[$t].' '.$gap_distance[$t];
				print OUT "$fill_parameter $fill_gap_region $contig_set\n";
				@temp = ("$fill_parameter $fill_gap_region $contig_set $kmer_step[$k] $large_fre[$k] >> ./$output_address/out.fa");
				`@temp`;
			}
			
			@temp = ("./GetFinalScaffoldFill $contig_set $fill_gap_region $scaffold_end_contig_index $scaffold_set_fill_gap");
			`@temp`;
			print OUT "real_iterative_count:$real_iterative_count\n";
			$real_iterative_count++;
			$iterative_count--;
		}
	}
	#unlink "graph.fa", "fill_gap_region.fa", "scaffold_end_contig_index.fa", "gap_information.fa", "end_contig_set.fa", "contig_set.fa";
	
	
	
	
	
	
	
	
	
