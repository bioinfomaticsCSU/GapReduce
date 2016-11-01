use strict;
	if(@ARGV != 2){
		die "\nPlease input correct command line:  
perl GapReduce.pl draft_sequence.fasta library.txt 
<draft_sequence.fasta>:
	The draft sequences of a genome.
<library.txt>:
	Each line represents one read library.
	The first column: 
		the first mate read file (*.fastq);
	The sencond column: 
		the second mate read file (*.fastq);
	The third column: 
		length of read;
	The fourth column: 
		insert size of read library;
	The fifth column: 
		standard deviation of insert size;
	The sixth column:  
		1 denotes paired-end reads, 0 denotes mate-paired reads;
	The seventh column: 
		a integer which should be shorter than read length;
	The eighth column: 
		a integer which should be shorter than read length and the 
		integer of the seventh column; This integer is the length of 
		the k-mers which are used for building de Bruijn graph; 
	The ninth column: 
		denotes which mapping tool will be used, it equals 
		bwa or bowtie2;";
	}
	
	my $library_number = 0;
	my @library_name;
	my @library_insertsize;
	my @library_readLength;
	my @library_std;

	my @library_orientation;
	my @library_sam;
	my @library_bam;
	
	my $scaffold_set = shift;
	my $library_information = shift;
	my $min_gap_distance = 0;
	my @large_kmer_length;
	my @kmer_length;
	my @mapping_tool;
	my $line;
	my $infor_output = "infor.txt";
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
		$mapping_tool[$library_number]=$infor[8];
		$library_number++;
	}
	
	my $folder;
	my $contig_set = "contig_set.fa";
	my $end_contig_set = "end_contig_set.fa";
	my $interval_length = 500;
	my $gap_information = "gap_information.fa";
	
	my $scaffold_end_contig_index = "scaffold_end_contig_index.fa";
	my $scaffold_set_fill_gap = "scaffold_set_fill_gap.fa";
	
	my @temp;
	
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
			
			print OUT "./splitScaffoldSet $scaffold_set $min_gap_distance $gap_information $contig_set $end_contig_length $interval_length $end_contig_set\n";
			@temp = ("./splitScaffoldSet $scaffold_set $min_gap_distance $gap_information $contig_set $end_contig_length $interval_length $end_contig_set");
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
			}else{
				if($previous_all_gap_length - $all_gap_length < 200){
					$iterative_count = 1;
				}
			}
			$previous_all_gap_length = $all_gap_length;
			
			$library_sam[2*$k] = "library_$k"."_left.sam";
			$library_sam[2*$k+1] = "library_$k"."_right.sam";
			$library_bam[2*$k] = "library_$k"."_left.bam";
			$library_bam[2*$k+1] = "library_$k"."_right.bam";
			if($mapping_tool[$k] eq "bowtie2"){
				@temp = ("bowtie2-build $end_contig_set contigs");
				system(@temp) == 0 or die "\nThe command 'bowtie2-build' can not be found! Please install the mapping tool BWA\n";;
				@temp = ("bowtie2 -x contigs $library_name[2*$k] -S $library_sam[2*$k]");
				`@temp`;
				@temp = ("bowtie2 -x contigs $library_name[2*$k+1] -S $library_sam[2*$k+1]");
				`@temp`;
			}
			if($mapping_tool[$k] eq "bwa"){
				@temp = ("bwa index $end_contig_set");
				system(@temp) == 0 or die "\nThe command 'bwa' can not be found! Please install the mapping tool BWA\n";
				@temp = ("bwa mem $end_contig_set $library_name[2*$k] > $library_sam[2*$k]");
				`@temp`;
				@temp = ("bwa mem $end_contig_set $library_name[2*$k+1] > $library_sam[2*$k+1]");
				`@temp`;
			}
			
			@temp = ("samtools view -Sb $library_sam[2*$k] > $library_bam[2*$k]");
			system(@temp) == 0 or die "\nThe command 'samtools' can not be found! Please install the tool Samtools\n";
			@temp = ("samtools view -Sb $library_sam[2*$k+1] > $library_bam[2*$k+1]");
			`@temp`;
			unlink glob $end_contig_set.".*";
			unlink glob $library_sam[2*$k];
			unlink glob $library_sam[2*$k+1];
			
			unlink glob "contigs.*";
			
			$max_insert_size = $library_insertsize[$k] + 3*$library_std[$k];
			$min_insert_size = $library_insertsize[$k] - 3*$library_std[$k];
			$folder = "gap_read_out_put_$k";
			if(-e $folder){
				unlink glob "$folder/*";
				unlink glob "$folder/*";
			}else{
				mkdir($folder);
			}
			print OUT "./fillGapRead $scaffold_set $contig_set $min_gap_distance $library_bam[2*$k] $library_bam[2*$k+1] $max_insert_size $min_insert_size $interval_length $library_orientation[$k] $k\n";
			@temp = ("./fillGapRead $scaffold_set $contig_set $min_gap_distance $library_bam[2*$k] $library_bam[2*$k+1] $max_insert_size $min_insert_size $interval_length $library_orientation[$k] $k");
			`@temp`;
			
			unlink $library_bam[2*$k], $library_bam[2*$k+1];
			
			my $gap_region = "fill_gap_region.fa";
			my $fill_gap_region = "fill_gap_region.fa";
			my $fill_parameter;
			
			if(-e $fill_gap_region){
				unlink $fill_gap_region;
			}
			print OUT "gap_cout:$gap_count\n";
			for($t=0;$t<$gap_count;$t++){
				$fill_parameter = "./fillGapRegion";
				my $left_gap_read = "leftReadSet_gapIndex_$t.fa";
				my $right_gap_read = "rightReadSet_gapIndex_$t.fa";
				
				$folder = "./gap_read_out_put_$k/";
				$left_gap_read = $folder.$left_gap_read;
				$right_gap_read = $folder.$right_gap_read;
				$fill_parameter = $fill_parameter.' '.$left_gap_read.' '.$right_gap_read.' '.$library_insertsize[$k].' '.$library_std[$k].' '.$library_orientation[$k];
				
				$fill_parameter = $fill_parameter.' '.$large_kmer_length[$k].' '.$kmer_length[$k].' '.$gap_left_contig_index[$t].' '.$gap_right_contig_index[$t].' '.$gap_distance[$t];
				print OUT "$fill_parameter\n";
				@temp = ("$fill_parameter");
				`@temp`;
			}
			
			@temp = ("./GetFinalScaffoldFill $contig_set $fill_gap_region $scaffold_end_contig_index $scaffold_set_fill_gap");
			`@temp`;
			print OUT "real_iterative_count:$real_iterative_count\n";
			$real_iterative_count++;
			$iterative_count--;
		}
	}
	unlink "graph.fa", "fill_gap_region.fa", "scaffold_end_contig_index.fa", "gap_information.fa", "end_contig_set.fa", "contig_set.fa";
	
	
	
	
	
	
	
	
	
