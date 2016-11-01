use strict;
	
	my $scaffold_file_count = @ARGV - 2;
	my $original_scaffold_file = shift;
	my $gap_file_reference = shift;
	my @scaffold_fill_file;
	for(my $k = 0; $k < $scaffold_file_count; $k++){
		$scaffold_fill_file[$k] = shift;
	}
	my $min_gap_distance = 0;
	my $alignment_identity = 98;
	my @temp;
	my $contig_set = "contig_set.fa";
	my $gap_information = "gap_contig_index.fa";
	@temp = ("./splitScaffoldSet $original_scaffold_file $min_gap_distance $gap_information $contig_set");
	`@temp`;
	
	my @gap_file_name;
	
	for(my $k = 0; $k < $scaffold_file_count; $k++){
		my $position = rindex($scaffold_fill_file[$k], ".");
		if($position != -1){
			$gap_file_name[$k] = substr($scaffold_fill_file[$k], 0, $position);
			$gap_file_name[$k] = $gap_file_name[$k]."_fill_gap_region.fa";
		}else{
			$gap_file_name[$k] = $scaffold_fill_file[$k]."_fill_gap_region.fa";
		}
		
		my $nucmer_out_prefix = "contig_set.nucmer";
		my $nucmer_out_delta = $nucmer_out_prefix.".delta";
		my $nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
		my $nucmer_out_coords = $nucmer_out_filter.".coords";
		
		@temp = ("nucmer -f -p $nucmer_out_prefix $scaffold_fill_file[$k] $contig_set");
		`@temp`;
		@temp = ("delta-filter -i $alignment_identity -q $nucmer_out_delta > $nucmer_out_filter");
		`@temp`;
		@temp = ("show-coords -dTlro $nucmer_out_filter > $nucmer_out_coords");
		`@temp`;
		
		@temp = ("./getGapRegionInScaffold $gap_information $nucmer_out_coords $contig_set $scaffold_fill_file[$k] $gap_file_name[$k]");
		`@temp`;
		unlink glob $nucmer_out_delta;
		unlink glob $nucmer_out_filter;
		unlink glob $nucmer_out_coords;
	}
	
	my $gap_count = 0;
	my $line;
	open(INPUT,$gap_file_reference)or die $!;
	while($line=<INPUT>){
		$gap_count++;
	}
	$gap_count = $gap_count/2;
	
	my @temp;
	my $fscore_file = "fscore.fa";
	for(my $k = 0; $k < $scaffold_file_count; $k++){
		@temp = ("./needleman_wunsch $gap_file_reference $gap_file_name[$k] $gap_count 1 $fscore_file");
		`@temp`;
	}
	
	unlink $contig_set, $gap_information, "scaffold_end_contig_index.fa";
	
	
	
	
