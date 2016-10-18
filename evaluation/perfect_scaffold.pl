use strict;
	
	my $contig_file = shift;
	my $reference_file = shift;
	my $perfect_scaffold_file = shift;
	my $gap_file = shift;
	my $nucmer_out_prefix = "temp.nucmer";
	my $nucmer_out_delta = $nucmer_out_prefix.".delta";
	my $nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
	my $nucmer_out_coords = $nucmer_out_filter.".coords";
	
	my @temp;
	@temp = ("nucmer -p $nucmer_out_prefix $reference_file $contig_file");
	`@temp`;
	@temp = ("delta-filter -i 98 -l 200 -q $nucmer_out_delta > $nucmer_out_filter");
	`@temp`;
	@temp = ("show-coords -dTlro $nucmer_out_filter > $nucmer_out_coords");
	`@temp`;
	
	@temp = ("./getPerfectScaffold $reference_file $nucmer_out_coords $perfect_scaffold_file $gap_file");
	`@temp`;
	
	unlink glob $nucmer_out_delta;
	unlink glob $nucmer_out_filter;
	unlink glob $nucmer_out_coords;
