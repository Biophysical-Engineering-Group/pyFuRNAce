use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;
use Encode;
use utf8;

my $number_of_designs;
my $number_start;

( $number_of_designs, $number_start ) = @ARGV;

if ( !defined $number_of_designs ) { $number_of_designs=1; }
if ( !defined $number_start ) { $number_start=1; }


my $output_file_path = "";

my $outputs = qx(cat target.txt);	 #### READ in target.txt  
my @results = split(' ', $outputs);	 #### Split the target.txt output into an array @results

my $name = "";
my $init_target = "";
my $constraint = "";
my $start_seq = "";
my $counter = 0;



foreach(@results){
	if($counter == 0){ $name = "$_";}
	if($counter == 1){ $init_target = "$_";}
	if($counter == 2){ $constraint = "$_";}
	if($counter == 3){ $start_seq = "$_";}
	$counter ++;
}



qx(mkdir data\_$name);
$output_file_path = "data\_$name";

open(my $output_spool, '>', "$output_file_path\/$name\_summary.txt"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

for(my $z=$number_start; $z<($number_of_designs+1); $z++){	
	open(my $file_out, '>', "target.txt");
		print $file_out "$name-$z\n";
		print $file_out "$init_target\n";
		print $file_out "$constraint\n";
		print $file_out "$start_seq";	
	close $file_out;
	
	print $output_spool "Design $name-$z started:\n\n";
	qx(perl revolvr.pl $output_file_path);
	
	my $read_sequence = qx(cat $output_file_path\/$name-$z\_design.txt);	 #### READ in the summary from each round 
	print $output_spool "$read_sequence\n";
}

#Done with designs, now fix the name in the target.txt
	open(my $file_out, '>', "target.txt");
		print $file_out "$name\n";
		print $file_out "$init_target\n";
		print $file_out "$constraint\n";
		print $file_out "$start_seq";	
	close $file_out;



print $output_spool "All Designs Completed.\n";
close($output_spool);