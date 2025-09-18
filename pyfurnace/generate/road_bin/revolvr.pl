use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;
use Encode;
use utf8;

#Copyright 2020  Cody Geary and Ebbe S. Andersen

#Permission is hereby granted, free of charge, to any person obtaining
#a copy of this software and associated documentation files (the
#"Software"), \ to deal in the Software without restriction, including
#without limitation the rights to use, copy, modify, merge, publish,
#distribute, sublicense, \ and/or sell copies of the Software, and to
#permit persons to whom the Software is furnished to do so, subject to
#the following conditions:

#The above copyright notice and this permission notice shall be
#included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY,\ FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\ LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALIN\ GS IN THE
# SOFTWARE.


##########
### User-Configurable Global Vars
##########

my $init_setup = 5;			  #Number of cycles of the initial setup
my $site_mutation_rate = 45;  #Rate of mutation in regions that are mispaired. Only mispaired nts are mutated.
							  #Mutations that improve the fold are retained, and become the next parent sequence.

my $mutate_loop = 2;		  #1-100, Rate at which mispaired loop sequence are mutated.
my $drift_period = 3;		  #Controls frequency of "drift", which allows mutatitions to be tested in properly folded areas.
							  #5 = drift every 5th round.  Targeted mutations do not occur during drift.

my $my_drift_rate = 2;		  #Rate for random point-mutations to the entire sequence.
my $gc_mut_ratio = 70;		  #How much GCs are preferred over AUs.
my $double_mutant = 99;		  #Percent of mutants that target base-pairs rather than single positions.
my $misfold_tolerance = 1.06; #factor of increase in distance that will be allowed during the mutation stage
my $output_threshold = 0;	  #1=verbose, 0=hide output for neutral mutations
my $max_GC=55;				  #set desired GC content of final design
my $KL_GC=45;				  #set GC for Kissing Loops
my $target_GU = 1;			  #Minimum GU, turned off since we are enforcing GU with the pattern
my $longrange = 250;		  #define the min distance for pairs to require AUC-limited alphabet in initilization.

my $duplicate_window = 10;    ### Length of Duplicates regions that Pattern Search looks for   (marked "D", default=10)
my $complement_window = 10;        #window-size to search for complement sequences - for GU placer
my $complement_zones=0;
my $duplication_zones=0;

my $max_failed_tries = 20;    #sets the number of times it will try to mutate the surrounding nts when the only remaining misfolds are in sequence-locked regions

my $generate_previews=0;       #runs trace to preview.txt every new round.  This is only helpful for debugging giant designs.

my $KL_min_delta_G = -6;		#minimum threshold for non-cognate KL interactions.

##Scoring function for KLs
my $MinKL = -7.2;    ##Default -7.2
my $MaxKL = -10.8;    ##Default -10.8


my $fav_rad_level=15;		  #'radiation' level, ~number of mutants per round, this scales up/down based on success.

my $output_file_path = "";

( $output_file_path ) = @ARGV;

$output_file_path = $output_file_path."\/";

my $start_time = time;       	 #start a timer
my $timeout_length = 7200; 		#length of timeout in seconds - default is set to 2h

my $full_seq = '';	#full sequence of the design, used for printing out the final design



##########
###Define Variables used in subroutines
##########

my $drift_rate = $my_drift_rate;
my $rad_level = $fav_rad_level;

#Initizalize global vars used in subroutines
my $fold = "";
my $distance = 0;
my @map = ();
my $randomletter = "A";
my $name = "";	#the name of the design is extracted from target.txt line 1 and stored here

###These store the folds of all the unique sequences tested
my $unique_sequence_id = 0;
my @sequence_catalog = ();
my @fold_catalog = ();
my $previously_tried=0;
my @duplex_catalog = ();
my @duplex_out_catalog = ();
my @duplex_parens_catalog =  ();
my $unique_duplex_id = 0;

#for the KL analysis sub
my $KL_A = "GCGCGC";
my $KL_B = "GUGCGU";
my $KLenergy = "";
my $parens = "";
my $numKL=0;
my @Alist = ();	 #array of KLA
my @Blist = ();	 #array of KLB
my @Elist = ();	 #array of Energies
my @NTlist = (); #store the nt position of each KL
my $KL_score = 0;
my $report_KL = 0;	#allows to toggle printed output from the KL analysis sub
my $full_KL_output =0;
my @KL_opt_target = ();	 #array of KL sites to target for optimization
my $KL_opt_target_text = "";  #the join of the array
my @parent_KL_opt_target = ();
my $numKL_sites=0;	 #total number of targets in @KL_opt_target (sites may be multiply targeted)
my $best_KL_score=1000;
my @KL_type = (); ##Store the type of KL each site is: (A)=180KL or (B)=120KL
my @KL_specified = (); ##Keep track if any KLs are user-specified, and ignore them for computing the energy calculation penalty
my $KL_notcounted = 0; ##Keep track if any KLs are user-specified, and ignore them for computing the energy calculation penalty
my $non_cognates = "";

my $failed_initial_design = 0;  #variable to track if the program needed to change the target structure to solve it

##########
###Define Subroutines for Calculation Functions
##########

my @lines;
sub read_file
{
    my ( $file ) = @_;

    my ( @lines, $content );

    if ( not open FILE, $file ) {
        die "Not readable: $file";
    }

    {
        # Reads the whole file by temporarily undefining the
        # record separator,

        local $/ = undef;
        $content = <FILE>;

        # Make Unix line ends from either Mac or DOS,

        $content =~ s/\r\n/\n/g;
        $content =~ s/\r/\n/g;
    }

    close FILE;

    @lines = split "\n", $content;

    return wantarray ? @lines : \@lines;
}


sub fold { #outputs to $fold

	###These var names store the folds of all the unique sequences tested
	# $unique_sequence_id = 0;
	# @sequence_catalog = ();
	# @fold_catalog = ();

	my $dont_fold = 0;
	$previously_tried=0;
	my $sequence_num_counter = 0;
	my $counter = 0;

	foreach(@sequence_catalog){
		if(defined $sequence_catalog[$sequence_num_counter]){
			if($_[0]eq$sequence_catalog[$sequence_num_counter]){
				$fold=$fold_catalog[$sequence_num_counter];
				$dont_fold=1;
				&printer("<found #$sequence_num_counter>\n");
				$previously_tried=1;
			}
		}
		$sequence_num_counter ++;
	}

	$sequence_num_counter = scalar(@sequence_catalog);

	if($dont_fold==0){

		$full_seq = $_[0];

		my $outputs = qx(python viennarna_funcs.py fold "$full_seq");
		chomp($outputs);         # Remove trailing newline
		$fold = "$outputs";  # Store output directly

		$sequence_catalog[$unique_sequence_id]=$_[0];
		$fold_catalog[$unique_sequence_id]=$fold;
		$unique_sequence_id ++;
	}
}

sub distance { #outputs to $distance
	my $outputs = qx(python viennarna_funcs.py distance "$_[0]" "$_[1]");
	chomp($outputs);         # Remove trailing newline
	$distance = "$outputs";  # Store output directly
}

sub map {
	@map = ();
	my @secondary = split('', $_[0]);
	my @last_bracket;
	my @last_brace;
	my @last_curlybrace;

	my @last_a_brace;
	my @last_b_brace;
	my @last_c_brace;
	my @last_d_brace;
	my @last_e_brace;
	my @last_f_brace;
	my @last_g_brace;
	my @last_h_brace;
	my @last_i_brace;
	my @last_j_brace;




	my $counter = 0;
	foreach(@secondary){
		if ($_ eq "."){
			$map[$counter] = $counter;	 ##Single strands map to themselves
		}
		elsif ($_ eq "("){
			push @last_bracket, $counter;
		}
		elsif ($_ eq ")"){
			my $partner = pop @last_bracket;
			$map[$counter]=$partner;
			$map[$partner]=$counter;
		}
		elsif ($_ eq "["){
			push @last_brace, $counter;
		}
		elsif ($_ eq "]"){
			my $buddy = pop @last_brace;
			$map[$counter]=$buddy;
			$map[$buddy]=$counter;
		}
		elsif ($_ eq "{"){
			push @last_curlybrace, $counter;
		}
		elsif ($_ eq "}"){
			my $pal = pop @last_curlybrace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "A"){    ##Outwards pointing KLs map as A-1 B-2 C-3.... I-9 J-0
			push @last_a_brace, $counter;
		}
		elsif ($_ eq "1"){
			my $pal = pop @last_a_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "B"){
			push @last_b_brace, $counter;
		}
		elsif ($_ eq "2"){
			my $pal = pop @last_b_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "C"){
			push @last_c_brace, $counter;
		}
		elsif ($_ eq "3"){
			my $pal = pop @last_c_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "D"){
			push @last_d_brace, $counter;
		}
		elsif ($_ eq "4"){
			my $pal = pop @last_d_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "E"){
			push @last_e_brace, $counter;
		}
		elsif ($_ eq "5"){
			my $pal = pop @last_e_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "F"){
			push @last_f_brace, $counter;
		}
		elsif ($_ eq "6"){
			my $pal = pop @last_f_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "G"){
			push @last_g_brace, $counter;
		}
		elsif ($_ eq "7"){
			my $pal = pop @last_g_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "H"){
			push @last_h_brace, $counter;
		}
		elsif ($_ eq "8"){
			my $pal = pop @last_h_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "I"){
			push @last_i_brace, $counter;
		}
		elsif ($_ eq "9"){
			my $pal = pop @last_i_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}
		elsif ($_ eq "J"){
			push @last_j_brace, $counter;
		}
		elsif ($_ eq "0"){
			my $pal = pop @last_j_brace;
			$map[$counter]=$pal;
			$map[$pal]=$counter;
		}


		$counter ++;
	}
}

sub randseq {  #returns $randomletter
	my $randompick = int(rand(4));
	if ($randompick == 0){$randomletter = "A";}
	if ($randompick == 1){$randomletter = "U";}
	if ($randompick == 2){$randomletter = "C";}
	if ($randompick == 3){$randomletter = "G";}
}

#arrays to track patterns and track which 150nt windows to target
my @pattern_zones = ();

##vars for the GC counter subroutine at the end
my $pairs=0;
my $GC_pairs=0;
my $AU_pairs=0;
my $GU_pairs=0;
my $KL_GC_cont=0;
my $KL_pairs=0;
my $KL_GC_pairs=0;
my $KL_repeats = 0;

##########
###Begin Program
##########


my $outputs = "target.txt";
my @results = &read_file( $outputs );

my $counter = 0;
my @target = ();
my @init_target = ();
my @constraint = ();
my @start_seq = ();
my @parent_fold = ();
my $targetstring = "";
my $targetstringpk = "";
my $generate_sequence = 'true';

foreach(@results){
	if($counter == 0){ $name = "$_";}
	if($counter == 1){ @init_target = split(//,"$_"); @target = split(//,"$_"); $targetstringpk = "$_";}
	if($counter == 2){ @constraint = split(//,"$_"); }
	if($counter == 3){ @start_seq = split(//,"$_"); $generate_sequence='false'; }
	$counter ++;
}

###############################
#Setup Output Spool
#################################
open(my $output_spool, '>', "$output_file_path$name\_spool.txt"); ##Output buffer for the screen/file outputs.  Spool will hold all screen outputs and send them to a file at the end.

&printer("Output path is $output_file_path \n");

&welcome;

if($generate_sequence eq 'false' ){&printer("I found a sequence!\n");}


&map($targetstringpk);

#######
## Parse the Target string into 2 arrays, 1 with PK and 1 without
#######

my $targetscrubbed = "";
foreach(@target){
	if ($_ eq "[" || $_ eq "]" || $_ eq "{" || $_ eq "}" ||
		$_ eq "A" || $_ eq "B" || $_ eq "C" || $_ eq "D" || $_ eq "E" || $_ eq "F" || $_ eq "G" || $_ eq "H" || $_ eq "I" || $_ eq "J" ||
		 $_ eq "1" || $_ eq "2" || $_ eq "3" || $_ eq "4" || $_ eq "5" || $_ eq "6" || $_ eq "7" || $_ eq "8" || $_ eq "9" || $_ eq "0") {
		$targetstring=$targetstring.".";
	}	## $targetstring has the PKs as .
	else {$targetstring=$targetstring.$_;}


	if ($_ eq "{"  || $_ eq "A" || $_ eq "B" || $_ eq "C" || $_ eq "D" || $_ eq "E" || $_ eq "F" || $_ eq "G" || $_ eq "H" || $_ eq "I" || $_ eq "J" ) { $targetscrubbed=$targetscrubbed."[";}	 ## $targetscrubbed has the { as [
	elsif ($_ eq "}" || $_ eq "1" || $_ eq "2" || $_ eq "3" || $_ eq "4" || $_ eq "5" || $_ eq "6" || $_ eq "7" || $_ eq "8" || $_ eq "9" || $_ eq "0") { $targetscrubbed=$targetscrubbed."]";}
	else {$targetscrubbed=$targetscrubbed.$_;}
}
@target = split(//,$targetscrubbed);   # @target has (.) and [.] but not {

my @targetnopk = split(//,$targetstring);	# @targetnopk has all the [ replaced with .

&printer("Designing $name:\n$targetstringpk\n");

my $strand_length = length $targetstring;

if (scalar(@constraint) != scalar(@target)){
	&printer("\nInvalid Constraint, using all N instead.\n");
	@constraint = ();
	for(my $i=0; $i<$strand_length; $i++){$constraint[$i]="N";}
}

&printer("The design is $strand_length nts long\n\n");

###BEGIN Design

my @best_sol = ();
my @mut_map = ();
my @defect_map = ();
my $mut_map_text = "";
my $best_dist = 1000;

my $i=0;

if($generate_sequence eq 'true'){

	##Set LOOP positions match the constraint
	for ( $i=0; $i<$strand_length; $i++){
		$best_sol[$i]="X";
		if ($target[$i] eq "."){
			if ($constraint[$i] =~ /[AUCG]/){$best_sol[$i]=$constraint[$i];}
		elsif ($constraint[$i] eq "R"){
			if(int(rand(2))==0){$best_sol[$i]="G";} else{$best_sol[$i]="A";}
		}
		elsif ($constraint[$i] eq "N"){
			&randseq;
			$best_sol[$i]=$randomletter;
		}
		elsif ($constraint[$i] eq "Y"){
			if(int(rand(2))==0){$best_sol[$i]="C";} else{$best_sol[$i]="U";}
		}
		else {$best_sol[$i]="A";} #in the odd case just use A, although this should never happen...
		}
	}

	my $GC_cont=100;
	my $bp_type=0;
	my $bp_direction=0;

	##Fill in the stems to match the constraints
	for ( $i=0; $i<$strand_length; $i++){
		if (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i] eq "A")){
			$best_sol[$i]="A";
			$best_sol[$map[$i]]="U";
		}
		elsif (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i]eq"C")){
			$best_sol[$i]="C";
			$best_sol[$map[$i]]="G";
		}
		elsif (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i]eq"U")){
			$best_sol[$i]="U";
			if ($constraint[$map[$i]]eq"A" || $constraint[$map[$i]]eq"N"){
				$best_sol[$map[$i]]="A";
			}
			elsif ($constraint[$map[$i]]eq"G" || $constraint[$map[$i]]eq"K"){
				$best_sol[$map[$i]]="G";
			}
		}
		elsif (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i]eq"G")){
			$best_sol[$i]="G";
			if ($constraint[$map[$i]]eq"C" || $constraint[$map[$i]]eq"N"){
				$best_sol[$map[$i]]="C";
			}
			elsif ($constraint[$map[$i]]eq"U" || $constraint[$map[$i]]eq"K"){
				$best_sol[$map[$i]]="U";
			}
		}
		elsif (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i]eq"K")){
			if ($constraint[$map[$i]]eq"G"){
				$best_sol[$map[$i]]="U";
			}
			elsif ($constraint[$map[$i]]eq"U"){
				$best_sol[$map[$i]]="G";
			}
			if ($constraint[$map[$i]]eq"K"){
				$bp_direction=int(rand(2));
				if( ($map[$i]-$i)>$longrange ){$bp_direction = 1;}
				if( $bp_direction==0 || $best_sol[$i]eq"G" || $best_sol[$map[$i]]eq"U"){
					$best_sol[$i]="G";
					$best_sol[$map[$i]]="U";
				}
				else {
					$best_sol[$i]="U";
					$best_sol[$map[$i]]="G";
				}
			}
		}
		elsif (($target[$i]eq("(") or $target[$i]eq("[")) && ($constraint[$i]eq"S")){
			if ($constraint[$map[$i]]eq"G"){
				$best_sol[$map[$i]]="C";
			}
			elsif ($constraint[$map[$i]]eq"C"){
				$best_sol[$map[$i]]="G";
			}
			if ($constraint[$map[$i]]eq"S"){
				$bp_direction=int(rand(2));
				if( ($map[$i]-$i)>$longrange ){$bp_direction = 1;}
				if( $bp_direction==0 || $best_sol[$i]eq"G" || $best_sol[$map[$i]]eq"C"){
					$best_sol[$i]="G";
					$best_sol[$map[$i]]="C";
				}
				else {
					$best_sol[$i]="C";
					$best_sol[$map[$i]]="G";
				}
			}
		}



		if (($target[$i]eq("(")) && ($constraint[$i]eq"N")){
			if ($constraint[$map[$i]]eq"A"){
				$best_sol[$i]="U";
				$best_sol[$map[$i]]="A";
			}
			elsif ($constraint[$map[$i]]eq"U"){
				$best_sol[$i]="A";
				$best_sol[$map[$i]]="U";
			}
			elsif ($constraint[$map[$i]]eq"C"){
				$best_sol[$i]="G";
				$best_sol[$map[$i]]="C";
			}
			elsif ($constraint[$map[$i]]eq"G"){
				$best_sol[$i]="C";
				$best_sol[$map[$i]]="G";
			}
			elsif ($constraint[$map[$i]]eq"N" && $best_sol[$i]eq"X" && $best_sol[$map[$i]]eq"X" ){
				$bp_direction= int(rand(2));
				$bp_type=(rand(100));
				if ($bp_type<(100-$max_GC-$target_GU/2)) {
					if ($bp_direction==0) {
						$best_sol[$i] = "A";
						$best_sol[$map[$i]] = "U";
					}
					else {
						$best_sol[$i] = "U";
						$best_sol[$map[$i]] = "A";
					}
				}
				else {
					if ($bp_direction==0) {
						$best_sol[$i] = "G";
						$best_sol[$map[$i]] = "C";
					}
					else {
						$best_sol[$i] = "C";
						$best_sol[$map[$i]] = "G";
					}
				}
			}
		}

		if (($target[$i]eq("[")) && ($constraint[$i]eq"N")){
			if ($constraint[$map[$i]]eq"A"){
				$best_sol[$i]="U";
				$best_sol[$map[$i]]="A";
			}
			elsif ($constraint[$map[$i]]eq"U"){
				$best_sol[$i]="A";
				$best_sol[$map[$i]]="U";
			}
			elsif ($constraint[$map[$i]]eq"C"){
				$best_sol[$i]="G";
				$best_sol[$map[$i]]="C";
			}
			elsif ($constraint[$map[$i]]eq"G"){
				$best_sol[$i]="C";
				$best_sol[$map[$i]]="G";
			}
			elsif ($constraint[$map[$i]]eq"N" && $best_sol[$i]eq"X" && $best_sol[$map[$i]]eq"X" ){
				$bp_direction= int(rand(2));
				$bp_type=(rand(100));
				if ($bp_type<$KL_GC) {
					if ($bp_direction==0) {
						$best_sol[$i] = "G";
						$best_sol[$map[$i]] = "C";
					}
					else {
						$best_sol[$i] = "C";
						$best_sol[$map[$i]] = "G";
					}
				}
				else {
					if ($bp_direction==0) {
						$best_sol[$i] = "U";
						$best_sol[$map[$i]] = "A";
					}
					else {
						$best_sol[$i] = "A";
						$best_sol[$map[$i]] = "U";
					}
				}
			}
		}

		if ($target[$i]eq("(") && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N") && (($map[$i]-$i)>$longrange)){	 ####find long-range pairs and make them A,U,C limited.
			$bp_type=int(rand(3));
			if($bp_type==0){
				$best_sol[$i] = "A";
				$best_sol[$map[$i]] = "U";
			}
			if($bp_type==1){
				$best_sol[$i] = "C";				 ##This puts AUC on the 5' most strand.	 Flip to have a AUG alphabet instead.
				$best_sol[$map[$i]] = "G";
			}
			if($bp_type==2){
				$best_sol[$i] = "U";
				$best_sol[$map[$i]] = "A";
			}
		}

		if (($init_target[$i]eq("{")) && ($constraint[$i]eq"N")){	##apply 100%GC to the dovetails specified with '{'
			if ($constraint[$map[$i]]eq"A"){
				$best_sol[$i]="U";
				$best_sol[$map[$i]]="A";
			}
			elsif ($constraint[$map[$i]]eq"U"){
				$best_sol[$i]="A";
				$best_sol[$map[$i]]="U";
			}
			elsif ($constraint[$map[$i]]eq"C"){
				$best_sol[$i]="G";
				$best_sol[$map[$i]]="C";
			}
			elsif ($constraint[$map[$i]]eq"G"){
				$best_sol[$i]="C";
				$best_sol[$map[$i]]="G";
			}
			elsif ($constraint[$map[$i]]eq"N"){
				$bp_direction= int(rand(2));
				if ($bp_direction==0) {
					$best_sol[$i] = "G";
					$best_sol[$map[$i]] = "C";
				}
				else {
					$best_sol[$i] = "C";
					$best_sol[$map[$i]] = "G";
				}
			}
		}
	}
}else{
	&printer("Inserting new Sequence!\n");
	@best_sol = @start_seq;
}

&printer("Target: \n");
&printer("$targetstring\n");
my $constraintstring = join "",@constraint;
&printer("\nConstraint\n");
&printer("$constraintstring\n\n\n");
my $best_sol_seq = join "", @best_sol;

&fold($best_sol_seq);
@parent_fold = split(//,$fold);
#&printer("parent fold =	@parent_fold\n");
my $parent_fold_string = $fold;
&distance($parent_fold_string, $targetstring);
my $parent_dist = $distance;

&printer("Initial Sequence: Distance = $parent_dist");
&countGC;
&printer("$best_sol_seq\n");
&printer($fold);

my @trial_sol = @best_sol;
my $new_trial_seq = $best_sol_seq;

my $bp_direction=0; my $bp_type=0; my $GC_cont=100;

################################################
##BEGIN first optimization loop.

my $attempts = 0;

do{	 ####

	if(($distance==0) or ($generate_sequence eq 'false')){$init_setup=0;}
		$attempts += 1;	 #Keep track of how many times this loops and fails
		my $helixend=0;
	for ( $i=0; $i<$strand_length-1; $i++){
		if ($constraint[$i]eq"N" && $constraint[$map[$i]]eq"N"){
			$bp_direction=int(rand(2));	#roll the dice to pick the direction of the bp
			$bp_type=(rand(100));
			if( $target[$i]eq"." && $parent_fold[$i]ne"." ){$trial_sol[$i]="A";}
			if( ($target[$i]eq"(" || $target[$i]eq"[") && $target[$i+1]ne"(" && $target[$i+1]ne"[" ){ $helixend=1;}
			if( ($target[$i]eq")" || $target[$i]eq"]") && $target[$i+1]ne")" && $target[$i+1]ne"]" ){ $helixend=1;}
			if($helixend==1){$bp_type=100;}
			if (($target[$i]eq"(") && ($parent_fold[$i]ne"(")){
				if ($bp_type<(100-$max_GC-$target_GU)){
					if ($bp_direction==0) {
						$trial_sol[$i] = "A";
						$trial_sol[$map[$i]] = "U";
					}
					else {
						$trial_sol[$i] = "U";
						$trial_sol[$map[$i]] = "A";
					}
				}
				else {
					if ($bp_direction==0) {
						$trial_sol[$i] = "G";
						$trial_sol[$map[$i]] = "C";
					}
					else {
						$trial_sol[$i] = "C";
						$trial_sol[$map[$i]] = "G";
					}
				}
			}
			if (($target[$i]eq")") && ($target[$i+1]eq"(") && ($parent_fold[$i]ne")") && ($parent_fold[$i+1]ne"(") ){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "C";
				$trial_sol[$i+1] = "G";
			}
			if (($target[$i]eq"(") && ($target[$i-1]eq")") && ($parent_fold[$i]ne"(") && ($parent_fold[$i-1]ne")")){
				$trial_sol[$map[$i]] = "C";
				$trial_sol[$i] = "G";
				$trial_sol[$map[$i-1]] = "G";
				$trial_sol[$i-1] = "C";
			}
			if (($target[$i]eq"(") && ($target[$i+1]eq"(") && (($map[$i]-$map[$i+1])>1) && ($parent_fold[$i]ne"(") && ($parent_fold[$i+1]ne"(")){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "G";
				$trial_sol[$i+1] = "C";
			}
			if (($target[$i]eq")") && ($target[$i+1]eq")") && (($map[$i]-$map[$i+1])>1) && ($parent_fold[$i]ne")") && ($parent_fold[$i+1]ne")")){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "G";
				$trial_sol[$i+1] = "C";
			}
		}	elsif ($constraint[$i]eq"S" && $constraint[$map[$i]]eq"S"){
			if (($target[$i]eq")") && ($target[$i+1]eq"(") && ($parent_fold[$i]ne")") && ($parent_fold[$i+1]ne"(") ){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "C";
				$trial_sol[$i+1] = "G";
			}
			if (($target[$i]eq"(") && ($target[$i-1]eq")") && ($parent_fold[$i]ne"(") && ($parent_fold[$i-1]ne")")){
				$trial_sol[$map[$i]] = "C";
				$trial_sol[$i] = "G";
				$trial_sol[$map[$i-1]] = "G";
				$trial_sol[$i-1] = "C";
			}
			if (($target[$i]eq"(") && ($target[$i+1]eq"(") && (($map[$i]-$map[$i+1])>1) && ($parent_fold[$i]ne"(") && ($parent_fold[$i+1]ne"(")){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "G";
				$trial_sol[$i+1] = "C";
			}
			if (($target[$i]eq")") && ($target[$i+1]eq")") && (($map[$i]-$map[$i+1])>1) && ($parent_fold[$i]ne")") && ($parent_fold[$i+1]ne")")){
				$trial_sol[$map[$i]] = "G";
				$trial_sol[$i] = "C";
				$trial_sol[$map[$i+1]] = "G";
				$trial_sol[$i+1] = "C";
			}
		}
	}

	##Compute the fold and distance

	$new_trial_seq = join "",@trial_sol;	#Make a string to hold the sequence of the new trial

	&fold($new_trial_seq);
	my $new_fold = $fold;			  #The fold of puzzle solution
	&distance($new_fold, $targetstring);	 #Distance the puzzle solution is away from the target
	my $trial_dist = $distance;
	#Decide if to keep the new sequence
	my $keep_sequence = 0;
	if ($trial_dist<($parent_dist*$misfold_tolerance)){
		$keep_sequence=1;
	}
	else{
		$keep_sequence=0;
	}

	##Process Ouputs
	if ($keep_sequence==1){		   #Check if the new mutant is better than the parent

			@mut_map = ();
			@defect_map = ();
		for ( $i=0; $i<$strand_length; $i++){
			if (($best_sol[$i]) eq ($trial_sol[$i])) {
					$mut_map[$i] = "-";
			} else {
				$mut_map[$i] = $trial_sol[$i];
			}
			if ($parent_fold[$i]ne$targetnopk[$i]){
				$defect_map[$i] = "X";
			} else {
				$defect_map[$i] = "_";
			}
			$best_sol[$i]= $trial_sol[$i];		#move the trial solution to best_sol array
		}

		my $mut_map_text = join "",@mut_map;
		my $defect_map_text = join "",@defect_map;
		$best_sol_seq = join "", @trial_sol;

		$parent_dist=$trial_dist;			  #copy the trial distance to be the new parent_dist
		my $best_dist=$trial_dist;
		$parent_fold_string = $fold;
		@parent_fold = split(//, $fold);

		&printer("\nIteration# $attempts, Distance = $trial_dist");	#Report outputs
		&countGC;
		if ($attempts>1){&printer("$mut_map_text\n");}
		&printer("$best_sol_seq\n");
		&printer("$parent_fold_string\n");
		&printer("$defect_map_text\n");

	}
	else {
		&printer("( $attempts )");
		for ( $i=0; $i<$strand_length; $i++){
			$trial_sol[$i]=$best_sol[$i];  #reset the trial solution
		}
	}
}while($attempts < $init_setup);

################################################
&printer("\n__________________________________________________________\n");
&printer("Begin Directed Mutation Cycle\n");

my $drift_event_count = 0;	#initialize the counter for mutation events
my $improved_parent_count = 0;
my $constraint_check = 0;
my $failed_tries =0;
do{

	my $keep_sequence=0;
	$attempts += 1;
	$improved_parent_count += 1;
	$drift_event_count += 1;

	my @new_fold_array = ();
	my $trial_dist=10000;

	my $mutate_now=0;
	my @mutation_sites = ();
	my $coin_flip = 0;
	my $num_mutation_spots =0;

	&printer("($attempts)");
	    if($drift_event_count==$drift_period){&printer("drift->");}
	    if($drift_event_count!=$drift_period){&printer("target->");}
	    if($improved_parent_count>1 && $drift_event_count!=$drift_period){&printer("expanded->");}
		for ( $i=2; $i<$strand_length-2; $i++){
			$mutate_now=0;
			$coin_flip=int(rand(2));
			if (($drift_event_count==$drift_period) && (rand(100)<$drift_rate) && ($constraint[$i]eq"N" || $constraint[$i]eq"K" || $constraint[$i]eq"S")) {	 #random drift events
				$mutate_now=1;
			}
			if (($parent_fold[$i]ne$targetnopk[$i]) && ($drift_event_count!=$drift_period) && (rand(100)<$site_mutation_rate) && ($constraint[$i]eq"N" || $constraint[$i]eq"K" || $constraint[$i]eq"S") ) {	 #targeted mutations
				$mutate_now=1;
			}
			if ($improved_parent_count>1) {
				if ( (($parent_fold[$i+1] ne $targetnopk[$i+1]) || ($parent_fold[$i-1] ne $targetnopk[$i-1]) || ($parent_fold[$i+2] ne $targetnopk[$i+2]) || ($parent_fold[$i-2] ne $targetnopk[$i-2]) ) &&
				($drift_event_count!=$drift_period) && (rand(100)<$site_mutation_rate) && $constraint[$i]eq"N") {	#expand mutational window by 2nts
					$mutate_now=1;
				}

				if($parent_fold[$i]ne$target[$i] && $target[$i]eq")" && $target[$i+1]eq"("){
					$mutation_sites[$i]=1; $mutation_sites[$map[$i]]=1;
					$mutation_sites[$i+1]=1; $mutation_sites[$map[$i+1]]=1;
					$mutation_sites[$map[$i]-1]=1;
					$mutation_sites[$map[$map[$i]-1]]=1;
					$mutate_now=1;
				}

				if($parent_fold[$map[$i]]ne$target[$map[$i]] || $parent_fold[$map[$i+1]]ne$target[$map[$i+1]] || $parent_fold[$map[$i-1]]ne$target[$map[$i-1]]){
					$mutation_sites[$map[$i]]=1; $mutation_sites[$map[$i+1]]=1; $mutation_sites[$map[$i-1]]=1;
					$mutate_now=1;
				}
			}

			if ($mutate_now==1){
				$mutation_sites[$i]=1; $num_mutation_spots++;}
			else{$mutation_sites[$i]=0;}
		}

	my $mut_scale = $rad_level/($num_mutation_spots+.001);	##scale back the odds of mutation so that there are always ~$radlevel mutations.
	&printer(" M[$num_mutation_spots] (*$rad_level) ");

	for ( $i=2; $i<$strand_length-2; $i++){
		my $randnum=rand(1);

		if(($mut_scale>$randnum) && $mutation_sites[$i]==1){

			if ( (($target[$i]eq"(") || ($target[$i]eq")") ) && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N")){
				&printer(" $i");
				if (rand(100)<$double_mutant){
					&printer(":");
					if (rand(100)<$gc_mut_ratio) {
						if (int(rand(2))==1) {
							$trial_sol[$i] = "G";
							$trial_sol[$map[$i]] = "C";
						}else{
							$trial_sol[$i] = "C";
							$trial_sol[$map[$i]] = "G";
						}
					}else{
						if (int(rand(2))==1) {
							$trial_sol[$i] = "A";
							$trial_sol[$map[$i]] = "U";
						}else{
							$trial_sol[$i] = "U";
							$trial_sol[$map[$i]] = "A";
						}
					}
				}else {
					&printer(".");
					if (($best_sol[$i] eq "G") && ($best_sol[$map[$i]] eq "C")){
						$trial_sol[$map[$i]]= "U";
					}
					elsif (($best_sol[$i] eq "C") && ($best_sol[$map[$i]] eq "G")){
						$trial_sol[$i]= "U";
					}
					elsif (($best_sol[$i] eq "A") && ($best_sol[$map[$i]] eq "U")){
						$trial_sol[$i]= "G";
					}
					elsif (($best_sol[$i] eq "U") && ($best_sol[$map[$i]] eq "A")){
						$trial_sol[$map[$i]]= "G";
					}
					elsif (($best_sol[$i] eq "G") && ($best_sol[$map[$i]] eq "U")){
						if ($coin_flip==0){$trial_sol[$map[$i]]= "C"; }
						else {$trial_sol[$i]= "A"; }
					}
					elsif (($best_sol[$i] eq "U") && ($best_sol[$map[$i]] eq "G")){
						if ($coin_flip==0){$trial_sol[$i]= "C"; }
						else {$trial_sol[$map[$i]]="A"; }
					}
				}
			}#endif for	"("
			elsif (($target[$i]eq".") && ($constraint[$i]eq"N")){
				&printer("L");
				if (rand(100) < ($mutate_loop)){
					&randseq();
					$trial_sol[$i]= $randomletter;
				}
			}
			elsif (($target[$i]eq".") && ($constraint[$i]eq"R")){
				&printer("Lr");
				if (int(rand(2))==0){$trial_sol[$i]= "A";}else{$trial_sol[$i]= "G";}
			}

			if ( (($target[$i]eq"[") || ($target[$i]eq"]")) && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N") ) {
				&printer(" $i");	 &printer("[:");
				if (rand(100)<($KL_GC)) {
					if (int(rand(2))==1) {
						$trial_sol[$i] = "G";
						$trial_sol[$map[$i]] = "C";
					}else{
						$trial_sol[$i] = "C";
						$trial_sol[$map[$i]] = "G";
					}
				}
				else{
					if (int(rand(2))==1) {
						$trial_sol[$i] = "A";
						$trial_sol[$map[$i]] = "U";
					}else{
						$trial_sol[$i] = "U";
						$trial_sol[$map[$i]] = "A";
					}
				}
			} #endif	 "["

		  if ( (($target[$i]eq"(") || ($target[$i]eq"[") || ($target[$i]eq")") || ($target[$i]eq"]")) && ($constraint[$i]eq"K") && ($constraint[$map[$i]]eq"K")){
			 &printer("K");
			 if ($best_sol[$i] eq "G"){
			   $trial_sol[$map[$i]]= "G";
			   $trial_sol[$i]="U";
			 }
			 elsif ($best_sol[$i] eq "U"){
			   $trial_sol[$i]= "G";
			   $trial_sol[$map[$i]]= "U";
			 }
		  }

		  if ( (($target[$i]eq"(") || ($target[$i]eq"[") || ($target[$i]eq")") || ($target[$i]eq"]")) && ($constraint[$i]eq"S") && ($constraint[$map[$i]]eq"S")){
			 &printer("S");
			 if ($best_sol[$i] eq "G"){
			   $trial_sol[$map[$i]]= "G";
			   $trial_sol[$i]="C";
			 }
			 elsif ($best_sol[$i] eq "C"){
			   $trial_sol[$i]= "G";
			   $trial_sol[$map[$i]]= "C";
			 }
		  }


		} #endif for mutscale

	} #end for loop


	if (($drift_event_count<=$drift_period)){  #reset drift counter
	   $drift_event_count=0;
	}

	##Compute Distances

	$new_trial_seq = join "",@trial_sol;	#Make a string to hold the sequence of the new trial
	&fold($new_trial_seq);
	my $new_fold = $fold;			  #The fold of the new_trial
	@new_fold_array = split(//,$fold);
	&distance($new_fold, $targetstring);
	$trial_dist = $distance;   #Distance the puzzle solution is away from the target

	#Decide if we keep the new sequence
	$keep_sequence=0;
	if ($trial_dist<$parent_dist) {$improved_parent_count = 0;}
	if ($trial_dist<($parent_dist*$misfold_tolerance)) {$keep_sequence=1;}

	#Process outputs

	if ($keep_sequence==1 && $previously_tried==0){		   #Check if the new mutant is better than the parent
		$rad_level += 1;
		for ( $i=0; $i<$strand_length; $i++){
			if (($best_sol[$i]) eq ($trial_sol[$i])) {
				$mut_map[$i] = "-";
			} else {
				$mut_map[$i] = $trial_sol[$i];
			}
			if ($new_fold_array[$i]ne$targetnopk[$i]){
				$defect_map[$i] = "X";
			} else {
				$defect_map[$i] = "-";
			}
			$best_sol[$i]= $trial_sol[$i];		#move the trial solution to best_sol array
		}

		$mut_map_text = join "",@mut_map;
		my $defect_map_text = join "",@defect_map;
		$best_sol_seq = join "",@best_sol;

		my $display_the_output = 0;

####
# New subroutine, if the defect map aligns with the constrained sequence, then move on in the design.

		$constraint_check = 0;
		for ( $i=0; $i<$strand_length; $i++){
			if ($defect_map[$i]eq"X" && ($constraint[$i]eq"N"||$constraint[$i]eq"K"||$constraint[$i]eq"S")){ $constraint_check += 1;}  #if there are no defects, or no defects that are not masked by constraints then $constraint_check will be 0
		}


		if($trial_dist<$parent_dist || $output_threshold==1 || $constraint_check==0){$display_the_output = 1;}

		$parent_dist=$trial_dist;			  #copy the trial distance to be the new parent_dist
		$best_dist=$trial_dist;
		@parent_fold = split(//,$fold);
		$parent_fold_string = $fold;

		if($display_the_output==1 && $previously_tried==0){
			&printer("\n\nIteration# $attempts, Distance = $trial_dist\n");  #Report outputs
			&countGC;
			&printer("\nMutations:\n");
			&printer("$mut_map_text\n");
			&printer("$best_sol_seq\n");
			&printer("$parent_fold_string\n");
			&printer("Defect Map:\n");
			&printer("$defect_map_text\n");

			if ($generate_previews==1){&preview;}

		} else {&printer(" keep -D$trial_dist");}
		$drift_event_count=0;
	}else { ##else if keep_sequence==0;
		$rad_level += -1;
		if ($rad_level <1){$rad_level=1;}
		&printer(" dump -D$trial_dist");
		for ( $i=0; $i<$strand_length; $i++){
			$trial_sol[$i]=$best_sol[$i];  #reset the trial solution
		}
	}
	if($previously_tried==1){$drift_event_count=($drift_period-1); &printer("Forcing Drift"); $rad_level++;} #force sequence drift if stuck in a rut
	&printer("\n");

	if ($constraint_check==0 && $distance>0 && $keep_sequence==1){
		$failed_tries++;
		$drift_period=1;
		&printer(" <Blocked by Constraints - $failed_tries > \n");
	}

	if ($constraint_check==0 && $distance>0 && $failed_tries>$max_failed_tries  && $keep_sequence==1){    ##Try to random drift at least $max_failed_tries times before giving up.
		&printer("All defects are locked by the sequence constraint, resetting target to current best defect mask\n");
		$distance=0;

		for ( $i=0; $i<$strand_length; $i++){
		    if( ($target[$i]ne"[") && ($target[$i] ne "]") ){$target[$i]=$new_fold_array[$i];}  #leave PKs in place
			$targetnopk[$i]=$new_fold_array[$i];
			$targetstring = join "",@new_fold_array;
		}
		$failed_initial_design = 1;
		&printer("WARNING: Target has been reset to a best-effort structure after $max_failed_tries attempts.\n");
	}


}while ($distance > 0);	#Cycle through sequence design again until the target is within the threshold distance from the fold.


&printer("Success!\nDesign $name folded.\n");

&pattern_prevent;

################################################

&printer("\n__________________________________________________________\n");
&printer("Begin GC-Reduction & GC/UA-Rich Reduction\n");  #report progress

$output_threshold = 1;
$rad_level = $fav_rad_level;
my $pattern_repeats = 1000;
my $poly_repeats = 1000;
my $restriction_sites = 1000;
my $happy=0;  #not happy!
my $last_pattern_repeats=$pattern_repeats;
my $last_poly_repeats=$poly_repeats;
my $last_restriction_sites =$restriction_sites;
my $last_complement_zones = $complement_zones;
my $last_duplication_zones = $duplication_zones;

###Pre-check if the cycle is satisfied

my	@repeat_map = ();
my	$repeat_map_text = "";

for ( $i=0; $i<$strand_length; $i++){
  $trial_sol[$i]=$best_sol[$i];	 #reset the trial solution
}

&countrepeats;

&printer("\nThere are $complement_zones Complementary Regions, $pattern_repeats Strong/Weak Regions, $poly_repeats Poly-N5 Repeats, $duplication_zones nt Duplicated sequence, and $restriction_sites common restriction sites present to start:\n");
$last_pattern_repeats=$pattern_repeats;
$last_poly_repeats=$poly_repeats;
$last_restriction_sites = $restriction_sites;
$last_complement_zones = $complement_zones;
$last_duplication_zones = $duplication_zones;

###END check

my $test_first =1;

$drift_rate=$my_drift_rate; #reset drift rate here

do{

	my $keep_sequence=0;

	my $mutate_now = 0;
	my $parent_fold = "";
	my $coin_flip =0;
	my $rand_num_one =0;
	my $rand_num_two =0;

	@repeat_map = ();
	$repeat_map_text = "";

    &countrepeats;

	my @au_bust=();
	my @gu_put =();
	my @mutation_sites =();  ##array storting the targeted mutations
	my $num_mutation_spots=0;


	if($test_first == 0){ ##test to see if this is the first time running this loop

		$attempts += 1;

		do{
			&printer("(Iteration $attempts)");
			$num_mutation_spots=0;


			for ( $i=0; $i<$strand_length; $i++){
				$trial_sol[$i]=$best_sol[$i];	     #reset the trial solution
				$gu_put[$i]=0;
				$mutation_sites[$i]=0;##zero out the arrays that track which type of mutation to put where
			}

			for ( $i=0; $i<$strand_length; $i++){

				$mutate_now=0;
				$coin_flip=int(rand(2));
				$rand_num_one=rand(100);
				$rand_num_two=rand(100);


				if( ($i>3) && ($i<($strand_length-3)) ){
					if ( ( (($target[$i-1]eq("(") || $target[$i-1]eq("["))) && (($target[$i]eq("(") || $target[$i]eq("["))) && (($target[$i+1]eq("(") || $target[$i+1]eq("["))) ) &&
					  (($trial_sol[$i-1]eq"C") || ($trial_sol[$i-1]eq"A") || (($trial_sol[$i-1]eq"G") && ($trial_sol[$map[$i-1]]eq"C")) || (($trial_sol[$i-1]eq"U") && ($trial_sol[$map[$i-1]]eq"A"))) &&
					  (($trial_sol[$i+1]eq"C") || ($trial_sol[$i+1]eq"A") || (($trial_sol[$i-1]eq"G") && ($trial_sol[$map[$i-1]]eq"C")) || (($trial_sol[$i-1]eq"U") && ($trial_sol[$map[$i-1]]eq"A")))	 ) {  #dont' put 2 GU next to each other
						if ( (int(rand(100/($target_GU/10+1)))==0) && ($GU_pairs<($target_GU*$pairs*.01)) ) {
							$gu_put[$i]=1;    ##this allows GUs
							$mutate_now=1;
							$mutation_sites[$i]+=1
						}

						if ($repeat_map[$i]eq"P" && rand($complement_zones)<15){   #the odds of placing a GU to be roughly 15 GU per total number of sites
							$gu_put[$i]=1;   ##this allows GUs for sites marked P
							$mutate_now=1;
							$mutation_sites[$i]+=1;
						}
					}
				}

				if ($repeat_map[$i]eq"-"){  ##if "-" then neutral drift stems
					if (($target[$i]eq("(") or $target[$i]eq("[")) && $rand_num_one<($drift_rate/2)) {  #any helix or KL is targeted for mutation, drift rate is 1/2 speed
						$mutate_now=1;
						$mutation_sites[$i]+=1
					}
				} elsif ($repeat_map[$i]eq"U"){
					$mutate_now=1;
					$mutation_sites[$i]+=2;  ##get rid of the poly-Us first.
				} elsif ($repeat_map[$i]eq"G" || $repeat_map[$i]eq"C" || $repeat_map[$i]eq"A"){
					$mutate_now=1;
					$mutation_sites[$i]+=1;
				} elsif ($repeat_map[$i]eq"S" || $repeat_map[$i]eq"W" || $repeat_map[$i]eq"X" || $repeat_map[$i]eq"D"){
					if ($coin_flip==0){
						$mutate_now=1;
						$mutation_sites[$i]+=1;
					}
				} elsif ($repeat_map[$i]eq"P") {  #the odds of random mutation to be roughly 5 per total number of sites
					if (rand($complement_zones)<5){
						$mutate_now=1;
						$mutation_sites[$i]+=1;
					}
				}

				if ($mutate_now==1){$num_mutation_spots++;}
			}

			if ($rad_level < 2){$rad_level = 2;} #set a min rad level
			if ($rad_level > 20){$rad_level = 20;} #set a max rad level

			my $mut_scale = $rad_level/($num_mutation_spots+.001);	##scale back the odds of mutation so that there are always ~radlevel mutations.
			&printer(" M[$num_mutation_spots] (*$rad_level) ");

			for ( $i=0; $i<$strand_length-1; $i++){

				if( rand(1) < ($mut_scale*$mutation_sites[$i]) ){
					if ($target[$i]eq"." && $constraint[$i]eq"N"){
						&printer(".");
						&randseq;
						$trial_sol[$i] = $randomletter;
					}

					if ($repeat_map[$i]eq"P"  && $gu_put[$i]==1 && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N")){   ##added GU_put checking to avoid too many tandem GUs, #bugfix, added constraint checking!
								&printer("+W ");	  #add GU
								if (rand(1)<0.5) {
									$trial_sol[$i] = "G";
									$trial_sol[$map[$i]] = "U";
								}else{
									$trial_sol[$i] = "U";
									$trial_sol[$map[$i]] = "G";
								}
					}
					if (($repeat_map[$i]eq"G" || $repeat_map[$i]eq"C") && rand(1)<0.5 &&
						 ( $target[$i]eq("(") || $target[$i]eq(")") || $target[$i]eq("[") || $target[$i]eq("]") ) && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N") ){
						&printer(":fS "); #flip strong
						if ($trial_sol[$i]eq"C") {
							$trial_sol[$i] = "G";
							$trial_sol[$map[$i]] = "C";
						}else{
							$trial_sol[$i] = "C";
							$trial_sol[$map[$i]] = "G";
						}
					}
					elsif(($target[$i]eq("(") ||$target[$i]eq(")") || $target[$i]eq("[") || $target[$i]eq("]") ) && ($constraint[$i]eq"N") && ($constraint[$map[$i]]eq"N") ){
						&printer(":"); #bp change
						if (($trial_sol[$i]eq"C" && $trial_sol[$map[$i]]eq"G") || ($trial_sol[$i]eq"G" && $trial_sol[$map[$i]]eq"C") ){   #GC Reducer
							if((($GC_pairs/$pairs*100)>$max_GC) || rand(1)<0.5){
								&printer("+A ");	   #add AU
								if ($coin_flip==0) {
									$trial_sol[$i] = "A";
									$trial_sol[$map[$i]] = "U";
								}else{
									$trial_sol[$i] = "U";
									$trial_sol[$map[$i]] = "A";
								}
							}
						}
						elsif (($trial_sol[$i]eq"A" && $trial_sol[$map[$i]]eq"U") || ($trial_sol[$i]eq"U" && $trial_sol[$map[$i]]eq"A") ){	#AU Reducer
							if((($GC_pairs/$pairs*100)<$max_GC) || rand(1)<0.5){
								&printer("+G ");	  #add GC
								if ($coin_flip==0) {
									$trial_sol[$i] = "G";
									$trial_sol[$map[$i]] = "C";
								}else{
									$trial_sol[$i] = "C";
									$trial_sol[$map[$i]] = "G";
								}
							}
						}
						elsif (($trial_sol[$i]eq"G" && $trial_sol[$map[$i]]eq"U") || ($trial_sol[$i]eq"U" && $trial_sol[$map[$i]]eq"G") ){	#GU Flipper/Deleter

							if(rand(1)<0.9){
								&printer(":fW "); #f, for flip wobble
								if ($trial_sol[$i] eq "G"){	##GU flipper
									$trial_sol[$i] = "U";
									$trial_sol[$map[$i]] = "G";
								} else {
									$trial_sol[$i] = "G";
									$trial_sol[$map[$i]] = "U";
								}
							} else {
								&printer(":xW "); #f, for delete wobble
								if ($trial_sol[$i] eq "G"){	##GU to AU
									$trial_sol[$i] = "A";
								} else {
									$trial_sol[$map[$i]] = "A";
								}
							}

						}

					} elsif (($target[$i]eq("(") || $target[$i]eq(")") || $target[$i]eq("[") || $target[$i]eq("]") ) && ($constraint[$i]eq"K") && ($constraint[$map[$i]]eq"K") ){
						&printer(":fW "); #f, for flip wobble
						if ($trial_sol[$i] eq "G"){	##GU flipper
							$trial_sol[$i] = "U";
							$trial_sol[$map[$i]] = "G";
						} else {
							$trial_sol[$i] = "G";
							$trial_sol[$map[$i]] = "U";
						}
					} elsif (($target[$i]eq("(") || $target[$i]eq(")") || $target[$i]eq("[") || $target[$i]eq("]") ) && ($constraint[$i]eq"S") && ($constraint[$map[$i]]eq"S") ){
						&printer(":fS "); #f, for flip strong
						if ($trial_sol[$i] eq "G"){	##GC flipper
							$trial_sol[$i] = "C";
							$trial_sol[$map[$i]] = "G";
						} else {
							$trial_sol[$i] = "G";
							$trial_sol[$map[$i]] = "C";
						}
					} #end mutate pairs
				}
			}
			#####Compute Distances

			###Count up the patterns and make a map

			&countrepeats;
	        &printer("(SW$pattern_repeats/N$poly_repeats/R$restriction_sites/W$complement_zones)");
			if(($pattern_repeats + $poly_repeats + $restriction_sites + $complement_zones/$complement_window + $duplication_zones) > ($last_pattern_repeats + $last_poly_repeats + $last_restriction_sites + $last_complement_zones/$complement_window + $last_duplication_zones )){&printer(" - bad roll, trying again...\n"); $rad_level += -1;}


		}while( (($pattern_repeats + $poly_repeats + $restriction_sites + $complement_zones/$complement_window + $duplication_zones) > ($last_pattern_repeats + $last_poly_repeats + $last_restriction_sites + $last_complement_zones/$complement_window + $last_duplication_zones)) && ($pattern_repeats + $poly_repeats + $restriction_sites + $complement_zones + $duplication_zones));

	}#end test_first comparison
	$test_first=0;

	$new_trial_seq = join "",@trial_sol;	#Make a string to hold the sequence of the new trial
	#&printer("(SW$pattern_repeats/N$poly_repeats/R$restriction_sites/W$complement_zones) ... ");
	&printer("SUCCESS folding...");
	&fold($new_trial_seq);
	my $new_fold = $fold;			  #The fold of puzzle solution
	&distance($new_fold, $targetstring);
	my $trial_dist = $distance;	  #Distance the puzzle solution is away from the target


	my $print_output=0;


	if ($distance == 0 && (($pattern_repeats <= $last_pattern_repeats) || ($poly_repeats <= $last_poly_repeats) || ($restriction_sites <= $last_restriction_sites)  || ($complement_zones > $last_complement_zones ) || ($duplication_zones > $last_duplication_zones ) )) {$keep_sequence=1;}

	if ($keep_sequence==1){		   #Check if the new mutant is better than the parent
		$last_pattern_repeats=$pattern_repeats;
		$last_poly_repeats=$poly_repeats;
		$last_restriction_sites = $restriction_sites;
		$last_complement_zones = $complement_zones;
		$last_duplication_zones = $duplication_zones;

		$rad_level += 2;	#ramp up confidence fast! things fail a lot at this stage.
		for ( $i=0; $i<$strand_length; $i++){
			if (($best_sol[$i]) eq ($trial_sol[$i])) {
				$mut_map[$i] = "-";
			} else {
				$mut_map[$i] = $trial_sol[$i];
			}
			$best_sol[$i]= $trial_sol[$i];		#move the trial solution to best_sol array
		}

		$mut_map_text = join "",@mut_map;
		$repeat_map_text = join "",@repeat_map;
		$best_sol_seq = join "",@best_sol;

		$parent_dist=$trial_dist;			  #copy the trial distance to be the new parent_dist
		$best_dist=$trial_dist;
		@parent_fold = split(//,$fold);
		$parent_fold_string = $fold;


		&countGC();
		&printer(" 8S|8W Repeats: $pattern_repeats	  5G|5C|5A|5U Repeats: $poly_repeats   Restriction: $restriction_sites  Complement: $complement_zones  Duplication: $duplication_zones\n");
		&printer("Mutation Map:\n$mut_map_text\n");
		&printer("$best_sol_seq\n");
		&printer("Patterns Remaining:\n$repeat_map_text\n\n");

		if ($generate_previews==1){&preview;}

		} else {	##else don't keep sequence
			$rad_level += -1;  #decrease confidence a little bit every time it fails
			if ($rad_level <3){$rad_level=3;}  ##set the minimum to 3, so that things don't slow down too much
			&printer(" Bad Fold: -D$trial_dist");
			for ( $i=0; $i<$strand_length; $i++){
				$trial_sol[$i]=$best_sol[$i];  #reset the trial solution
		}
	}

	if ($GC_cont<$max_GC && $pattern_repeats==0 && $distance==0 && $poly_repeats==0  && $restriction_sites==0 && $complement_zones==0 && $duplication_zones==0) {
		$happy = 1;
	} elsif($GC_cont>$max_GC){&printer( "GC content still too high \n");
	}

	if ($keep_sequence==1){
		my $unmasked_pattern = 0;
		my $masked_pattern =0;
		for($i=0;$i<$strand_length; $i++){
			if( ($constraint[$i]eq"N" || $constraint[$i]eq"S" || $constraint[$i]eq"K" ) && ( $repeat_map[$i]ne"-" ) ){ $unmasked_pattern+=1;  }
			if( ($constraint[$i]eq"A" || $constraint[$i]eq"U" || $constraint[$i]eq"C" || $constraint[$i]eq"G" ) && ( $repeat_map[$i]ne"-" ) ){ $masked_pattern+=1;  }
		}
		if($unmasked_pattern==0 && $masked_pattern>0){
			&printer("\n\n Remaining patterns are masked by sequence constraints.  Passing to next phase. \n\n");
			$happy=1;
		}
	}


}while ($happy==0);	#Cycle through sequence design again until the target is within the threshold distance from the fold.

&printer("\nTarget conditions met.\n");
&printer("Total attempts: $attempts\n");
&printer("$best_sol_seq\n");

$report_KL=1;
&pattern_prevent;

##############################################################################################

&printer("\n__________________________________________________________\n");
&printer("Removing KL Repeats\n");	 #report progress

$happy=0;

	do{

	my $keep_sequence=0;
	$attempts += 1;

	&printer("($attempts)");

	my $mutate_now = 0;
	my $parent_fold = "";
	my $coin_flip =0;
	my $rand_num_one =0;
	my $rand_num_two =0;
	my $old_KL_repeats=$KL_repeats;

	for ( $i=0; $i<$strand_length; $i++){	 ##look for the repeated KLs and mutate them

		if (($target[$i]eq"[" || $target[$i]eq"]") && $constraint[$i]eq"N" && $pattern_zones[$i]ne"-"){
			$rand_num_one = rand(100);
			$bp_direction= int(rand(2));
			$bp_type=(rand(100));
			if($rand_num_one<100*(1/(4*$KL_repeats))){
				&printer("!");
				if ($bp_type<$KL_GC) {
					if ($bp_direction==0) {
						$trial_sol[$i] = "G";
						$trial_sol[$map[$i]] = "C";
					}
					else {
						$trial_sol[$i] = "C";
						$trial_sol[$map[$i]] = "G";
					}
				}
				else {
					if ($bp_direction==0) {
						$trial_sol[$i] = "A";
						$trial_sol[$map[$i]] = "U";
					}
					else {
						$trial_sol[$i] = "U";
						$trial_sol[$map[$i]] = "A";
					}
				}
			}
		}
	}

	#####Compute Distances

	$new_trial_seq = join "",@trial_sol;	#Make a string to hold the sequence of the new trial

	&fold($new_trial_seq);
	my $new_fold = $fold;			  #The fold of puzzle solution
	&distance($new_fold, $targetstring);
	my $trial_dist = $distance;	  #Distance the puzzle solution is away from the target


	#Decide if we keep the new sequence

	#Process outputs

	my $print_output=0;

	$report_KL=1;
	&pattern_prevent;
	&countrepeats;  #search for repeated sequences


	if ($distance == 0 && ($KL_repeats<=$old_KL_repeats) && $GC_cont<$max_GC && $pattern_repeats==0 && $poly_repeats==0  && $restriction_sites==0 && $complement_zones==0 && $duplication_zones==0) {$keep_sequence=1;}    ##ADD PATTERN PREVENT CHECKING HERE.

	if ($keep_sequence==1){		   #Check if the new mutant is better than the parent
		for ( $i=0; $i<$strand_length; $i++){
			if (($best_sol[$i]) eq ($trial_sol[$i])) {
				$mut_map[$i] = "-";
			} else {
				$mut_map[$i] = "*";
			}
			$best_sol[$i]= $trial_sol[$i];		#move the trial solution to best_sol array
		}
		$mut_map_text = join "",@mut_map;
		$best_sol_seq = join "",@best_sol;

		$parent_dist=$trial_dist;			  #copy the trial distance to be the new parent_dist
		$best_dist=$trial_dist;
		@parent_fold = split(//,$fold);
		$parent_fold_string = $fold;

		&printer("\nIteration# $attempts\n");  #Report outputs
		&countGC();
		&printer("KL Repeats: $KL_repeats\n");
		&printer("$best_sol_seq\n\n");

	} else {
		&printer(" X \n");
		for ( $i=0; $i<$strand_length; $i++){
			$trial_sol[$i]=$best_sol[$i];  #reset the trial solution
		}
	}

	if ($distance==0 && $KL_repeats==0) {$happy = 1;}

}while ($happy==0);	#Cycle through sequence design again until the target is within the threshold distance from the fold.

&printer("\nInitial design of $name is complete!\n\n");

$full_KL_output =0;
$report_KL=1;
&analyzeKL;

@parent_KL_opt_target = @KL_opt_target;	 ##store the parent NTlist

##############################################################################################

&printer("\n\n__________________________________________________________\n");
&printer("Optimizing KL Sequences\n");	 #report progress

$happy=0;

do{

	$attempts += 1;

	&printer("\n(Iteration $attempts)");

	my $mutate_now = 0;
	my $parent_fold = "";
	my $coin_flip =0;
	my $rand_num_one =0;
	my $rand_num_two =0;
	my $old_KL_repeats=$KL_repeats;

	&printer(" [M$numKL_sites] ");

	for ( $i=0; $i<$strand_length; $i++){	 ##look for the repeated KLs and mutate them
		if ((($target[$i]eq"[") || ($target[$i]eq"]") )  && $constraint[$i]eq"N"){
			$rand_num_one = rand(100);
			$bp_direction= int(rand(2));
			$bp_type=(rand(100));
			if(rand(1)<(6*$parent_KL_opt_target[$i]/($numKL_sites+.01) ) ){    #~6 mutations per round.
				&printer("$i: ");
				if ($bp_type<$KL_GC) {
					if ($bp_direction==0) {
						$trial_sol[$i] = "G";
						$trial_sol[$map[$i]] = "C";
					}else {
						$trial_sol[$i] = "C";
						$trial_sol[$map[$i]] = "G";
					}
				}else {
					if ($bp_direction==0) {
						$trial_sol[$i] = "A";
						$trial_sol[$map[$i]] = "U";
					}else {
						$trial_sol[$i] = "U";
						$trial_sol[$map[$i]] = "A";
					}
				}
			}
		}
	}

	#####Compute Distances
	$new_trial_seq = join "",@trial_sol;	#Make a string to hold the sequence of the new trial

	$report_KL=0;
	&pattern_prevent;  ##search for KL
	&countrepeats;  #search for repeated sequences

	my $new_fold = " ";
	my $trial_dist = 1;


	my $unmasked_pattern = 0;
	my $masked_pattern =0;
	for($i=0;$i<$strand_length; $i++){
		if( ($constraint[$i]eq"N" || $constraint[$i]eq"S" || $constraint[$i]eq"K" ) && ( $repeat_map[$i]ne"-" ) ){ $unmasked_pattern+=1;  }
		if( ($constraint[$i]eq"A" || $constraint[$i]eq"U" || $constraint[$i]eq"C" || $constraint[$i]eq"G" ) && ( $repeat_map[$i]ne"-" ) ){ $masked_pattern+=1;  }
	}
	if($unmasked_pattern==0 && $masked_pattern>0){  ##set the patterns to zero once we verify that the only patterns are masked ones.
		$poly_repeats=0;
		$restriction_sites=0;
		$pattern_repeats =0;
		$complement_zones=0;
		$duplication_zones=0;
	}

	my $pattern_prevent_failure =0;
	if($KL_repeats==0 && $poly_repeats==0  && $restriction_sites==0 && $pattern_repeats ==0 && $complement_zones==0 && $duplication_zones==0){	  #only fold if the new sequences passes KL Repeats test
		&printer("\nFolding...\n");
		&fold($new_trial_seq);
		$new_fold = $fold;			   #The fold of puzzle solution
		&distance($new_fold, $targetstring);
		$trial_dist = $distance;   #Distance the puzzle solution is away from the target

		$last_pattern_repeats=$pattern_repeats;
		$last_poly_repeats=$poly_repeats;
		$last_restriction_sites = $restriction_sites;

	} else {
		$trial_dist = 1;  ##force a fail
		$pattern_prevent_failure=1;
	}

	my $print_output=0;

	my $kl_analyzed = 0;
	if($trial_dist ==0 && $KL_repeats==0 && $previously_tried==0){	 #only analyze the KLs if it folded right
		&printer("\nFolding succesful.\nAnalyzing KLs...\n\n");
		$report_KL=1;
		&analyzeKL;
		$kl_analyzed=1;
	}else{$kl_analyzed =0;}

	if ($KL_score<=$best_KL_score && $trial_dist == 0 && $KL_repeats==0 && $previously_tried==0 && $poly_repeats==0  && $restriction_sites==0 && $pattern_repeats ==0 && $complement_zones==0 && $duplication_zones==0){  #Check if the design performs better or not

		$best_KL_score=$KL_score;

		for ( $i=0; $i<$strand_length; $i++){
			$parent_KL_opt_target[$i] = $KL_opt_target[$i];
			if (($best_sol[$i]) eq ($trial_sol[$i])) {
				$mut_map[$i] = "-";
			} else {
				$mut_map[$i] = $trial_sol[$i];
			}
			$best_sol[$i]= $trial_sol[$i];		#move the trial solution to best_sol array
		}
		$mut_map_text = join "",@mut_map;
		$best_sol_seq = join "",@best_sol;
		$parent_dist=$trial_dist;			  #copy the trial distance to be the new parent_dist
		$best_dist=$trial_dist;
		@parent_fold = split(//,$fold);
		$parent_fold_string = $fold;
		&printer("\nIteration# $attempts\n");  #Report outputs
		&countGC();
		&printer("KL Repeats: $KL_repeats\n");
		&printer("KL Score: $KL_score\n");
		$KL_opt_target_text = join "",@parent_KL_opt_target;
		&printer("\nKL Target Mask:\n$KL_opt_target_text\n\n");
		&printer("$best_sol_seq\n\n");
		&printer("$mut_map_text\n\n");

		if ($generate_previews==1){&preview;}

	} else {

		&printer("\nRound $attempts Failed:   ");
		##print error summary here:
		if($KL_repeats>0){&printer("   $KL_repeats KL repeats were found\n");}
		if($previously_tried>0){&printer( "   The Sequence was previously tested\n");}
		if($KL_score>$best_KL_score && $kl_analyzed==1){&printer("   The KL score $KL_score is greater than the best score $best_KL_score\n")}
		if($pattern_prevent_failure==1){&printer("   Pattern Prevent Failed.\n");}

		for($i=0; $i<$strand_length; $i++){
			$KL_opt_target[$i] = $parent_KL_opt_target[$i];
			$trial_sol[$i]=$best_sol[$i];   #reset the trial solution
		}
	}
	if ($distance==0 && $KL_score==0) {$happy = 1;}
}while ($happy==0);	#Cycle through sequence design again until the target is within the threshold distance from the fold.







&victory;
&printer("\nOptimization of design $name is complete!\n\n");

&printer("\n\n\nIteration# $attempts\n");	#Report outputs
&countGC();
&printer("KL Repeats: $KL_repeats\n");
&printer("KL Score: $KL_score\n");
&printer("$best_sol_seq\n\n");

if ($generate_previews==1){&preview;}

$full_KL_output =1;
$report_KL=0;
&analyzeKL;

&countrepeats;

# qx(perl trace_analysis.pl pattern.txt seq.txt > $output_file_path$name\_trace.txt);

&export;

close $output_spool;

##############################################################################################
# ANALYSIS Subroutines
################################################

sub countGC {#Count GCs
   $pairs=0;
   $GC_pairs=0;
   $AU_pairs=0;
   $GU_pairs=0;
   $KL_pairs=0;
   $KL_GC_pairs=0;

   for ( $i=0; $i<$strand_length; $i++){
	 if ($target[$i]eq"("){
		$pairs += 1;
		if ((($best_sol[$i]eq"C") && ($best_sol[$map[$i]]eq"G")) || (($best_sol[$i]eq"G") && ($best_sol[$map[$i]]eq"C"))) {
		  $GC_pairs+=1;
		}
		if ((($best_sol[$i]eq"A") && ($best_sol[$map[$i]]eq"U")) || (($best_sol[$i]eq"U") && ($best_sol[$map[$i]]eq"A"))) {
		  $AU_pairs+=1;
		}
		if ((($best_sol[$i]eq"U") && ($best_sol[$map[$i]]eq"G")) || (($best_sol[$i]eq"G") && ($best_sol[$map[$i]]eq"U"))) {
		  $GU_pairs+=1;
		}
	 }
	 elsif ($init_target[$i]eq"["){	 ##only count KLs, not '{' hinted bps.
		$KL_pairs += 1;
		if ((($best_sol[$i]eq"C") && ($best_sol[$map[$i]]eq"G")) || (($best_sol[$i]eq"G") && ($best_sol[$map[$i]]eq"C"))) {
		  $KL_GC_pairs+=1;
		}
	 }
   }

   if ($pairs>0){$GC_cont=int(1000*($GC_pairs)/($pairs))/10;}else{$GC_cont=-1;}
   if($KL_pairs>0){$KL_GC_cont= int(1000*$KL_GC_pairs/$KL_pairs)/10;} else {$KL_GC_cont=0;}

   &printer("\n GC Content: $GC_cont	TotalPairs: $pairs\n GC pairs: $GC_pairs	 AU pairs: $AU_pairs	 GU pairs: $GU_pairs	KL: $KL_GC_cont %GC\n");
}

sub countrepeats {
	for ($i=0; $i<$strand_length; $i++){$repeat_map[$i]="-";}  ##zero the map

    my $trial_sol_seq = join "",@trial_sol;
    my $pattern_sense = "AAAAAAA";
    my $pattern_antisense = "UUUUUUU";



	for ( $i=0; $i<($strand_length-$complement_window); $i++){		#search for antisense matches

		$pattern_sense = substr($trial_sol_seq, $i, $complement_window);  #grab the substring pattern here
		$pattern_antisense = "";  #reset the pattern here
		for (my $j=0; $j<$complement_window+1; $j++){ #generate the complement sequence
			if (substr($pattern_sense,$complement_window-$j,1) eq "A" ){$pattern_antisense=$pattern_antisense."U";}
			if (substr($pattern_sense,$complement_window-$j,1) eq "U" ){$pattern_antisense=$pattern_antisense."A";}
			if (substr($pattern_sense,$complement_window-$j,1) eq "C" ){$pattern_antisense=$pattern_antisense."G";}
			if (substr($pattern_sense,$complement_window-$j,1) eq "G" ){$pattern_antisense=$pattern_antisense."C";}
		}
		##  print "$pattern_sense - $pattern_antisense\n";   #for bugchecking

		my $test_result = index($trial_sol_seq, $pattern_antisense, 0);

		for (my $j=$i; $j<$strand_length-$complement_window; $j++){
			$test_result = index($trial_sol_seq, $pattern_antisense, $j);   #search for antisense matches
			if ($test_result != -1){
				for (my $k=0; $k<$complement_window; $k++){
					$repeat_map[$i+$k]="P";
					$repeat_map[$test_result+$k]="P";
				}
			}
		}

	}

	$complement_zones=0;
	for ($i=0; $i<$strand_length; $i++){
		if ($repeat_map[$i]eq"P" || $repeat_map[$i]eq"p"){$complement_zones+=1;}
	}


	for ( $i=0; $i<($strand_length-$duplicate_window); $i++){		#search for sense duplications

		$pattern_sense = substr($trial_sol_seq, $i, $duplicate_window);  #grab the substring pattern here
		my $test_result = index($trial_sol_seq, $pattern_sense, 0);

		for (my $j=$i; $j<$strand_length-$duplicate_window; $j++){
			$test_result = index($trial_sol_seq, $pattern_sense, $j);    	#search for sense duplications
			if ($test_result != -1 && $test_result > $j ){
				for (my $k=0; $k<$duplicate_window; $k++){ ##bugfix try $k=1 rather than 0..
					$repeat_map[$i+$k]="D";
					$repeat_map[$test_result+$k]="D";
				}
			}
		}

	}

	$duplication_zones=0;
	for ($i=0; $i<$strand_length; $i++){
		if ($repeat_map[$i]eq"D"){$duplication_zones+=1;}
	}


	$pattern_repeats=0;
	for ( $i=0; $i<$strand_length-7; $i++){	   #count up stretches of A/U or G/C that are 8nts
		if ( ($trial_sol[$i]eq"A" || $trial_sol[$i]eq"U") &&
		 ($trial_sol[$i+1]eq"A" || $trial_sol[$i+1]eq"U") &&
		 ($trial_sol[$i+2]eq"A" || $trial_sol[$i+2]eq"U") &&
		 ($trial_sol[$i+3]eq"A" || $trial_sol[$i+3]eq"U") &&
		 ($trial_sol[$i+4]eq"A" || $trial_sol[$i+4]eq"U") &&
		 ($trial_sol[$i+5]eq"A" || $trial_sol[$i+5]eq"U") &&
		 ($trial_sol[$i+6]eq"A" || $trial_sol[$i+6]eq"U") &&
		 ($trial_sol[$i+7]eq"A" || $trial_sol[$i+7]eq"U") ) {  #all stretch of GC_rich are targeted
			$pattern_repeats+=1;
			$repeat_map[$i]="W"; $repeat_map[$i+1]="W";
			$repeat_map[$i+2]="W"; $repeat_map[$i+3]="W";
			$repeat_map[$i+4]="W"; $repeat_map[$i+5]="W";
			$repeat_map[$i+6]="W"; $repeat_map[$i+7]="W";
		}
		if ( ($trial_sol[$i]eq"C" || $trial_sol[$i]eq"G") &&
		 ($trial_sol[$i+1]eq"C" || $trial_sol[$i+1]eq"G") &&
		 ($trial_sol[$i+2]eq"C" || $trial_sol[$i+2]eq"G") &&
		 ($trial_sol[$i+3]eq"C" || $trial_sol[$i+3]eq"G") &&
		 ($trial_sol[$i+4]eq"C" || $trial_sol[$i+4]eq"G") &&
		 ($trial_sol[$i+5]eq"C" || $trial_sol[$i+5]eq"G") &&
		 ($trial_sol[$i+6]eq"C" || $trial_sol[$i+6]eq"G") &&
		 ($trial_sol[$i+7]eq"C" || $trial_sol[$i+7]eq"G") ) {  #all stretch of GC_rich are targeted
		 	if ($constraint[$i]eq"S" && $constraint[$i+1]eq"S" && $constraint[$i+2]eq"S" && $constraint[$i+3]eq"S" &&
		 		$constraint[$i+4]eq"S" && $constraint[$i+5]eq"S" && $constraint[$i+6]eq"S" && $constraint[$i+7]eq"S" ){  #only allow 8S in a row in the rare case that they have been user-specified
				&printer ("\nLocked contraint of 8S found at position $i\n");
			}else {$pattern_repeats+=1;}
			$repeat_map[$i]="S"; $repeat_map[$i+1]="S";
			$repeat_map[$i+2]="S"; $repeat_map[$i+3]="S";
			$repeat_map[$i+4]="S"; $repeat_map[$i+5]="S";
			$repeat_map[$i+6]="S"; $repeat_map[$i+7]="S";
		}
	}

	$poly_repeats=0;

		for ( $i=0; $i<$strand_length-6; $i++){
		if ( ($trial_sol[$i]eq$trial_sol[$i+1]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+2]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+3]) &&
			 ($trial_sol[$i]eq$trial_sol[$i+4])) {  #all stretch of 5G,5C,%U,5A strand are counted
			$poly_repeats+=1;
			$repeat_map[$i]=$trial_sol[$i]; $repeat_map[$i+1]=$trial_sol[$i+1];
			$repeat_map[$i+2]=$trial_sol[$i+2]; $repeat_map[$i+3]=$trial_sol[$i+3];
			$repeat_map[$i+4]=$trial_sol[$i+4];
		}
	}


	$restriction_sites =0;

	for ( $i=0; $i<$strand_length-6; $i++){
		if 	(	( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"G") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## Bsa1 GGUCUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"G") && ($trial_sol[$i+3]eq"A") && ($trial_sol[$i+4]eq"C") && ($trial_sol[$i+5]eq"C") ) ||    ## Bsa1 GAGACC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"A") && ($trial_sol[$i+3]eq"G") && ($trial_sol[$i+4]eq"A") && ($trial_sol[$i+5]eq"C") ) ||    ## BBs1 GAAGAC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"U") && ($trial_sol[$i+2]eq"C") && ($trial_sol[$i+3]eq"U") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## BBs1 GUCUUC
				( ($trial_sol[$i]eq"C") && ($trial_sol[$i+1]eq"G") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"C") ) ||    ## BsMB1 CGUCUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"G") && ($trial_sol[$i+3]eq"A") && ($trial_sol[$i+4]eq"C") && ($trial_sol[$i+5]eq"G") ) ||    ## BsMB1 GAGACG
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"C") && ($trial_sol[$i+2]eq"U") && ($trial_sol[$i+3]eq"C") && ($trial_sol[$i+4]eq"U") && ($trial_sol[$i+5]eq"U") && ($trial_sol[$i+6]eq"C")) ||   ## Sap1 GCUCUUC
				( ($trial_sol[$i]eq"G") && ($trial_sol[$i+1]eq"A") && ($trial_sol[$i+2]eq"A") && ($trial_sol[$i+3]eq"G") && ($trial_sol[$i+4]eq"A") && ($trial_sol[$i+5]eq"G") && ($trial_sol[$i+6]eq"C")) ||   ## Sap1 GAAGAGC
				( ($trial_sol[$i]eq"A") && ($trial_sol[$i+1]eq"U") && ($trial_sol[$i+2]eq"C") && ($trial_sol[$i+3]eq"U") && ($trial_sol[$i+4]eq"G") && ($trial_sol[$i+5]eq"U") && ($trial_sol[$i+6]eq"U")) ) {    ## PTH  AUCUGUU

			$restriction_sites+=1;
			$repeat_map[$i]="X"; $repeat_map[$i+1]="X";
			$repeat_map[$i+2]="X"; $repeat_map[$i+3]="X";
			$repeat_map[$i+4]="X"; $repeat_map[$i+5]="X";
		}
	}
}



sub pattern_prevent {
$KL_repeats=0;
for ( $i=0; $i<$strand_length; $i++){
   $pattern_zones[$i] = "-";  ##clear out the pattern array
}

for ( $i=0; $i<$strand_length-10; $i++){	#seek any repeated KL sequence and mark them.
	for (my $j=($i+1); $j<$strand_length-9; $j++){
		if(($trial_sol[$j]eq$trial_sol[$i]) &&   ##look for 180KL signature
		($trial_sol[$j+1]eq$trial_sol[$i+1]) &&
		($trial_sol[$j+2]eq$trial_sol[$i+2]) &&
		($trial_sol[$j+3]eq$trial_sol[$i+3]) &&
		($trial_sol[$j+4]eq$trial_sol[$i+4]) &&
		($trial_sol[$j+5]eq$trial_sol[$i+5]) &&
		($trial_sol[$j+6]eq$trial_sol[$i+6]) &&
		($trial_sol[$j+7]eq$trial_sol[$i+7]) &&
		($targetnopk[$j]eq".") &&
		($targetnopk[$j+1]eq".") &&
		($targetnopk[$j+2]eq".") &&
		($targetnopk[$j+3]eq".") &&
		($targetnopk[$j+4]eq".") &&
		($targetnopk[$j+5]eq".") &&
		($targetnopk[$j+6]eq".") &&
		($targetnopk[$j+7]eq".") &&
		($targetnopk[$i]eq".") &&
		($targetnopk[$i+1]eq".") &&
		($targetnopk[$i+2]eq".") &&
		($targetnopk[$i+3]eq".") &&
		($targetnopk[$i+4]eq".") &&
		($targetnopk[$i+5]eq".") &&
		($targetnopk[$i+6]eq".") &&
		($targetnopk[$i+7]eq".") &&
               ($constraint[$j+4]eq"N")){
			$pattern_zones[$i]=$trial_sol[$i];	 $pattern_zones[$j]=$trial_sol[$i];
			$pattern_zones[$i+1]=$trial_sol[$i+1]; $pattern_zones[$j+1]=$trial_sol[$i+1];
			$pattern_zones[$i+2]=$trial_sol[$i+2]; $pattern_zones[$j+2]=$trial_sol[$i+2];
			$pattern_zones[$i+3]=$trial_sol[$i+3]; $pattern_zones[$j+3]=$trial_sol[$i+3];
			$pattern_zones[$i+4]=$trial_sol[$i+4]; $pattern_zones[$j+4]=$trial_sol[$i+4];
			$pattern_zones[$i+5]=$trial_sol[$i+5]; $pattern_zones[$j+5]=$trial_sol[$i+5];
			$pattern_zones[$i+6]=$trial_sol[$i+6]; $pattern_zones[$j+6]=$trial_sol[$i+6];
			$pattern_zones[$i+7]=$trial_sol[$i+7]; $pattern_zones[$j+7]=$trial_sol[$i+7];
			$KL_repeats += 1;
		}
		elsif (($trial_sol[$j+1]eq$trial_sol[$i+1]) &&   ##look for 120KL signature
		($trial_sol[$j+2]eq$trial_sol[$i+2]) &&
		($trial_sol[$j+3]eq$trial_sol[$i+3]) &&
		($trial_sol[$j+4]eq$trial_sol[$i+4]) &&
		($trial_sol[$j+5]eq$trial_sol[$i+5]) &&
		($trial_sol[$j+6]eq$trial_sol[$i+6]) &&
		($trial_sol[$j+7]eq$trial_sol[$i+7]) &&
		($targetnopk[$j]eq"(") &&
		($targetnopk[$j+1]eq".") &&
		($targetnopk[$j+2]eq".") &&
		($targetnopk[$j+3]eq".") &&
		($targetnopk[$j+4]eq".") &&
		($targetnopk[$j+5]eq".") &&
		($targetnopk[$j+6]eq".") &&
		($targetnopk[$j+7]eq".") &&
		($targetnopk[$j+8]eq")") &&
		($targetnopk[$i]eq"(") &&
		($targetnopk[$i+1]eq".") &&
		($targetnopk[$i+2]eq".") &&
		($targetnopk[$i+3]eq".") &&
		($targetnopk[$i+4]eq".") &&
		($targetnopk[$i+5]eq".") &&
		($targetnopk[$i+6]eq".") &&
		($targetnopk[$i+7]eq".") &&
		($targetnopk[$i+8]eq")") &&
               ($constraint[$j+4]eq"N") ){
			$pattern_zones[$i+1]=$trial_sol[$i+1]; $pattern_zones[$j+1]=$trial_sol[$i+1];
			$pattern_zones[$i+2]=$trial_sol[$i+2]; $pattern_zones[$j+2]=$trial_sol[$i+2];
			$pattern_zones[$i+3]=$trial_sol[$i+3]; $pattern_zones[$j+3]=$trial_sol[$i+3];
			$pattern_zones[$i+4]=$trial_sol[$i+4]; $pattern_zones[$j+4]=$trial_sol[$i+4];
			$pattern_zones[$i+5]=$trial_sol[$i+5]; $pattern_zones[$j+5]=$trial_sol[$i+5];
			$pattern_zones[$i+6]=$trial_sol[$i+6]; $pattern_zones[$j+6]=$trial_sol[$i+6];
			$pattern_zones[$i+7]=$trial_sol[$i+7]; $pattern_zones[$j+7]=$trial_sol[$i+7];
			$KL_repeats += 1;
		}
	}
}

my $pattern_zones_text = join "",@pattern_zones;
if ($report_KL==1){&printer("\n$KL_repeats repeated KL sequences:\n$pattern_zones_text\n\n");}
else {&printer("\n$KL_repeats KL repeats found.\n\n");}

}

sub duplex { #outputs to $distance
	my $duplex_num_counter=0;
	my $dont_fold = 0;

	foreach my $duplex_name (@duplex_catalog){
		if(("$_[0].$_[1]" eq $duplex_name) ||
		   ("$_[1].$_[0]" eq $duplex_name) ){
			$KLenergy=$duplex_out_catalog[$duplex_num_counter];
			$parens=$duplex_parens_catalog[$duplex_num_counter];
			$dont_fold=1;
		}
		$duplex_num_counter ++;
	}

	if($dont_fold==0){
		# open(my $duin, '>', 'KLtest.txt');
		# 	print $duin "$_[0]\n";
		# 	print $duin "$_[1]\n";
		# close $duin;

		# $outputs = qx(rnaduplex < KLtest.txt);
		$outputs = qx(python viennarna_funcs.py duplex "$_[0]" "$_[1]");
		@results = split(' ', $outputs);
		$parens = $results[0];
		$KLenergy = $results[1];
		chomp($KLenergy);

		# $counter = 0; $KLenergy = "";

		# foreach(@results){
		# 	if($counter == 0){ $parens = "$_"; }
		# 	if($counter == 4){ $KLenergy = "$_"; }
		# 	if($counter == 5){ $KLenergy = "$KLenergy$_"; }
		# 	$counter ++;
		# }

		# my $KLenergy_clean = $KLenergy =~ qr/(-||' ')(\d+).(\d+)/; ##parse out the number from the (-X.XX) in parenthesis
		# $KLenergy=$&;  ##store the regex part back into the variable

		$duplex_catalog[$unique_duplex_id]="$_[0].$_[1]";
		$duplex_out_catalog[$unique_duplex_id]=$KLenergy;
		$duplex_parens_catalog[$unique_duplex_id]=$parens;

		$unique_duplex_id ++;
	}

}

sub analyzeKL {

	my $threshold = $KL_min_delta_G;

	####reset the global variables
	$numKL=0;
	$numKL_sites=0;
	$KL_notcounted =0;
	$KL_score=0;
	@Alist = ();  #array of KLA
	@Blist = ();  #array of KLB
	@Elist = ();  #array of Energies
	@NTlist = ();
	for ( $i=0; $i<$strand_length; $i++){  #fill the KL-target array with -
		$KL_opt_target[$i]=0;
	}


	if($report_KL==1){ &printer("Mapping Kissing Loops:\n");}

	for ( $i=0; $i<$strand_length-8; $i++){	  #Look for ".]]]]]]." and extract the KLs as $KL_A and $KL_B
		if (  ($target[$i]eq(".") || $target[$i]eq(")") ) &&   #added ")" to identify bKLs as KL motifs
			$target[$i+1]eq("]") &&  ##search these in order of closing
			$target[$i+2]eq("]") &&
			$target[$i+3]eq("]") &&
			$target[$i+4]eq("]") &&
			$target[$i+5]eq("]") &&
			$target[$i+6]eq("]") &&
			$target[$i+7]eq(".")){

			my $BKLflag=0;
			if($target[$i]eq(")")){
				$BKLflag=1;
				$KL_B="$trial_sol[$i+1]$trial_sol[$i+2]$trial_sol[$i+3]$trial_sol[$i+4]$trial_sol[$i+5]$trial_sol[$i+6]$trial_sol[$i+7]";
			} else {
				$KL_B="$trial_sol[$i]$trial_sol[$i+1]$trial_sol[$i+2]$trial_sol[$i+3]$trial_sol[$i+4]$trial_sol[$i+5]$trial_sol[$i+6]$trial_sol[$i+7]";
			}

			if($target[$map[$i+1]-1]eq("(")){
			$BKLflag=2;
			$KL_A="$trial_sol[$map[$i+6]-1]$trial_sol[$map[$i+6]]$trial_sol[$map[$i+5]]$trial_sol[$map[$i+4]]$trial_sol[$map[$i+3]]$trial_sol[$map[$i+2]]$trial_sol[$map[$i+1]]";
			} else {
			$KL_A="$trial_sol[$map[$i+6]-1]$trial_sol[$map[$i+6]]$trial_sol[$map[$i+5]]$trial_sol[$map[$i+4]]$trial_sol[$map[$i+3]]$trial_sol[$map[$i+2]]$trial_sol[$map[$i+1]]$trial_sol[$map[$i+1]+1]";
			}



			$numKL++;
			if($numKL<10){&printer(" ");}#add a space to make it look nice
			if($report_KL==1 && $BKLflag==0){&printer("$numKL,180 KL,");}
			if($report_KL==1 && $BKLflag==1){&printer("$numKL,KL-BKL,");}
			if($report_KL==1 && $BKLflag==2){&printer("$numKL,BKL-KL,");}


			&duplex($KL_A, $KL_B);
			if($report_KL==1){&printer("$KL_A,$KL_B,$parens,$KLenergy");}

			$Alist[$numKL]=$KL_A;
			$Blist[$numKL]=$KL_B;
			$NTlist[$numKL]=$i;
			$Elist[$numKL]=$KLenergy;

			if ($constraint[$i+1]eq"N" && $constraint[$i+2]eq"N" && $constraint[$i+3]eq"N" &&   #ignore KL calculation for KLs with ANY constraints in them
			$constraint[$i+4]eq"N" && $constraint[$i+5]eq"N" && $constraint[$i+6]eq"N"){
			$KL_specified[$numKL]="N";}else{$KL_specified[$numKL]="Y"; $KL_notcounted++; if($report_KL==1){&printer(" (ignored)");}}

			if($report_KL==1){&printer("\n");}

			$KL_type[$numKL]="A";

			if((($KLenergy> $MinKL) || ($KLenergy< $MaxKL)) && $KL_specified[$numKL]eq"N"){    ##Scoring function for 180KLs
			   $KL_score+=2;
			   $numKL_sites += 6;
			   $KL_opt_target[$i+1]+=1; $KL_opt_target[$i+2]+=1;
			   $KL_opt_target[$i+3]+=1; $KL_opt_target[$i+4]+=1;
			   $KL_opt_target[$i+5]+=1; $KL_opt_target[$i+6]+=1;
			}

		}
		elsif ($target[$i]eq("[") &&
			$target[$i+1]eq("[") &&
			$target[$i+2]eq("[") &&
			$target[$i+3]eq("[") &&
			$target[$i+4]eq("[") &&
			$target[$i+5]eq("[") &&
			$target[$i+6]eq("[")  ){

			$KL_A="$trial_sol[$i]$trial_sol[$i+1]$trial_sol[$i+2]$trial_sol[$i+3]$trial_sol[$i+4]$trial_sol[$i+5]$trial_sol[$i+6]";

			$KL_B="$trial_sol[$map[$i+6]]$trial_sol[$map[$i+5]]$trial_sol[$map[$i+4]]$trial_sol[$map[$i+3]]$trial_sol[$map[$i+2]]$trial_sol[$map[$i+1]]$trial_sol[$map[$i]]";
			$numKL++;



			if($numKL<10){&printer(" ");}#add a space to make it look nice
			if($report_KL==1){&printer("$numKL,120 KL,");}
			&duplex($KL_A, $KL_B);
			if($report_KL==1){&printer("$KL_A,$KL_B,$parens,$KLenergy");}

			if ($constraint[$i]eq"N" && $constraint[$i+1]eq"N" && $constraint[$i+2]eq"N" && $constraint[$i+3]eq"N" &&   #ignore KL calculation for KLs with ANY constraints in them
			$constraint[$i+4]eq"N" && $constraint[$i+5]eq"N" && $constraint[$i+6]eq"N"){
			$KL_specified[$numKL]="N";}else{$KL_specified[$numKL]="Y"; $KL_notcounted++; if($report_KL==1){&printer(" (ignored)");}}

			if($report_KL==1){&printer("\n");}

			$Alist[$numKL]=$KL_A;
			$Blist[$numKL]=$KL_B;
			$NTlist[$numKL]=$i;
			$Elist[$numKL]=$KLenergy;
			$KL_type[$numKL]="B";


			if( (($KLenergy> $MinKL) || ($KLenergy< $MaxKL)) && $KL_specified[$numKL]eq"N"){  ##Scoring function for 120KLs
				$KL_score+=2;
				$numKL_sites += 7;
				$KL_opt_target[$i]+=1; $KL_opt_target[$i+1]+=1;
				$KL_opt_target[$i+2]+=1; $KL_opt_target[$i+3]+=1;
				$KL_opt_target[$i+4]+=1; $KL_opt_target[$i+5]+=1;
				$KL_opt_target[$i+6]+=1;
			}
		}

	}
	my $calculation = $KL_score/2;
	&printer("\n\n$calculation KLs with energy lower than $MinKL kcal or greater than $MaxKL kcal were counted.\n\n");
	if ($KL_notcounted>0){&printer("$KL_notcounted KLs were not counted towards penalties because they were user specified.\n\n");}
	$non_cognates = "";

		if($report_KL==1){&printer("\n\nMapping Non-Cognate Pairings (non-cognate at least as strong as on-target, or stronger than $KL_min_delta_G kcal):\n");}
		if($full_KL_output==1){&printer("\n\nMapping Non-Cognate Pairings (all non-cognate at least as strong as $KL_min_delta_G kcal):\n");}
		$non_cognates=$non_cognates."\n\nMapping Non-Cognate Pairings (non-cognate at least as strong as on-target, or stronger than $KL_min_delta_G kcal):\n";

		for ($i=1; $i<($numKL); $i++){
			for (my $j=($i+1); $j<($numKL+1); $j++){

				&duplex($Alist[$i], $Alist[$j]);
				if( ($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G) ){
					if($KL_specified[$i]eq"N" || $KL_specified[$j]eq"N"){$KL_score++;}

					if($KL_type[$i]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$i]+1]+=1; $KL_opt_target[$NTlist[$i]+2]+=1;
						$KL_opt_target[$NTlist[$i]+3]+=1; $KL_opt_target[$NTlist[$i]+4]+=1;
						$KL_opt_target[$NTlist[$i]+5]+=1; $KL_opt_target[$NTlist[$i]+6]+=1;
					}
					elsif($KL_type[$i]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$i]]+=1; $KL_opt_target[$NTlist[$i]+1]+=1;
						$KL_opt_target[$NTlist[$i]+2]+=1; $KL_opt_target[$NTlist[$i]+3]+=1;
						$KL_opt_target[$NTlist[$i]+4]+=1; $KL_opt_target[$NTlist[$i]+5]+=1;
						$KL_opt_target[$NTlist[$i]+6]+=1;
					}

					if($KL_type[$j]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$j]+1]+=1; $KL_opt_target[$NTlist[$j]+2]+=1;
						$KL_opt_target[$NTlist[$j]+3]+=1; $KL_opt_target[$NTlist[$j]+4]+=1;
						$KL_opt_target[$NTlist[$j]+5]+=1; $KL_opt_target[$NTlist[$j]+6]+=1;
					}
					elsif($KL_type[$j]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$j]]+=1; $KL_opt_target[$NTlist[$j]+1]+=1;
						$KL_opt_target[$NTlist[$j]+2]+=1; $KL_opt_target[$NTlist[$j]+3]+=1;
						$KL_opt_target[$NTlist[$j]+4]+=1; $KL_opt_target[$NTlist[$j]+5]+=1;
						$KL_opt_target[$NTlist[$j]+6]+=1;
					}
				}
				if(($full_KL_output==1 && $KLenergy<$threshold) || ($report_KL==1 && ($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G)) ){
				   &printer("$i(A) + $j(A),$Alist[$i],$Alist[$j],$parens,$KLenergy\n");
				   $non_cognates=$non_cognates."$i(A) + $j(A),$Alist[$i],$Alist[$j],$parens,$KLenergy\n";
				}


				&duplex($Alist[$i], $Blist[$j]);
				if(($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G ) ){
					if($KL_specified[$i]eq"N" || $KL_specified[$j]eq"N"){$KL_score++;}

					if($KL_type[$i]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$i]+1]+=1; $KL_opt_target[$NTlist[$i]+2]+=1;
						$KL_opt_target[$NTlist[$i]+3]+=1; $KL_opt_target[$NTlist[$i]+4]+=1;
						$KL_opt_target[$NTlist[$i]+5]+=1; $KL_opt_target[$NTlist[$i]+6]+=1;
					}
					elsif($KL_type[$i]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$i]]+=1; $KL_opt_target[$NTlist[$i]+1]+=1;
						$KL_opt_target[$NTlist[$i]+2]+=1; $KL_opt_target[$NTlist[$i]+3]+=1;
						$KL_opt_target[$NTlist[$i]+4]+=1; $KL_opt_target[$NTlist[$i]+5]+=1;
						$KL_opt_target[$NTlist[$i]+6]+=1;
					}

					if($KL_type[$j]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$j]+1]+=1; $KL_opt_target[$NTlist[$j]+2]+=1;
						$KL_opt_target[$NTlist[$j]+3]+=1; $KL_opt_target[$NTlist[$j]+4]+=1;
						$KL_opt_target[$NTlist[$j]+5]+=1; $KL_opt_target[$NTlist[$j]+6]+=1;
					}
					elsif($KL_type[$j]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$j]]+=1; $KL_opt_target[$NTlist[$j]+1]+=1;
						$KL_opt_target[$NTlist[$j]+2]+=1; $KL_opt_target[$NTlist[$j]+3]+=1;
						$KL_opt_target[$NTlist[$j]+4]+=1; $KL_opt_target[$NTlist[$j]+5]+=1;
						$KL_opt_target[$NTlist[$j]+6]+=1;
					}
				}
				if(($full_KL_output==1 && $KLenergy<$threshold) || ($report_KL==1 && ($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G)) ){
				   &printer("$i(A) + $j(B),$Alist[$i],$Blist[$j],$parens,$KLenergy\n");
				   $non_cognates=$non_cognates."$i(A) + $j(B),$Alist[$i],$Blist[$j],$parens,$KLenergy\n";
				}


				&duplex($Blist[$i], $Alist[$j]);
				if(($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G ) ){
					if($KL_specified[$i]eq"N" || $KL_specified[$j]eq"N"){$KL_score++;}

					if($KL_type[$i]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$i]+1]+=1; $KL_opt_target[$NTlist[$i]+2]+=1;
						$KL_opt_target[$NTlist[$i]+3]+=1; $KL_opt_target[$NTlist[$i]+4]+=1;
						$KL_opt_target[$NTlist[$i]+5]+=1; $KL_opt_target[$NTlist[$i]+6]+=1;
					}
					elsif($KL_type[$i]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$i]]+=1; $KL_opt_target[$NTlist[$i]+1]+=1;
						$KL_opt_target[$NTlist[$i]+2]+=1; $KL_opt_target[$NTlist[$i]+3]+=1;
						$KL_opt_target[$NTlist[$i]+4]+=1; $KL_opt_target[$NTlist[$i]+5]+=1;
						$KL_opt_target[$NTlist[$i]+6]+=1;
					}

					if($KL_type[$j]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$j]+1]+=1; $KL_opt_target[$NTlist[$j]+2]+=1;
						$KL_opt_target[$NTlist[$j]+3]+=1; $KL_opt_target[$NTlist[$j]+4]+=1;
						$KL_opt_target[$NTlist[$j]+5]+=1; $KL_opt_target[$NTlist[$j]+6]+=1;
					}
					elsif($KL_type[$j]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$j]]+=1; $KL_opt_target[$NTlist[$j]+1]+=1;
						$KL_opt_target[$NTlist[$j]+2]+=1; $KL_opt_target[$NTlist[$j]+3]+=1;
						$KL_opt_target[$NTlist[$j]+4]+=1; $KL_opt_target[$NTlist[$j]+5]+=1;
						$KL_opt_target[$NTlist[$j]+6]+=1;
					}
				}
				if(($full_KL_output==1 && $KLenergy<$threshold) || ($report_KL==1 && ($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G)) ){
				   &printer("$i(B) + $j(A),$Blist[$i],$Alist[$j],$parens,$KLenergy\n");
				   $non_cognates=$non_cognates."$i(B) + $j(A),$Blist[$i],$Alist[$j],$parens,$KLenergy\n";
				}

				&duplex($Blist[$i], $Blist[$j]);
				if(($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G)){
					if($KL_specified[$i]eq"N" || $KL_specified[$j]eq"N"){$KL_score++;}

					if($KL_type[$i]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$i]+1]+=1; $KL_opt_target[$NTlist[$i]+2]+=1;
						$KL_opt_target[$NTlist[$i]+3]+=1; $KL_opt_target[$NTlist[$i]+4]+=1;
						$KL_opt_target[$NTlist[$i]+5]+=1; $KL_opt_target[$NTlist[$i]+6]+=1;
					}
					elsif($KL_type[$i]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$i]]+=1; $KL_opt_target[$NTlist[$i]+1]+=1;
						$KL_opt_target[$NTlist[$i]+2]+=1; $KL_opt_target[$NTlist[$i]+3]+=1;
						$KL_opt_target[$NTlist[$i]+4]+=1; $KL_opt_target[$NTlist[$i]+5]+=1;
						$KL_opt_target[$NTlist[$i]+6]+=1;
					}

					if($KL_type[$j]eq"A") {
						$numKL_sites += 6;
						$KL_opt_target[$NTlist[$j]+1]+=1; $KL_opt_target[$NTlist[$j]+2]+=1;
						$KL_opt_target[$NTlist[$j]+3]+=1; $KL_opt_target[$NTlist[$j]+4]+=1;
						$KL_opt_target[$NTlist[$j]+5]+=1; $KL_opt_target[$NTlist[$j]+6]+=1;
					}
					elsif($KL_type[$j]eq"B") {
						$numKL_sites += 7;
						$KL_opt_target[$NTlist[$j]]+=1; $KL_opt_target[$NTlist[$j]+1]+=1;
						$KL_opt_target[$NTlist[$j]+2]+=1; $KL_opt_target[$NTlist[$j]+3]+=1;
						$KL_opt_target[$NTlist[$j]+4]+=1; $KL_opt_target[$NTlist[$j]+5]+=1;
						$KL_opt_target[$NTlist[$j]+6]+=1;
					}
				}
				if(($full_KL_output==1 && $KLenergy<$threshold) || ($report_KL==1 && ($KLenergy<=$Elist[$i]) || ($KLenergy<=$Elist[$j]) || ($KLenergy<=$KL_min_delta_G)) ){
				   &printer("$i(B) + $j(B),$Blist[$i],$Blist[$j],$parens,$KLenergy\n");
					$non_cognates=$non_cognates."$i(B) + $j(B),$Blist[$i],$Blist[$j],$parens,$KLenergy\n";
				}
			}
		}
		&printer("\nFinal KL Score = $KL_score\n");
}



sub export {

		# my $ensemblediversity = qx(rnafold --noPS -p < seq.txt);
		my $ensemblediversity = qx(python viennarna_funcs.py diversity "$full_seq");
		chomp($ensemblediversity);

		# my $ensemblediversity_clean = $ensemblediversity =~ qr/;\s*(.+)$/; ##parse out everything after the semicolon
		# $ensemblediversity=$&;  ##store the regex part back into the variable

		if ($failed_initial_design == 1){
			&printer(" WARNING: Misfolding within the locked-sequence required altering the target structure. \n");
			&printer("          Sequence design failed for the inputted target structure.\n\n");
		}


		&printer("Exporting design to $name\_design\n\n");
		open(my $file_out, '>', "$output_file_path$name\_design.txt");
		print $file_out "$name\n";
		print $file_out "$targetstringpk\n";
		print $file_out "$best_sol_seq\n\n";

		print $file_out " GC Content: $GC_cont	TotalPairs: $pairs\n GC pairs: $GC_pairs	 AU pairs: $AU_pairs	 GU pairs: $GU_pairs	KL: $KL_GC_cont %GC\n\n\n";

		print $file_out "$ensemblediversity \n\n ";

		if ($failed_initial_design == 1){
			print $file_out " WARNING: Misfolding within the locked-sequence required altering the target structure. \n";
			print $file_out "          Sequence design failed for the inputted target structure.\n\n";
		}

		print $file_out "Kissing Loop List:\n";
		for (my $j=0; $j<($numKL); $j++){
			$i=$j+1;
			print $file_out "$i A,$Alist[$i],$i B,$Blist[$i],$Elist[$i]\n";
		}
		print $file_out "\n$non_cognates";

		print $file_out "\n\n\n";

		close $file_out;
}


# sub preview {
# 	qx(perl trace_analysis.pl pattern.txt seq.txt > $output_file_path\_Preview.txt);
# 	qx(echo '$name Iteration $attempts\n' >> $output_file_path\_Preview.txt);
# }

##############################################################################################
# Banner Graphics Subroutines
################################################

sub welcome { #prints welcome screen

&printer("\n____    ____ _  _ ____ _    _  _ ____  \n");
&printer("|__/ __ |___ |  | |  | |    |  | |__/ \n");
&printer("|  \\    |___  \\/  |__| |___  \\/  |  \\  v1.93\n\n\n");

&printer("       0                                                                                               \n");
&printer("      1111111                             1111  0                 1    0                            111\n");
&printer("     11 010111                        11111000 0              11111 1111                      1    011\n");
&printer("  0111  1111110110                 111111111   0110          00 0 111111110                   000111 01\n");
&printer("1 10111 1110  11111              000011111111 101001 1    00  1000 11111100111110            10 011100111\n");
&printer("11011001101 111111110           00001     111 100 01000  111100 00111     0111110       000011   1  11\n");
&printer("10101 1000  1111100110         0000 00    1 11111 11101  011100001 11     0 111111    010101      11   \n");
&printer("1000 1110   1110 011111      0000   00     1 1111111111  111 000101 1     00  1110111 0110000      11   \n");
&printer("010 0000    1     01111     0000   000     11 111111000 1110 0010  1     000    1011110010000 1     11  \n");
&printer("00000001   11      10111   000      00     1  111111000 1110000  11     00       0111111    000     1  \n");
&printer("010011     1        10111 00001      00     1   0111000 011111   1     000       001111      00        \n");
&printer(" 10001    11    1111 101111111111      0         01111  11111           0       01010111     00     00 \n");
&printer("1000      11     111 1 0110111 111             0  0111011011    0             01010110111    000     0 \n");
&printer("001       1     1111 0 11010 11 11     1     000   0111000     00             0000011100101  00      0 \n");
&printer("011                1  0001000 0        11     00  1101110 11  000     1     101111010111000  00      00\n");
&printer("11        0     0    0000 0010   11     11    000 1000111010 000     11     00111 111  11000          0\n");
&printer("        00      0   00100 10100   1     11      000000001111 00      1     11001  0111 011011   11    0\n");
&printer("       000     00 000000110001110 11     1       00001 111100       11    1111111 0111 0 01110  11     \n");
&printer("1      00      0 00000111011001110 1     11     00000   001010      1    11111 111 1110111111000  1     \n");
&printer("110   0        1001  100001111 011111     111110000      000010    11   1111  111  1100001110001  1   1\n");
&printer(" 11101      11100 1  11 1 1  1   111 1     1100000         010111111  11100   111       11  1111  11111\n");
&printer("  111110011110 0  1 1 1   0  00  111111111100001             11110011110001  111        11   0111111111\n");
&printer("      0100001   11           100 1111111100                    101010  101111            1  001 111011\n");
&printer("         11111   1               0  111111                         101011 0                 1101110    \n");
&printer("         0                       0     1                             010 0                   0 110     \n");
&printer("\n\n");
&printer("Written by Cody Geary. Copyright 2020.\n");

}

sub victory { #print fireworks

	&printer("                                 .''.                    \n");
	&printer("       .''.             *''*    :_\\/_:     .            \n");
	&printer("      :_\\/_:   .    .:.*_\\/_*   : /\\ :  .'.:.'.       \n");
	&printer("  .''.: /\\ : _\\(/_  ':'* /\\ *  : '..'.  -=:o:=-       \n");
	&printer(" :_\\/_:'.:::. /)\\*''*  .|.* '.\\'/.'_\\(/_'.':'.'      \n");
	&printer(" : /\\ : :::::  '*_\\/_* | |  -= o =- /)\\    '  *       \n");
	&printer("  '..'  ':::'   * /\\ * |'|  .'/.\\'.  '._____           \n");
	&printer("      *        __*..* |  |     :      |.   |' .---'|     \n");
	&printer("       _*   .-'   '-. |  |     .--'|  ||   | _|    |     \n");
	&printer("    .-'|  _.|  |    ||   '-__  |   |  |    ||      |     \n");
	&printer("    |' | |.    |    ||       | |   |  |    ||      |     \n");
	&printer(" ___|  '-'     '    ''       '-'   '-.'    '`      |____ \n");
	&printer("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n");
	&printer("\n\n");

	my $elapsed_time = time-$start_time;
	&printer("Total Design Time: $elapsed_time seconds\n\n");
	$start_time = time;  #reset the elapsed timer

	if ($failed_initial_design == 1){
		&printer(" WARNING: Misfolding within the locked-sequence required altering the target structure. \n");
		&printer("          Sequence design failed for the inputted target structure.\n\n");
	}
}

sub printer {
	if ( (time-$start_time)>$timeout_length ){die "\nComputation Halted - Timeout after $timeout_length seconds\n";}
	print $output_spool "$_[0]";
	print ("$_[0]");
}
