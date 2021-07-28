
# ++++++++++ SECTION 0 - SATISTICS AND MATH

sub mk_median {

my(@number) = @_;
my($count);
$count = @number;
my(@sort);
@sort = sort {$a <=> $b}(@number);
@number = @sort;
my($mid) = int($count/2);
my($median) = $number[$mid];
return ($count,$median);

}


sub mk_N50 {

my(@number) = @_;
my($count); 
my($sum);
my($i);
my(@sort);
my($half);
my($n50);

@sort = sort {$b <=> $a}(@number);
@number = @sort;
$count = @number;

for ($i=0;$i<@number;$i++) {
	$sum += $number[$i];
	}
$half = $sum/2;
$sum =0;

LOOPN50: for ($i=0;$i<@number;$i++) {
        $sum += $number[$i];
	if ($sum > $half) {
		$n50 = $number[$i];
		last LOOPN50;
		} 
        }
return ($count,$n50);

}



# calculate standard deviation
sub mk_STD {

my(@number) = @_;
my($count);
my($sum);
my($i);
my($avg);
my($variance);
my($std);
my($sq);
my(@sort);

# calc avg
$sum = 0;

for($k=0;$k<@number;$k++) {
        $sum += $number[$k];
        #print "$j\t$hit[$k][0]\t-$hit[$k][1]\t$sizes[$k]\n";
        }
$avg = (int($sum/@number*10))/10;

# calc standard deviation
$sum = 0;
for($k=0;$k<@number;$k++) {
        $sq = ($number[$k]-$avg)**2;
        $sum+=$sq;
        }
$variance = $sum/@number;
$std = (int((sqrt($variance))*100))/100;

return ($avg,$std);





}
















# ++++++++++ SECTION I - SEQUENCE FORMAT CONVERSIONS



# &&&&&&&&&&&&& NEW SUBROUTINE - make fasta format
sub mk_fasta {

my($seq) = @_;
my($i);
my($back);
my($sub);

$seq =~ s/ //g;

for($i=0; $i <= length($seq);$i+=50) {
	$sub = substr($seq,$i,50);
	$back .= $sub;
	$back .= "\n";
	}

return ($back);

} # end mk_fasta



# &&&&&&&&&&&&& NEW SUBROUTINE - make fasta format with 60 bp lines

sub mk_fasta60 {

my($seq) = @_;
my($i);
my($back);
my($sub);

for($i=0; $i <= length($seq);$i+=60) {
        $sub = substr($seq,$i,60);
        $back .= $sub;
        $back .= "\n";
        }

return ($back);

} # end mk_fasta


# compare two sequences and make nice alignment output
sub nice_alignment {

my(@seq) = @_;
my($seq1) = $seq[0];
my($seq2) = $seq[1];
my($k);
my($out_align) = '';
my($sub1);
my($sub2);
my($sub3);
my($line) = 100;

my($align) = &compare_2_seqs($seq1,$seq2);

for ($k =0; $k <= length($align); $k+=$line) {
	$sub1 = substr($seq1,$k,$line);
	$sub2 = substr($align,$k,$line);
	$sub3 = substr($seq2,$k,$line);
	$out_align .=  "SEQ1 $sub1\nALN  $sub2\nSEQ2 $sub3\n";
	}

return($out_align);

}


# translate inseq

sub translate {

my($seq) = @_;
my($i);
my($protein) = '';
my($codon);

for ($i = 0; $i < (length($seq) -2); $i+=3) {
	$codon = substr($seq, $i, 3);
	$protein .= codons($codon);
	}

return($protein);

}






# +++++++++++++++++++++  pure sequnece +++++++++++++++++++++++++++
# &&&&&&&&&&&&& NEW SUBROUTINE - PURE_SEQUENCE
# extracts ony the seqeunce part from any sequence file
# in EMBL, GenBank or FASTA format and returns it as 1 variable
# for the seqeunce and optional others fro DE line etc.
# probably meets problems with seqeunces >> 1 000 000 bp 

sub pure_sequence {

my ($infile) = @_;

# suck file in one variable

open (INPURE, "<$infile");
my(@temp) = <INPURE>;
my($temp) = join('',@temp);
close (INPURE);

# distinguish fasta, EMBL and genbank

if ($temp =~ /^LOCUS/) {
        (@pure) = &pure_gb_sequence($infile);
        } elsif ($temp =~ /^>/) {
        (@pure) = &pure_fasta_sequence($infile);
        } else {
        (@pure) = &pure_embl_sequence($infile);
        }
return(@pure);

}

# &&&&&&&&&&&&& SUB-SUBROUTINE - PURE_GENBANK_SEQUENCE
# extracts seqeunce from GenBank file

sub pure_gb_sequence {

my ($infile) = @_;

open (INGB,"<$infile");
my (@temp) = (<INGB>);
my ($temp) = join('',@temp);
close (INGB);

my ($de) = ($temp =~ /^DEFINITION\s+(\S+.+)\n^ACCESSION/ms);
$de =~ s/\n//;
my ($seq) = ($temp =~ /ORIGIN(.+)\/\//s);

$seq =~ s/\d//g;
$seq =~ s/\s+//g;
$seq =~ tr/[a-z]/[A-Z]/;

# return sequence and definition line

return($seq, $de);

}

# &&&&&&&&&&&&& SUB-SUBROUTINE - PURE_FASTA_SEQUENCE

sub pure_fasta_sequence {

my ($infile) = @_;

open (INFASTA,"<$infile");
my (@temp) = (<INFASTA>);
my ($temp) = join('',@temp);
close (INFASTA);

my ($de) = ($temp =~ /^>([^\n]+)/);
my ($seq) = ($temp =~ /^>[^\n]+\n(.+)/s);

$seq =~ s/\s+//g;
$seq =~ tr/[a-z]/[A-Z]/;

# return sequence and definition line

return($seq, $de);

}

# &&&&&&&&&&&&& SUB-SUBROUTINE - PURE_FASTA_SEQUENCE with quality data

sub pure_fasta_sequence_qual {

my ($infile) = @_;

open (INFASTA,"<$infile");
my (@temp) = (<INFASTA>);
my ($temp) = join('',@temp);
close (INFASTA);

my ($de) = ($temp =~ /^>([^\n]+)/);
my ($seq) = ($temp =~ /^>[^\n]+\n(.+)/s);

$seq =~ s/\s+//g;

# return sequence and definition line

return($seq, $de);

}


# &&&&&&&&&&&&& SUB-SUBROUTINE - PURE EMBL SEQUENCE

sub pure_embl_sequence {

my ($infile) = @_;

open (INEMBL,"<$infile");
my (@temp) = (<INEMBL>);
my ($temp) = join('',@temp);
close (INEMBL);

my ($de) = ($temp =~ /^DE\s+(\S.+)$/m);
my ($seq) = ($temp =~ /\.\.(.+)/s);


$seq =~ s/\d//g;
$seq =~ s/\s+//g;
$seq =~ s/\n//g;

$seq =~ tr/[a-z]/[A-Z]/;
my ($len) = length($seq);

return($seq,$de);

}
# +++++++++++++ end pure_sequence ++++++++++++++++++++++++++++++++




# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++ SECTION II - PARSERS






# &&&&&&&&&&&&&& NEW SUBROUTINE - make size list of chromsomes in /data/dir_chromsomes

sub get_chr_size {

my(@data) = @_;
my($chr_query) = $data[0];
my(@names);
my(@size);
my(%chr);
my($chr_count);
my($name);
my($file) = "/home/wicker/data/dir_chromosomes/size_list";

#`ls -l /home/wicker/data/dir_chromosomes/  > chr_size_list`;

print "file = $file\t$chr_query\n";

open(CHR_SIZE,"$file");
while (<CHR_SIZE>) {
	($len,$name) = /^\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+.+\s+.+($chr_query\S+)/;
	if ($len) {
		$len = int($len*59/60);
		

		$entry = "$name-$len";
		push(@names,$entry);
		}
	}

close(CHR_SIZE);	
return(@names);

}




# &&&&&&&&&&&&&& NEW SUBROUTINE - PARSE FOR REPEAT KEYS

sub parse_repeat {

my($string) = @_;
my($test) = 0;

my(@reps) = ("transpos", "reverse", "integrase", "gag-pol", "polyprotein", 
"gypsy", "copia", "mutator", "repetitive", "cr1", "cr-1", "mariner", 
"TNP2", "CACTA", "retroelement", "retrovir", " alu");

foreach $rep (@reps) {
	if ($string =~ /$rep/i) {
		$test =1;
		}
	}

if ($test ==1) {
	$string =~ s/Similar to/sim/i;
        $string =~ s/putative//i;
        $string =~ s/hypothetical protein/HyP/i;
        $string =~ s/hypothetical//i;
        $string =~ s/unknown protein/UP/i;
        $string =~ s/with similarity to/sim/i;
        $string =~ s/probable//;
	$string =~ s/reverse transcriptase/RT/i;
        $string =~ s/integrase/INT/i;
	}

if ($test ==0) {
	my(@hyps) = ("hypothetical protein", "unknown", "putative protein", "unnamed");
	foreach $rep (@hyps) {
	        if ($string =~ /$rep/i) {
			$test =2;
	                $string = $rep;
	                }
	        }
	}
if ($test ==2) {
	$string =~ s/unnamed/UnP/i;
        $string =~ s/putative protein/PuP/i;
        $string =~ s/hypothetical protein/HyP/i;
        $string =~ s/unknown/UkP/i;
        $string =~ s/Similar to//i;
        $string =~ s/probable//;
        $string =~ s/with similarity to//i;
        }

return($test,$string);

}


# &&&&&&&&&&&&& NEW SUBROUTINE - PARSE SMITH WATERMAN

sub parse_smith {

my($file) = @_;

my $rethit = '';
my ($seq1,$seq2,$ide,$ev,$beg1,$beg2,$query,$subj) = '';
my (@seq) = ();
my (@file) =();
my ($smith) = '';
my ($query) = '';
my ($subj) = '';

open (IN,"<$file");

@file = <IN>;
$smith = join('',@file);

# get seqeunce names and create identifier

($seq1) = ($smith =~ /#[^\n]+sequences: 2\n#\s+1:\s+([^\n]+)/s);
($seq2) = ($smith =~ /ned_sequences: 2\n[^\n]+\n#\s+2:\s+([^\n]+)/s);
$ID = "${seq1} x ${seq2}";

# it only works if both seqeunces have a name

if (($seq1 eq '')||($seq2 eq '')) {
        $ID = "WARNING - sequence without name!!";
        }

# extract info

($ide) = ($smith =~ /#\s+Similarity:\s+\S+\s+[^\d]+(\d+.+)%/);
($ev) = ($smith =~ /#\s+Score:\s+([\d\.]+)/);
($beg1) = ($smith =~ /Score.+#[\s\n]+#[=\s\n]+\S+\s+(\d+)/sm);
($beg2) = ($smith =~ /#[=\s\n]+\S+\s+\d+[^\|]+[\|\.\:\s\n]+\S+\s+(\d+)/s);
($end1) = ($smith =~ /\s+(\d+)[\s\n]+[^\n]+[^\|]+[\s\n]+#-+\n#-+/s);
($end2) = ($smith =~ /(\d+)[\s\n\|]+#-+\n#-+/s);
(@seq) = ($smith =~ /^\S+\s+\d+\s+(\S+)\s+\d+[\s\n]+/smg);

$line =@seq;

# create the two sequences from @seq
for ($i =0; $i <$line; $i+=2){
        $query .= $seq[$i];
        }

for ($i =1; $i <$line; $i+=2){
        $subj .= $seq[$i];
        }

# create return variable

# hitindex can be provided by the main program (optional);

if ($hitindex eq '') { $hitindex = 1; }

$rethit = "ID = ${ID}\ndata=${hitindex};${ev};${ide};${beg1};${end1};";
$rethit .= "${beg2};${end2};${query};${subj}\n";

return ($rethit);

}

                                                                                

# ++++++++++++++++++++++ parse_blast
# &&&&&&&&&&&& NEW SUBROUTINE - PARSE_BLAST


sub parse_blast {

my ($file) = @_;
my (@hit)=();
my (@onehit) =();
my ($count_hit) =0;
my ($check) = 0;
my ($blastx) = 0;
my ($score);
my ($len);
my ($eval);
my ($ID);
my ($rem);
my ($this_hit);
my ($id);
my ($beg_check)=0;
my ($qu);
my ($be);
my ($en);
my ($su);
my ($query);
my ($subj);
my ($beg1);
my ($end1);
my ($beg2);
my ($end2);
my ($re);
my ($de_second) = 0;

if ($file =~ /blastx/) {$blastx =1;}

# walk through blastfile
open (IN_BLAST, "<$file");
                                                                                
while (<IN_BLAST>){

        # extract IDs from DE line, initiate multiline description
        ($ID) = /^>\s{0,1}(\S+)/;

        if ($ID) {
		# replace pipes
		$ID =~ s/\|/_/g;

		# grab everything after the ID
		($rem) = /$ID\s+(\S.+)/;

		$rem =~ s/;//g;

		if ($rem eq ''){
			$rem = $ID;
			}

                if ($query ne '') {
                        $this_hit .= "$beg1;$end1;$beg2;$end2;$query;$subj\n";
			$query = $subj = '';
                        }

                push(@hit,$this_hit);
		$de_second=1;
                $this_hit = "ID=$ID;$rem";
                }


	# consider multiple lines of descriptions
        if ($de_second ==1) {
                ($string) = /^([^>].+)$/;
                unless ($string =~ /Length/) {
                        $string=~ s/;//g;
                        $this_hit .= " $string";
                        }
                }


        ($len) = /Length\s*=\s*(\d+)/;
        if ($len) {
                $this_hit =~ s/\s+/ /g;
                $this_hit =~ s/\s$//g;
                $this_hit .= ";";
                $de_second =0;
                $this_hit .="$len";
        	}        



	# extract score, E-value
	($score,$eval) = /^\s+Score\s+=\s+(\S+)\s+.+Expect.+=\s+(\S+)/;
	if ($score) {
		# add the preivious alignment
		if ($query ne '') {
			$this_hit .= "$beg1;$end1;$beg2;$end2;$query;$subj";
			$query = $subj = '';
			}
		$eval =~ s/,//;
		$this_hit .= "\ndata=$score;$eval;";
		}

	($id) = /^\s+Identities\s+=\s+\d+\/\d+\s+\((\d+)%\)/;
	if ($id) {
		$this_hit .=  "$id;";
		$beg_check =0;
		}

	# extract query alignemnt data
	($be,$qu,$en) = /Query[:]{0,1}\D+(\d+)\s+(\S+)\s+(\d+)/;
	if ($be) {
		# check if beginnin of alignment
		if ($beg_check ==0) { 
			$beg1 = $be; 
			}
		$query .= "$qu";
		$end1 = $en;
		}

	# extract subject alignment data
	($be,$su,$en) = /Sbjct[:]{0,1}\D+(\d+)\s+(\S+)\s+(\d+)/;
        if ($be) {
                # check if beginnin of alignment
                if ($beg_check ==0) {
                        $beg2 = $be;
			$beg_check=1;
                        }
                $subj .= "$su";
                $end2 = $en;
                }

	}


# add the last alignment and the last hit
$this_hit .= "$beg1;$end1;$beg2;$end2;$query;$subj\n";
push(@hit,$this_hit);

# get rid of the first (empty) element
shift(@hit);

return (@hit);

} # ---------- end parse_blast 














# &&&&&&&&&&& NEW SUBROUTINE - PARSE_ANNOTATION
# pulls out data from annotated files and reproduces them 
# in simple table form

sub parse_annotation {

my ($infile) = @_;

# distinguish EMBL and genbank and
# call subroutines to make tmp file with only annotation

open (IN_ANN, "<$infile");
@temp = <IN_ANN>;
$temp = join('',@temp);
close (IN_ANN);

if ($temp =~ /^LOCUS/) {
        &mk_genbank_array;
        } else {
        &mk_EMBL_array;
        }


# &&&&&&&&& catch annotation part from EMBL files

sub mk_EMBL_array {

($annot) = ($temp =~ /FT\s+source(.+\.\.)/s);
$annot =~ s/^XX.+\.\.//ms;

open (OUT_ARRAY, ">tempfile");
print OUT_ARRAY "$annot";
close (OUT_ARRAY);

# +++++++ parse annotation, put features in one array

open (IN, "<tempfile");

$count = 0;

while (<IN_EMBL>) {
        ($go) = /^FT\s{2,4}\S+/;
        if ($go) { $count ++; }
        $feature[$count] .= $_;
        }

close (IN_EMBL);


# get rid of the first element "source etc."
shift(@feature);

}
# %%%%%%%% end of mk_EMBL_array

# &&&&&&&&& subroutine specific for genbank files

sub mk_genbank_array {

my ($infile) = @_;

# make tmp file with only annotation

open (IN_ARRAY, "<$infile");
@temp = <IN_ARRAY>;
$temp = join('',@temp);
close (IN_ARRAY);

($annot) = ($temp =~ /FEATURES.+source(.+)BASE\sCOUNT/s);
$annot =~ s/^XX.+\.\.//ms;

open (OUT_TEMP, ">tempfile");
print OUT_TEMP "$annot";
close (OUT_TEMP);

# +++++++ parse annotation, put features in one array

open (IN_GENB, "<tempfile");

$count = 0;

while (<IN_GENB>) {
        ($go) = /^\s{4,6}\S+/;
        if ($go) {
                $count ++;
                }

        $feature[$count] .= "FT";
        $feature[$count] .= $_;
        }

close (IN_GENB);
shift(@feature);

}

# %%%%%%%% end of mk_genbank_array

# +++++++++ PART II - process data and write outfile
# process features from array one by one

$out = "${infile}_visual";

open (OUT_TABLE, ">$out");

foreach (@feature) {
        ($feat, $rem) = &process_feature($_);
        }

# &&&&&&&& subroutine - extract feature data

sub process_feature {

my ($feat) = @_;
my ($rem_add) = '';

if ($feat =~ /mRNA/) { return;}
if ($feat =~ /intron/) { return;}

if ($feat =~ /complement/) {
        $orient = "-";} else { $orient = "+"; }

# proces info from /rpt_type= which is not between ".."

($rem_add) = ($feat =~ /rpt_type=(.+)$/);
$rem_add =~ tr/[A-Z]/[a-z]/;
if ($rem_add ne '') {
        $rem_add = "${rem_add} rep;";
        }

# compress annotation notes (kick out blah-blah)

$feat =~ s/\/protein_id="[^"]+"//;
$feat =~ s/SPTREMBL://;
$feat =~ s/SWISS-PROT://;
$feat =~ s/\/organism="[^"]+"//;
$feat =~ s/^FT\s+//gm;
$feat =~ s/similar to/sim/g;
$feat =~ s/similar/sim/g;
$feat =~ s/putative//g;
$feat =~ s/\/translation="[^"]+"//;
# collect all remaining notes between ".." in $rem
# and add $rem_add at the beginning

(@rem) = ($feat =~ /"([^"]+)"/gm);
unshift (@rem, $rem_add);
$rem = join(' ',@rem);
$rem =~ s/\n/ /sg;
$rem =~ s/\s+/ /g;

# collect all position info in @pair

($pos) = ($feat =~/([^\/]+)\//sm);
$pos =~ s/\n//g;
(@pair) = ($pos =~ /(\d+\.\s{0,1}\.\d+)/g);

# write table to outfile

$level =1;
$frag_count = 1;

foreach (@pair) {
        if (length($rem) > 60) {
                ($rem_temp) = substr($rem,0,60);
                } else {
                $rem_temp = $rem;
                }
        $rem_temp = "${rem_temp}\.${frag_count}";
        $_ =~ s/\./ /g;
        print OUT_TABLE "$_  $level  $orient $rem_temp\n";
        $frag_count++;
        }


return ($pos, $rem);

}

close (OUT_TABLE);

}
# end of PARSE_ANNOTATION



# +++++++++ SECTION III - SEQUENCE PROCESSING SUBROUTINES

# &&&&&&&&&&&&&&& NEW SUBROUTINE - fasta_from_flat
# extracts one sequence from a large flat file
# returns pure sequnece, not fasta

sub fasta_from_flat {

my ($file,$name) =@_;
my ($flag) =0;
my ($seq) = '';

print "\n=====================================================\n";
print "==========  WICKERsoft database search ==============\n";
print "===============  \"DNA_from_flat\" ==================\n\n";
print "search file $file\nfor $name\n";

open (IN_TEMP,"<$file");

while (<IN_TEMP>){

	if (/^>/) { 
		$flag =0; 
		print "$_";
		}

	if (/^>.*$name/) {$flag =1;}

	if (($flag ==1) && ($_ !~ /^>/)) {
		$seq .= $_;
		}

	if (($flag ==0) && ($seq ne '')){
		last;
		}

	}

$seq =~ s/\n//g;

my ($len) = length($seq);

print "extracted $name, size = $len\n";
print "=====================================================\n\n";


return ($seq);

}



# &&&&&&&&&&&&&&& NEW SUBROUTINE - collapse array
# looks for redundancy in an array of elements and gets rid of everything
# that occurs in duplicate

sub collapse_array {

my (@sort1);
my (@sort2);
my(@list) = @_;
my(@collect) =();
my($i);


my (@sort1) = sort (@list);
my (@sort2) = sort {$a<=>$b} (@sort1);

for ($i = 0; $i <@sort2;$i++) {
	if ($sort2[$i] ne $sort2[$i+1]) {
		push(@collect,$sort2[$i]);
		}
	}

return (@collect);

} 



# &&&&&&&&&&&&&&& NEW SUBROUTINE - collapse array
# looks for redundancy in an array of elements and gets rid of everything
# that occurs in duplicate

sub collapse_array_ascii {

my (@sort);
my(@list) = @_;
my(@collect) =();
my($i);

#my (@sort) = sort {$a<=>$b} (@list);

my (@sort) = sort (@list);

for ($i = 0; $i <@sort;$i++) {
        if ($sort[$i] ne $sort[$i+1]) {
                push(@collect,$sort[$i]);
                }
        }

return (@collect);

}





# &&&&&&&&&&&&&&& NEW SUBROUTINE - COMPARE 2 SEQUENCES
# simple seqeunce comparison (not alignement!) of 2 sequences
# which ARE already aligned
# returns a variable it " " or "|" for un/matching bases

sub compare_2_seqs {

my($seq1,$seq2) = @_;
$seq1 =~ tr/a-z/A-Z/;
$seq2 =~ tr/a-z/A-Z/;
my (@seq1) = split('',$seq1);
my (@seq2) = split('',$seq2);
my ($comp) = '';
my ($count) =0;
my ($match) = 0;
my ($i) = 0;
my ($align) = 0;
my ($rel) = 0;
my ($len) = length($seq1);

for($i=0;$i<@seq1;$i++) {

	if ($seq1[$i] eq $seq2[$i]) {
		$comp .= "|";
		}

	if ($seq1[$i] ne $seq2[$i]) {
                $comp .= " ";
                }

	if (($seq1[$i] eq $seq2[$i]) && ($seq1[$i] ne '-') && ($seq2[$i] ne '-') && ($seq1[$i] ne 'n') && ($seq2[$i] ne 'n')) {
		$match++;
		}

	if (($seq1[$i] ne '-') && ($seq2[$i] ne '-')) {
		$align++;
		}

        }

#print "$len\t$match/$align\n";

if ($align > 0) {
	$rel = (int($match/$align*10000))/100;
	}

return ($comp,$rel);

}

# &&&&&&&&&&&&&&& NEW SUBROUTINE - INVERT SEQUENCE
# produces the reverse complement of inpu seqeunce

sub invert_sequence {

my ($inseq) = @_;

#$inseq =~ tr/acgt/ACGT/;

my(@inseq) = split('',$inseq);
my(@rev) = reverse(@inseq);
my($rev) = join('',@rev);

$rev =~ tr/ACGT/TGCA/;
$rev =~ tr/acgt/tgca/;

return($rev);

}

# &&&&&&&&&&&&& NEW SUBROUTINE - GC content

sub gc_content {

my (@dna) = ();
my ($dna) = @_;
my ($len) = 0;
my ($a) = 0;
my ($c) = 0;
my ($g) = 0;
my ($t) = 0;
my ($gc) =0;
my ($GC) =0;
my ($n) = 0;
my (@n) =0;

(@n) = ($dna =~ /N/g);
$n = @n;

@dna = split('',$dna);
$len = (length($dna))-$n;

foreach $base (@dna) {
        if ($base eq "A") { $a++; }
        if ($base eq "C") { $c++; }
        if ($base eq "G") { $g++; }
        if ($base eq "T") { $t++; }
	}

$GC = int(($g+$c)/$len*10000);
$gc = $GC/100;

return ($gc);

}

# &&&&&&&&&&&&& NEW SUBROUTINE - base_composition

sub base_composition {

my (@dna) = ();
my ($dna) = @_;
my ($len) = 0;
my ($a) = 0;
my ($c) = 0;
my ($g) = 0;
my ($t) = 0;
my ($A) = 0;
my ($C) = 0;
my ($G) = 0;
my ($T) = 0;
my ($gc) =0;
my ($GC) =0;

$len = length($dna);
@dna = split('',$dna);

foreach $base (@dna) {
        if ($base eq "A") { $A++; }
        if ($base eq "C") { $C++; }
        if ($base eq "G") { $G++; }
        if ($base eq "T") { $T++; }
        }

if ($len >0){
	$a = int(($A/$len)*100);
	$c = int(($C/$len)*100);
	$g = int(($G/$len)*100);
	$t = int(($T/$len)*100);
	} else {
	print "$infile is zero length*****\n";
	}

return ($a,$c,$g,$t);

}

# &&&&&&&&&&&&& NEW SUBROUTINE - species code

sub species_code {

my ($species) = @_;

my (%species_code) = (

'Ggal' => 'Gallus gallus',
'Grai' => 'Gossypium raimondii',
'Garb' => 'Gossypium arboreum',
'Gste' => 'Gossypium stertianum',
'Gster' => 'Gossypium stertianum', 
'Gher' => 'Gossypium herbaceum',
'Gano' => 'Gossypium anomalum',
'Gtri' => 'Gossypium triphyllum',
'Gsto' => 'Gossypium stocksii',
'Gsom' => 'Gossypium somalense',
'Glon' => 'Gossypium longicalix',
'Gbic' => 'Gossypium bickii',
'Osat' => 'Oryza sativa',

);


if (exists $species_code{$species}) {
        return $species_code{$species};
        } else {
        return "undefined";
        }

}


# &&&&&&&&&&&&& NEW SUBROUTINE - TREP classifictaions

sub TREP_class {

my ($class) = @_;
my ($test);
my ($guy);

my (@guys) = qw(TRIM MITE non-LTR gypsy copia athila LINE SINE stowaway stowaway tourist cacta mutator lite);

foreach $guy (@guys) {
	if ($class =~ /$guy/i) {
		$test = $guy;
		}
	}

$test = lc $test;

my (%classes) = (

'gypsy' => 'retrotransposon, LTR, gypsy',
'copia' => 'retrotransposon, LTR, copia',
'athila' => 'retrotransposon, LTR, athila',
'line' => 'retrotransposon, non-LTR, LINE',
'sine' => 'retrotransposon, non-LTR, SINE',
'stowaway' => 'foldback element, MITE, stowaway',
'tourist' => 'foldback element, MITE, tourist',
'tourist' => 'foldback element, MITE, tourist',
'cacta' => 'transposon, CACTA',
'mutator' => 'transposon, mutator',
'lite' => 'foldback element, LITE',
'mite' => 'foldback element, MITE',
'trim' => 'retrotransposon, TRIM',

);

if (exists $classes{$test}) {
        return $classes{$test};
        } else {
	return "unclassified";
	}

}




# TREP code, convertd 3-letter code into human readable classification
sub TREP_code {

my ($code) = @_;
my($class) = '';

my (%classes) = (

'RLG' => 'retrotransposon, LTR, Gypsy',
'RLC' => 'retrotransposon, LTR, Copia',
'RLH' => 'retrotransposon, LTR, Halcyon',
'RLF' => 'retrotransposon, LTR, Echo',
'RLX' => 'retrotransposon, LTR, unknown',

'RPA' => 'retrotransposon, LTR, Pan',

'RIC' => 'retrotransposon, LTR, Chronos',
'RIR' => 'retrotransposon, non-LTR (LINE), R2',
'RIJ' => 'retrotransposon, non-LTR (LINE), Jokey',
'RIL' => 'retrotransposon, non-LTR (LINE), L1',
'RII' => 'retrotransposon, non-LTR (LINE), I',
'RIX' => 'retrotransposon, non-LTR (LINE), unknown',

'RSX' => 'retrotransposon, non-LTR (SINE), unknown',

'DTC' => 'DNA transposon, TIR, CACTA',
'DTH' => 'DNA transposon, TIR, Harbinger',
'DTT' => 'DNA transposon, TIR, Mariner',
'DTM' => 'DNA transposon, TIR, Mutator',
'DTA' => 'DNA transposon, TIR, hAT',
'DTX' => 'DNA transposon, TIR, unknown',

'DHH' => 'DNA transposon, Helitron, Helitron',

'DXX' => 'DNA transposon, unknown, unknown',

'XXX' => 'unknown, unknown, unknown',

);

if (exists $classes{$code}) {
        return $classes{$code};
        } else {
        return "NO CODE";
        }

}








# &&&&&&&&&&&&& NEW SUBROUTINE - amino acid names

sub aa_similarity {


my ($aa) = @_;

# these are aa exchanges of Blumeria

my (%aa_names) = (

'A' => 'TVS',
'C' => 'SYR',
'D' => 'EN',
'E' => 'DKG',
'F' => 'LYS',
'G' => 'SEDA',
'H' => 'QYN',
'I' => 'VTL',
'K' => 'RNE',
'L' => 'IFSP',
'M' => 'ILT',
'N' => 'SDK',
'P' => 'SLT',
'Q' => 'HEKR',
'R' => 'KQGS',
'S' => 'PTNA',
'T' => 'AIS',
'V' => 'IA',
'W' => 'RLCS',
'Y' => 'HFC',



);

if (exists $aa_names{$aa}) {
        return $aa_names{$aa};
        } else {
        return "x";
        }


}


# &&&&&&&&&&&&& NEW SUBROUTINE - amino acid names

sub aa_name {

my ($aa) = @_;

my (%aa_names) = (

'A' => 'Ala',
'C' => 'Cys',
'D' => 'Asp',
'E' => 'Glu',
'F' => 'Phe',
'G' => 'Gly',
'H' => 'His',
'I' => 'Ile',
'K' => 'Lys',
'L' => 'Leu',
'M' => 'Met',
'N' => 'Asn',
'P' => 'Pro',
'Q' => 'Gln',
'R' => 'Arg',
'S' => 'Ser',
'T' => 'Thr',
'V' => 'Val',
'W' => 'Trp',
'Y' => 'Tyr',

);



if (exists $aa_names{$aa}) {
        return $aa_names{$aa};
        } else {
        return "x";
        }



}



# &&&&&&&&&&&&& NEW SUBROUTINE - genetic code

sub codons {

my ($codon) = @_;

$codon = uc $codon;

my (%genetic_code) = (

'TCA' => 'S',
'TCC' => 'S',
'TCG' => 'S',
'TCT' => 'S',
'TTC' => 'F',
'TTT' => 'F',
'TTA' => 'L',
'TTG' => 'L',
'TAC' => 'Y',
'TAT' => 'Y',
'TAA' => '*',
'TAG' => '*',
'TGC' => 'C',
'TGT' => 'C',
'TGA' => '*',
'TGG' => 'W',
'CTA' => 'L',
'CTC' => 'L',
'CTG' => 'L',
'CTT' => 'L',
'CCA' => 'P',
'CCC' => 'P',
'CCG' => 'P',
'CCT' => 'P',
'CAC' => 'H',
'CAT' => 'H',
'CAA' => 'Q',
'CAG' => 'Q',
'CGA' => 'R',
'CGC' => 'R',
'CGG' => 'R',
'CGT' => 'R',
'ATA' => 'I',
'ATC' => 'I',
'ATT' => 'I',
'ATG' => 'M',
'ACA' => 'T',
'ACC' => 'T',
'ACG' => 'T',
'ACT' => 'T',
'AAC' => 'N',
'AAT' => 'N',
'AAA' => 'K',
'AAG' => 'K',
'AGC' => 'S',
'AGT' => 'S',
'AGA' => 'R',
'AGG' => 'R',
'GTA' => 'V',
'GTC' => 'V',
'GTG' => 'V',
'GTT' => 'V',
'GCA' => 'A',
'GCC' => 'A',
'GCG' => 'A',
'GCT' => 'A',
'GAC' => 'D',
'GAT' => 'D',
'GAA' => 'E',
'GAG' => 'E',
'GGA' => 'G',
'GGC' => 'G',
'GGG' => 'G',
'GGT' => 'G'

);

if (exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
        } else {
        #print "Bad codon: $codon!!\n";
        return "x";
        }

}




# SUBROUTINE - 3 letter code

sub TE_code {

my ($code) = @_;

my (%TE_code) = (

'RLC' => 'Copia',
'RLG' => 'Gypsy',
'RLX' => 'RLX',
'RIX' => 'LINE',
'RSX' => 'SINE',
'DTT' => 'Mariner',
'DTA' => 'hAT',
'DTM' => 'Mutator',
'DTP' => 'P',
'DTH' => 'Harbinger',
'DTC' => 'CACTA',
'DHH' => 'Helitron',
'DTX' => 'DTX',
'XXX' => 'unclassified',
);

return $TE_code{$code};

}

sub seq2date { 

 (my $s, my $t, my $v) = @_;
 
 #the substitution rate used for div time calculations
my $subst_rate = 1.3*(10**-8); # for grasses (Ma and Bennetzten, 2004)
# my $subst_rate = 2.22*(10**-10); # for 25S rDNA (Wicker)
# my $subst_rate = 6.8*(10**-11); # for chloroplast C8 (Wicker)
#my $subst_rate = 1.5*(10**-10); # for triticeae chloroplast (Wicker)
#my $subst_rate = 6.5*(10**-9); # for grasses (Gaut)
#my $subst_rate = 6.5*(10**-10); # for fungi 28S (Takamatsu)
#my $subst_rate = 6.8*(10**-10); # for pines (BMC Evol Biol. 2012;12:8)

#my $subst_rate = 1.5*(10**-8); # for Atal (Koch et al. 2000)
#my $subst_rate = 5.75*(10**-9); # for Solanaceae   (Mol. Biol. Evol. 2017, 34:1363–1377)
#my $subst_rate = 6.5*(10**-9); # for Brassicaceae (Mol. Biol. 2017, Evol. 2017, 34:1363–1377)
#my $subst_rate = 6.3*(10**-9); # for Grasses      (Mol. Biol. Evol. 2017, 34:1363–1377)
#my $subst_rate = 6.3*(10**-9); # for Pinaceae     (Mol. Biol. Evol. 2017, 34:1363–1377)

$ratio = 1-2*($v/$s);

if (($ratio > 0) && ((1-2*($t/$s)-($v/$s)) > 0)) {
 
 #the K2P formula
 my $k2p = -0.5*log((1-2*($t/$s)-($v/$s)) * sqrt(1-2*($v/$s)));
    $k2p =  sprintf("%.6f", $k2p);
 #print "$k2p\n";

 #the STD of K2p
 #split the whole moster-formula in four parts for easier editing
 my $t1 = 1/(1-2*($t/$s)-($v/$s))**2*($t/$s);
 my $t2 = (0.5 * ((1/((1-2*($t/$s)-($v/$s))))+(1/(1-2*($v/$s)))))**2*($v/$s);
 my $t3 = 1/(1-2*($t/$s)-($v/$s))*($t/$s);
 my $t4 = (0.5*((1/((1-2*($t/$s)-($v/$s))))+(1/(1-2*($v/$s)))))*($v/$s);

 my $k2p_std = sqrt((1/$s)*($t1 + $t2 -($t3 + $t4)**2));
    $k2p_std = sprintf("%.6f", $k2p_std);
 
 #the actual divergence time
 my $div_time     = sprintf("%.0f", $k2p / (2*($subst_rate)));
 my $div_time_std = sprintf("%.0f", $k2p_std / (2*($subst_rate)));
 my $out  = "-------\nDiv.time (MYA): $div_time    +/- $div_time_std\n--------\n";
    $out .= "K2P:\t $k2p\nK2P STD: $k2p_std\n";

my @out = ($div_time, $div_time_std);


    return(@out);
} else {
my $div_time = "NO DATE: 1-2*(TV/SITES) < 0";
my $div_time_std = "NO DATE: 1-2*(TV/SITES) < 0";
my @out = ($div_time, $div_time_std);
return(@out);
}

}



# compare 2 seqgs with WATER and write reformatted output +++++++++++++++++++++++++++++++++++++

sub bestfit {

my(@seq) = (@_);

my ($seq1) = $seq[2];
my ($seq2) = $seq[3];
my ($name1) = $seq[0];
my ($name2) = $seq[1];


open(WATER,">$name1");
($fasta) = &mk_fasta($seq1);
print WATER ">$name1\n$fasta\n";
close(WATER);

open(WATER,">$name2");
($fasta) = &mk_fasta($seq2);
print WATER ">$name2\n$fasta\n";
close(WATER);

$outfile = "$name1\__x__$name2\.pair";
$outfile2 = "$infile1\__x__$infile2\.pair_tmp";

`water $name1 $name2 -gapopen=50 -gapextend=0.1  -outfile=$outfile`;

`rm $name1 $name2`;

print "\n\nwrote $outfile\n\n";


# reformat output file +++++++++++++++++++++++++++++++++

open (OUT,">$outfile2");
open (IN2, "$outfile");

$which = 0;
$space = "                     ";


while (<IN2>){

	($nam1) = /# 1:\s+(\S+)/;
	if ($nam1) {
		$name1 = $nam1;
		$name1 .= "                         ";
		$name1 = substr($name1,0,30);
		}

	($nam2) = /# 2:\s+(\S+)/;
        if ($nam2) {
                $name2 = $nam2;
		$name2 .= "                         ";
                $name2 = substr($name2,0,30);
                }

	if (/#/) {
		print OUT;
		}

        ($go) = /^$/;
        if ($go) {
		$which = 0;
		print OUT;
		}

        ($name,$be,$seq,$en) = /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/;

        ($al) = /^$space([\s\|\.:\s]+)$/;
        if ($al) {
		chomp($al);
		print OUT "                                         $al\n";
                }

        if ($seq && $which ==1) {
                $be += $beg1;
                $en += $beg1;

		$pos = "$be          ";
		$be = substr($pos,0,9);

		print OUT "$name1 $be $seq $en\n";

                }

        if ($seq && ($which ==3)) {
                $be += $beg2;
                $en += $beg2;
                $which = 0;

                $pos = "$be          ";
                $be = substr($pos,0,9);

                print OUT "$name2 $be $seq $en\n";

                }

        $which++;

        }

close(IN2);
close(OUT);

`mv $outfile2 $outfile`;



} # end bestfit ---------------------------------------------------------------





# Tk widget functions ++++++++++++++++++++++++++++++++++

# write postscript and copy to JPG
sub write_ps_to_jpg {

my(@data) = @_;
my($ps_name) = $data[0];

my($ps_out) = "$ps_name.ps";

$canvas -> postscript(
        -file, "$ps_out",
        -rotate,1,
        -pageheight, '500',
        -pagewidth, '750'
        );

print "\n\nwrote $ps_out\n\n";

my($jpg_out) = "$ps_name.png";
my($dens) = 600;

print "convert $ps_out\nto      $jpg_out\n";

print "convert -density $dens $ps_out -quality 100 $jpg_out\n";
`convert -density $dens $ps_out -quality 100 $jpg_out`;

print "rotate jpg\n";
my($jpg_out2) = "$ps_name-1.jpg";
my($rot) = 90;

`convert $jpg_out -rotate $rot $jpg_out2`;
print "\nconvert $jpg_out -rotate $rot $jpg_out2\n";

print "\nwrote $jpg_out2\n\n";

`rm $jpg_out`;


} # end save postscript -----------------------------------





# make text ++++++++++++++++++++++++++++++++++++++++++++

sub TkText {


my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($text) = $data[2];
my($anchor) = $data[3];

my($color) = "black";
if (exists $data[4]) {
        $color = $data[4];
        }

my($font) = "arial {8}";
if (exists $data[5]) {
        $font = $data[5];
        }

$canvas -> createText(
        $x,$y,
        -text, $text,
	-font, $font,
        -anchor, $anchor,
	-fill, $color,
        );
}


# make Box ++++++++++++++++++++++++++++++++++++++++++++

sub TkBox {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];
my($color) = $data[4];

my($outline) = "black";
if (exists $data[5]) {
        $outline = $data[5];
        }

my($wide) = 1;
if (exists $data[6]) {
        $wide = $data[6];
        }

$canvas -> createRectangle(
        $x,$y,$x2,$y2,
        -fill, $color,
	-outline, $outline,
	-width, $wide,
        );
}

# make empty Box with colored outline ++++++++++++++++++++++++++++++++++++++++++++

sub TkBoxEmpty {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];
my($outline) = $data[4];

my($wide) = 1;
if (exists $data[5]) {
        $wide = $data[5];
        }

$canvas -> createRectangle(
        $x,$y,$x2,$y2,
        -outline, $outline,
        -width, $wide,
        );
}





# make Polygon with 4 corners ++++++++++++++++++++++++++++++++++++++++++++

sub TkPoly {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];
my($x3) = $data[4];
my($y3) = $data[5];
my($x4) = $data[6];
my($y4) = $data[7];
my($color) = $data[8];

my($outline) = $color;
if (exists $data[9]) {
        $outline = $data[9];
        }

$canvas -> createPolygon(
        $x,$y,$x2,$y2,$x3,$y3,$x4,$y4,
        -fill, $color,
        -outline, $outline,
        );
}


# make Polygon with 3 corners ++++++++++++++++++++++++++++++++++++++++++++

sub TkPoly3 {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];
my($x3) = $data[4];
my($y3) = $data[5];
my($color) = $data[6];

my($outline) = $color;
if (exists $data[7]) {
        $outline = $data[7];
        }

$canvas -> createPolygon(
        $x,$y,$x2,$y2,$x3,$y3,
        -fill, $color,
        -outline, $outline,
        );
}



# make Line ++++++++++++++++++++++++++++++++++++++++++++

sub TkLine {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];

my($color) = "black";
if (exists $data[4]) {
        $color = $data[4];
        }

my($wide) = 1;
if (exists $data[5]) {
	$wide = $data[5];
	}


$canvas -> createLine(
        $x,$y,$x2,$y2,
        -fill, $color,
	-width, $wide,
        );
} # end mkae line ---------------------------------------

# make arrow ++++++++++++++++++++++++++++++++++++++++++++

sub TkArrow {

my(@data) = @_;
my($x) = $data[0];
my($y) = $data[1];
my($x2) = $data[2];
my($y2) = $data[3];

my($color) = "black";
if (exists $data[4]) {
        $color = $data[4];
        }

my($wide) = 1;
if (exists $data[5]) {
        $wide = $data[5];
        }


$canvas -> createLine(
        $x,$y,$x2,$y2,
        -fill, $color,
        -width, $wide,
	-arrow, 'last',
        );
} # end mkae line ---------------------------------------




# make heatmap colors RED ++++++++++++++++++++++++++++++++++
# make red color spec
sub mk_red {

my($i);
my($hexval);
my($hexval2);
my($integer);


my(@heat_red) = ();

for($i=0;$i<32;$i++) {

        # calculate color code in hexadec system +++++

        $modulo = $i%16;
        $down = 15-$modulo;
        $integer = int($i/16);
        $hexval = sprintf("%x", $modulo);
        $hexval =~ tr/a-z/A-Z/;
        $hexval2 = sprintf("%x", $down);
        $hexval2 =~ tr/a-z/A-Z/;

        # go from black to red
        if ($integer == 0) {
        $red = "$hexval$hexval";
        $green = "00";
        $blue = "00";
        }

        # go from red to yellow 
        if ($integer == 1) {
        $red = "FF";
        $green = "$hexval$hexval";
        $blue = "00";
        }

        $color = "#$red$green$blue";
        push(@heat_red,$color);

}

return (@heat_red);

}




 





# merge overlap +++++++++++++++++++++++++++++++++++++++++++++
# check if 2 sequences overlap with a minimum given primer
# bad quality matches at both ends are ignores as long as the required primer is found

sub merge_OL {

my(@seq) = @_;
my($seq1) = $seq[0];
my($seq2) = $seq[1];
my($primer) = $seq[2];
my($i);
$seq1 =~ s/-//g;
$seq2 =~ s/-//g;
my($len1) = length($seq1);
my($len2) = length($seq2);
my($merge) = '';
my($rest1_1);
my($rest1_1);
my($rest1_2);
my($rest2_1);
my($rest2_2);
my($len_rest1_1);
my($len_rest1_2);
my($len_rest2_1);
my($len_rest2_2);


MERGE: for($i=0;$i<=($len1-$primer);$i++) {

         $sub = substr($seq1,$i,$primer);

        ($rest2_1,$rest2_2) = ($seq2 =~ /^(.*)$sub(.*)$/);

        # if primer is found +++++++
        if ($rest2_1) {

                ($rest1_1,$rest1_2) = ($seq1 =~ /^(.*)$sub(.*)$/);

                # determine which overjang is longer and create longest possible fusion product
                $len_rest1_1 = length($rest1_1);
                $len_rest1_2 = length($rest1_2);
                $len_rest2_1 = length($rest2_1);
                $len_rest2_2 = length($rest2_2);

                if (($len_rest1_1 >= $len_rest2_1) && ($len_rest1_2 >= $len_rest2_2)) {
                        $merge = "$rest1_1$sub$rest1_2";
                        }

                if (($len_rest2_1 >= $len_rest1_1) && ($len_rest1_2 >= $len_rest2_2)) {
                        $merge = "$rest2_1$sub$rest1_2";
                        }

                if (($len_rest2_1 >= $len_rest1_1) && ($len_rest2_2 >= $len_rest1_2)) {
                        $merge = "$rest2_1$sub$rest2_2";
                        }

                if (($len_rest1_1 >= $len_rest2_1) && ($len_rest2_2 >= $len_rest1_2)) {
                        $merge = "$rest1_1$sub$rest2_2";
                        }

                print "merged\n$seq1\n$seq2\n$merge\n\n";
                last MERGE;
                } # end if primer is found
        }

return($merge);
}# -----------------------------------------------------------------------------------------


# extraxt sequneces from pair file produced by water
sub sequences_from_pair {

my(@pair) = @_;
my($infile) = $pair[0];
my($seq1) = '';
my($seq2) = '';
my($aln);
my($se);
my($be1);
my($en1);
my($be2);
my($en2);
my($be);
my($en);
my($flag) = 1;
my($first1) = 1;
my($first2) = 1;
my($space) = '';
my($len_space) = 0;
my($al) = '';

open(INPAIRFILE,"$infile");
while (<INPAIRFILE>) {

        # extract sequences ++++++++++++++++++++++++++
        ($be,$se,$en) = /^\S+\s+(\d+)\s+(\S+)\s+(\d+)/;
        if ($be) {

                # collect sequence 1 +++++++++
                if ($flag ==1) {
                        $seq1 .= $se;
                        $en1 = $en;
                        if ($first1 ==1) {
                                $be1 = $be;
                                $first1 =0;
				# dtermine laenth of spacer
				($space) = /^(\S+\s+\d+\s+)\S+/;
				$len_space = length($space);
                                }
                        $flag = 2;
                        next;
                        }

                # collect sequence 2 +++++++++
                if ($flag ==2) {
                        $seq2 .= $se;
                        $en2 = $en;
                        if ($first2 ==1) {
                                $be2 = $be;
                                $first2 =0;
                                }
                        $flag = 1;
                        next;
                        }
                }

	# extract alignment string
	($al) = /^\s{$len_space}(.+)/;
	if ($al) {
		if ($first1==0) {
			$aln .= $al;
			}
		}


        }

close(INPAIRFILE);

# return sequnece, beg and end points
return($seq1,$seq2,$be1,$en1,$be2,$en2,$aln);

} # end sequences from pair ------------------------------------- 


# fetch a sequence from a flatfile 
sub fetch_seq {

my($query,$db) = @_;
my($len);
my($go);
my($check);
my($stop);
my($len);
my($found) =0;
my($out) = $query;
my($count);

#print "fetch $query from $db\n";

open (FETCH_SEQ,">$query");
open (FETCH_DB, "$db");
FETCH: while (<FETCH_DB>) {
        ($stop) = /^(>)/;
        if ($stop) {
                $check = 0;
                if ($found ==1) {
                        last FETCH;
                        }
                }

        ($go) = /^>($query)[^\S]/;
	#($go) = /^>($query)/;
        if ($go) {
                $found =1;
                print FETCH_SEQ ">$query\n";
                $check =1;
                $count++;
                }

        if ($check ==1) {
                unless ($_ =~ />/) {
                        print FETCH_SEQ "$_";
                        chomp($_);
                        $len += length($_);
                        }
                }

        }

close(FETCH_DB);
close(FETCH_SEQ);

if ($found ==1 ) {
        #print "$out length = $len\n";
        } else {
        print "$query not found in $db!!!\n";
        `rm $out`;
        }

} # end fetch ---------------------------------------










1

