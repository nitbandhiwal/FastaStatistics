#!/usr/bin/perl


my $version_number = "1.0";
my $help = help_message($version_number);

# if there are no arguments, return the help message and quit
unless($ARGV[0]) {
    die "$help";
}
# if -h is arguments, return the help message and quit
if($ARGV[0] eq '-h')
		{ die "$help"; }

$fasta_file=pop @ARGV;

my $n_arg = scalar(@ARGV);
foreach (my $i = 0;$i < $n_arg;$i++) {
    	my $arg = $ARGV[$i];
	if($arg eq '-n')
		{my $number = number_seq($fasta_file); 
		 print "number of sequences in fasta file : $number\n";
		}
	  if($arg eq '-l')
		{ length_seq($fasta_file); }
	  if($arg eq '-GC')
		{ gc_seq($fasta_file); }
	  if($arg eq '-BP')
		{ my $base_pair = base_pair_seq($fasta_file);
		  print "Total number of base pair in fasta file : $base_pair\n";
		}
	  if($arg eq '-P')
		{ pseudogenome($fasta_file); } 
	  if($arg eq '-F')
		{ $min_len = $ARGV[$i+1];
		filter_seq($fasta_file,$min_len);
		}
	  if($arg eq '-T')
		{ open F1,"$fasta_file" or die $!;
		  @arr = <F1>;
		  if($arr[1] =~ m/M/)
			{ print "Translation can not be done. Sequences are Protein already.\n"; }
		  else	{ translation($fasta_file); }
		}
	  if($arg eq '-CP')
		{ updated_coordinates($fasta_file); }
	  if($arg eq '-N50')
		{ N50($fasta_file); }
	}

####################
# Here be sub-routines

sub help_message {
    my($version) = @_;
    my $message = "
Fasta_stats.pl : Calculating statistics of fasta file and give some related file like Pseudogenome, 6 frame translation, updated coordinates of new pseudogenome etc.

Version: $version

Usage: Fasta_stats.pl \[options\] fasta_file > output_file

Options:
-h  	    : Print help message and quit
-n  	    : Prints number of sequences in fasta file.
-BP 	    : Total base pair in Pseudogenome.
-l  	    : Prints length of every sequences of fasta file in output file.
-T  	    : 6 frame translation (in case fasta file is nucleotide.)
-P  	    : Combine all sequences (make Pseudogenome).
-CP 	    : Updated coordiantes according to Pseudogenome.
-GC 	    : Prints GC content of every sequences in fasta file.
-F [number] : Discard sequences which have sequence length below F.
-N50        : Print assembly statistics such as N50, Number of Sequneces, total base pair, largest and smallest contig 

";
    return $message;
}

sub N50 {
	my ($fasta_file) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{			
		push @numbers,$seq_length;
		}
	use List::Util qw/sum min max/;
	my $min_C = min @numbers;
	my $max_C = max @numbers;
	close($fh);
	$minlen ||= 0;
	$n__    ||= 50;

	my @len = ();
	open SEQ, "<", $fasta_file or die "Cannot open file: $fasta_file: $!\n";
	while(<SEQ>){
   		if(/^>/){
      			push @len, 0;
   				}else{
      				next if /^;/;
      				chomp;
      				s/\W//g;
      				$len[-1]+=length $_;
   				}
		  }
	close SEQ;
	@len = sort { $a <=> $b } map { $_>=$minlen?$_:() } @len;
	my $tot = (sum(@len) || 0);
	my $thr = $n__*$tot/100;
	my $pos = 0;
	for(@len){
		$pos+= $_;
		if($pos>=$thr){
      		#print "N$n__: $_\n";
		$N50="$_";
      		last;
   		}
		}

		#print "Sequences: ".scalar(@len)."\n";
		#print "Total length: $tot\n";
		print "N50\tNumberOfSeq\tTotalLength\tMinContigLength\tLargestContigLength\n";
		print "$N50\t".scalar(@len)."\t$tot\t$min_C\t$max_C\n";


}

sub  number_seq {
	my ($fasta_file) = @_;
	open F1,"$fasta_file" or die $!;
	$num_count =0;
	while(<F1>)
		{ if($_ =~ m/>/)
			{ $num_count++; }
		}
	close(F1);
   return $num_count;
}

sub length_seq {
	my ($fasta_file) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{			
		print "$sequence_data{header}\t$seq_length\n";
		}
	close($fh);
}

sub updated_coordinates {
my ($fasta_file) = @_;
	open($fh,$fasta_file) or die "can't open $fasta_file: $!\n";
	open (F1,">tmp_length_file") or die $!;
	my %sequence_data;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{			
		print F1 "$sequence_data{header}\t$seq_length\n";
		}
	close($fh);
	close(F1);
	open (F2,"tmp_length_file") or die "can't open tmp_length_file: $!\n";
	$start=1;
	$count=0;
	$end=0;
	while(<F2>)
        { $count++;
          chomp; @arr = split /\t|\s+/,$_; 
          $end=$start+$arr[1];
          if($count == 1){
        print "Pseudogenome\t$start\t$arr[1]\t$arr[0]\n"; $start = 1+$arr[1];
                }
          else{
        print "Pseudogenome\t$start\t$end\t$arr[0]\n"; $start = $end+1;
                }
          
        }
close(F2);
system("rm -rf tmp_length_file");
}

sub gc_seq {
	my ($fasta_file) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data;
	my $A=$C=$T=$G=0;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{
		$A=$sequence_data{seq} =~ s/[Aa]/[Aa]/g;
		$C=$sequence_data{seq} =~ s/[Cc]/[Cc]/g;
		$T=$sequence_data{seq} =~ s/[Tt]/[Tt]/g;
		$G=$sequence_data{seq} =~ s/[Gg]/[Gg]/g;
		$cnt = $G+$C;
		$total = $A+$T+$C+$G;
		$gc= ($cnt/$total)*100;
		print "$sequence_data{header}\t$seq_length\t$gc\n";
		}
	close($fh);
}

sub base_pair_seq {
	my ($fasta_file) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data; my $seq='';
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{
		$seq=$seq.$sequence_data{seq};
		}
		$len=length($seq);
	close($fh);
	return $len;
}

sub pseudogenome {
	my ($fasta_file) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data; my $seq=''; @x=split /\_/,$fasta_file; $x[-1]=~s/\.fa//g;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{
		$seq=$seq.$sequence_data{seq};
		}
		$len=length($seq);
		print ">$x[-1]\n";
		my @arr2 = unpack("(A60)*", $seq);
	  		foreach $j(@arr2)
				{ print "$j\n"; }
		close($fh);
}

sub filter_seq {
my ($fasta_file,$min_len) = @_;
	open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
	my %sequence_data;
  	while (read_fasta_sequence($fh, \%sequence_data)) 
		{
		if( $seq_length >= $min_len)
			{ print ">$sequence_data{header} $seq_length\n";
			  my @arr2 = unpack("(A60)*", $sequence_data{seq});
	  		  foreach $j(@arr2)
				{ print "$j\n"; }
			}
		}
		close($fh);
}

sub translation {
my ($fasta_file) = @_;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
my $seq='';
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G');
while(<$fh>)
	{ $seq=$seq.$_;}
my @arr = split />/,$seq;
foreach $i (@arr)
	{  $f='';my @arr1 = split /\n/,$i;
	  $name = shift @arr1;	#contain name of contig or transcript
	  $DNA='';
	  foreach $p(@arr1){ $DNA=$DNA.$p;}
	  #$DNA = shift @arr1;	#contain DNA seq
	  $DNA1 = $DNA;
my $protein1=$protein2=$protein3=$protein4=$protein5=$protein6='';
# for first reading frame####
@arr2 = unpack("(A3)*", $DNA);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein1.=$g{$key}; }
			}
		}
### for second reding frame####
$s=reverse($DNA);
chop($s);
$DNA=reverse($s);
@arr2 = unpack("(A3)*", $DNA);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein2.=$g{$key}; }
			}
		}
#### for third reading frame###
$s=reverse($DNA);
chop($s);
$DNA=reverse($s);
@arr2 = unpack("(A3)*", $DNA);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein3.=$g{$key}; }
			}
		}
#### for fourth reading frame####
$revdna=reverse($DNA1);
$revdna =~ tr/ACGTacgt/TGCAtgca/;
@arr2 = unpack("(A3)*", $revdna);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein4.=$g{$key}; }
			}
		}
#### for fifth reading frame####
$s=reverse($revdna);
chop($s);
$revdna=reverse($s);
@arr2 = unpack("(A3)*", $revdna);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein5.=$g{$key}; }
			}
		}
##### for sixth reading frame#####
$s=reverse($revdna);
chop($s);
$revdna=reverse($s);
@arr2 = unpack("(A3)*", $revdna);
	foreach $j(@arr2)
		{ foreach $key(keys %g)
			{ if($key eq $j)
				{ $protein6.=$g{$key}; }
			}
		}

#### comparison for selecting largest ######
my @frag1 = split /_/,$protein1;
my @frag2 = split /_/,$protein2;
my @frag3 = split /_/,$protein3;
my @frag4 = split /_/,$protein4;
my @frag5 = split /_/,$protein5;
my @frag6 = split /_/,$protein6;
$len=0; foreach $q(@frag1) { $len1=length($q); if($len1 >$len){ $len=$len1; $f1=$q;} } push(@final,$f1); # $f1 contain largest from 1st reading frame...
$len=0; foreach $q(@frag2) { $len1=length($q); if($len1 >$len){ $len=$len1; $f2=$q;} } push(@final,$f2); # $f2 contain largest from 2nd reading frame...
$len=0; foreach $q(@frag3) { $len1=length($q); if($len1 >$len){ $len=$len1; $f3=$q;} } push(@final,$f3); # $f3 contain largest from 3rd reading frame...
$len=0; foreach $q(@frag4) { $len1=length($q); if($len1 >$len){ $len=$len1; $f4=$q;} } push(@final,$f4); # $f4 contain largest from 4th reading frame...
$len=0; foreach $q(@frag5) { $len1=length($q); if($len1 >$len){ $len=$len1; $f5=$q;} } push(@final,$f5); # $f5 contain largest from 5th reading frame...
$len=0; foreach $q(@frag6) { $len1=length($q); if($len1 >$len){ $len=$len1; $f6=$q;} } push(@final,$f6); # $f6 contain largest from 6th reading frame...
$len=0; foreach $k(@final) { $len1=length($k); if($len1 >$len){ $len=$len1; $f=$k;} } # $f contain largest from all reading frame...

my $i=0;
my $flag=0;
my $out_trans_seq="";
my @splited= split //,$f;
foreach (@splited)
{
	if( $splited[$i] eq "M" && $flag==0)
	{
	$flag=1;
	}
	if($flag==1)
	{
	$out_trans_seq="$out_trans_seq"."$splited[$i]";
	}
$i++;
}
if($out_trans_seq ne "")
{
print ">$name\n$out_trans_seq\n";
}
### to empty final array ####
#undef @final;
$#final = -1;
}
close($fh);
}


sub read_fasta_sequence {
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/^>//;  
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         $seq_info->{seq} .= $_;

		#find length of each sequence
	 $seq_length = length($seq_info->{seq});
	 
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}

