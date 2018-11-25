#!/usr/bin/perl

# Histogram 2D
# 
# D. Gonze
# 27/2/2004
#
#########################################################

&ReadArguments;

&CreateDataVector;

&Histo;

################################################################
### Read arguments from the command line

sub ReadArguments {

$verbo=0;
$infile = "";
$outfile = "";
$col1=0; # first column
$col2=1; # second column
$min1="";
$min2="";
$max1="";
$max2="";
$int1="";
$int2="";

foreach my $a (0..$#ARGV) {

    ### help
    if ($ARGV[0] eq "-h") {
    	&PrintHelp;
    }
            
    ### input file
    elsif ($ARGV[$a] eq "-i") {
    	$infile = $ARGV[$a+1];
    }

    ### output file 
    elsif ($ARGV[$a] eq "-o") {
	$outfile = $ARGV[$a+1];
    }

    ### column with the data
    elsif ($ARGV[$a] eq "-col") {
    	$col1 = $ARGV[$a+1]-1;
    	$col2 = $ARGV[$a+2]-1;
    }
	
    ### maximums 
    elsif ($ARGV[$a] eq "-max1") {
    	$max1 = $ARGV[$a+1];
    }
    elsif ($ARGV[$a] eq "-max2") {
    	$max2 = $ARGV[$a+1];
    }

    ### minimums
    elsif ($ARGV[$a] eq "-min1") {
    	$min1 = $ARGV[$a+1];
    }
    elsif ($ARGV[$a] eq "-min2") {
    	$min2 = $ARGV[$a+1];
    }

    ### intervals
    elsif ($ARGV[$a] eq "-int1") {
    	$int1 = $ARGV[$a+1];
    }
    elsif ($ARGV[$a] eq "-int2") {
    	$int2 = $ARGV[$a+1];
    }

#    ### frequency
#    elsif ($ARGV[$a] eq "-freq") {
#    	$freq = 1;
#    }

    ### verbosity 
    elsif ($ARGV[$a] eq "-v") {
	$verbo=1;
    }
	
}
	
if ($infile eq "") {
   die "STOP! You have to give the name of the input file!\n";
}

if ($outfile eq "") {
   $outfile = "$infile.histo";
}

} # End of ReadArgument


##########################################################################################
### Print help

sub PrintHelp {
  open HELP, "| more";
  print <<EndHelp;
NAME
        histogram-2d.pl

DESCRIPTION
	Prepare data for 2D-histrogram.

AUTHOR
	Didier Gonze (dgonze\@ucmb.ulb.ac.be)  

UPDATED
	26/2/2004

OPTIONS
	-i in_file_name
		Specify the input file containing the data. 
		This argument is obligatory (except if using option -h).
	       
	-o out_file_name
		Specify the output file. Default name is in_file_name.histo

	-col # #
		Specify the columns containing the data. Default: 1 2.

	-max1 #
		Specify the max of x-data for the histogram.
		Default: max of x-data.

	-max2 #
		Specify the max of y-data for the histogram.
		Default: max of y-data.
		
	-min1 #
		Specify the min of x-data for the histogram.
		Default: min of x-data.

	-min2 #
		Specify the min of y-data for the histogram.
		Default: min of y-data.

	-int1 #
		Specify the interval of x-data for the histogram. 
		Default=(max1-min1)/10.

	-int2 #
		Specify the interval of y-data for the histogram. 
		Default=(max2-min2)/10.
		
#	-freq
#		To have the frequency instead of the occurrence.
	
	-v 
		Verbosity: print detailed informations during the 
		process.
	       
	-h 
		Give help (print this message). This argument must be 
		the first.

EXAMPLE
        perl histogram-2d.pl -i datafile -col 2 3 -o results.out

EndHelp
close HELP;
die "\n";

} # Enf of Help


##########################################################################################
### Read the data and fill the data vector

sub CreateDataVector {

my $i=0;

open inf, $infile or die "STOP! File $infile not found";
if ($verbo==1) {print "Open input file: $infile\n";}

$dmin1=99999999;
$dmax1=-99999999;
$dmin2=99999999;
$dmax2=-99999999;

foreach $line (<inf>){
   if ($line !~ /^#/ and $line ne ""){
      chomp $line;
      @line=split /\t/,$line;      
      $x[$i]=$line[$col1];
      $y[$i]=$line[$col2];
      
      if ($x[$i] <= $dmin1) {$dmin1=$x[$i]}
      elsif ($x[$i] >= $dmax1) {$dmax1=$x[$i]}
      if ($y[$i] <= $dmin2) {$dmin2=$y[$i]}
      elsif ($y[$i] >= $dmax2) {$dmax2=$y[$i]}

      $i=$i+1;
   }
}

if ($min1 eq ""){$min1=$dmin1;}    # automatic determination of min1
if ($max1 eq ""){$max1=$dmax1;}    # automatic determination of max1
if ($min2 eq ""){$min2=$dmin2;}    # automatic determination of min2
if ($max2 eq ""){$max2=$dmax2;}    # automatic determination of max2
if ($int1 eq ""){$int1=($max1-$min1)/10;}    # automatic determination of int1
if ($int2 eq ""){$int2=($max2-$min2)/10;}    # automatic determination of int2

$nx=$i;

if ($verbo==1) {print "Total number of data = $nx\n";}

close inf;

} # End of CreateDataVector

##########################################################################################
### Histo

sub Histo {

my ($i,$j1,$j2)=(0,0);
my $jmax1, $jmax2;
my @z=0;

open ouf, ">$outfile";
if ($verbo==1) {print "Creating output file: $outfile\n";}

if ($verbo==1) {print "Generate data for 2D-histogram...\n";}

$jmax1=($max1-$min1)/$int1;
$jmax2=($max2-$min2)/$int2;

if ($verbo==1) {print "For x-data: min=$min1, max=$max1, int=$int1 ($jmax1 steps)\n";}
if ($verbo==1) {print "For y-data: min=$min2, max=$max2, int=$int2 ($jmax2 steps)\n";}

for ($i=0;$i<=$nx-1;$i++){
  for ($j1=0;$j1<=$jmax1;$j1++){
    $cmin1=$min1+$j1*$int1;
    $cmax1=$cmin1+$int1;
  for ($j2=0;$j2<=$jmax2;$j2++){
    $cmin2=$min2+$j2*$int2;
    $cmax2=$cmin2+$int2;
    if(($x[$i]>=$cmin1) and ($x[$i]<$cmax1) and ($y[$i]>=$cmin2) and ($y[$i]<$cmax2)){
      $z[$j1][$j2]=$z[$j1][$j2]+1;
    }
  }
  }
}

print ouf "#";
for ($j1=0;$j1<=$jmax1;$j1++){
   $cmin1=$min1+$j1*$int1;
   $cmin1=sprintf("%.2f",$cmin1);
   print ouf "\t$cmin1";
}
print ouf "\n";

for ($j2=0;$j2<=$jmax2;$j2++){
   $cmin2=$min2+$j2*$int2;
   $cmin2=sprintf("%.2f",$cmin2);
   print ouf "$cmin2";
for ($j1=0;$j1<=$jmax1;$j1++){
   $zz=sprintf("%.0f",$z[$j1][$j2]);
   print ouf "\t$zz";
}
print ouf "\n";
}

close ouf;

} # End of Histo

