#!/usr/bin/perl -w
#
# take dna fasta file
# run phylip, then paml
#
if (!@ARGV || !-f $ARGV[0]) {
  print "Usage: $0 <dna fasta file>\n";
  exit;
}

my $file = $ARGV[0];
my $root = 0;
my $kappa = 0;

use slchen;

open F, $file;
@dna = <F>;
close F;

my @tree = slchen::do_ml_tree($root, $kappa, @dna);

open T, ">tree.paml";
print T $tree[0], "\n";
close T;

system "fasta2phy.pl $file | phy2paml > dna.paml";
&reset_paml_control;

$M1 =~ s/__OUTFILE__/$file.M1/g;
open C, ">M1.ctl";
print C $M1;
close C;
system "codeml M1.ctl";
unlink "M1.ctl";
system "mv rst rst.M1";

$M2 =~ s/__OUTFILE__/$file.M2/g;
open C, ">M2.ctl";
print C $M2;
close C;
system "codeml M2.ctl";
unlink "M2.ctl";
system "mv rst rst.M2";

$bsA =~ s/__OUTFILE__/$file.bsA/g;
open C, ">bs-A.ctl";
print C $bsA;
close C;
system "slctree.pl tree.paml";
system "convert out.svg /home/slchen/public_html/out.png";
unlink "out.svg";
print "Edit tree.paml to put in foreground branches\n";
print "Tree is at http://hultgren.wustl.edu/~slchen/out.png\n";
print "Press Enter when done: \n";
#<STDIN>;
#system "codeml bs-A.ctl";
#unlink "bs-A.ctl";

foreach $temp qw(2NG.dN 2NG.dS 2NG.t 4fold.nuc lnf infile rst1 rub) {
  unlink $temp;
}

sub reset_paml_control {

$M1 = 'seqfile = dna.paml
treefile = tree.paml
outfile = __OUTFILE__
noisy = 0
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = /usr/local/src/paml4/dat/wag.dat
model = 0
NSsites = 1
icode = 0
fix_kappa = 0
kappa = 3
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
clock = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 1
';

$M2 = 'seqfile = dna.paml
treefile = tree.paml
outfile = __OUTFILE__
noisy = 0
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = /usr/local/src/paml4/dat/wag.dat
model = 0
NSsites = 2
icode = 0
fix_kappa = 0
kappa = 3
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
clock = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 1
';

$bsA = 'seqfile = dna.paml
treefile = tree.paml
outfile = __OUTFILE__
noisy = 0
verbose = 1
runmode = 0
seqtype = 1
CodonFreq = 2
aaDist = 0
aaRatefile = /usr/local/src/paml4/dat/wag.dat
model = 2
NSsites = 2
icode = 0
fix_kappa = 0
kappa = 3
fix_omega = 0
omega = 1
fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 10
clock = 0
getSE = 0
RateAncestor = 0
Small_Diff = .5e-6
cleandata = 1
';
}
