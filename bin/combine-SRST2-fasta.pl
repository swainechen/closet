#!/usr/bin/perl -w
#
# combine fasta files formatted for SRST2
# just make sure clusterUniqueIdentifier is unique
#
use File::Basename;
use Getopt::Long;
Getopt::Long::Configure("pass_through");
my $SRST2 = "/usr/local/lib/SRST2";
my $VERSION = "0.2.0";
GetOptions (
  'srst2=s' => \$SRST2,
  'version=s' => \$VERSION
);
if (!-f "combine.db") {
  print "Please make sure there is a combine.db file in the current directory\n";
  print "This file should specify all the SRST2 databases to combine\n";
  print "DB name in the final combined file will be based on the individual filenames\n";
  print "Output goes to STDOUT\n";
  exit;
}
my $id = ();
my $i = 1;
my @f;
my $file;
my $label;
open DB, "combine.db";
while ($file = <DB>) {
  next if $file =~ /^#/;
  next if $file =~ /^$/;
  chomp $file;
  $file =~ s/\$SRST2/$SRST2/g;
  $file =~ s/\$VERSION/$VERSION/g;
  if (!-f $file) {
    print STDERR "Can't find file $file in combine.db - ignoring...\n";
    next;
  }
  $label = File::Basename::basename($file);
  # get rid of hopefully standard file extensions
  $label =~ s/\.fa$//;
  $label =~ s/\.fasta$//;
  $label =~ s/\.fna$//;
  $label =~ s/\.ffn$//;
  if ($file =~ /\/VFDB\//) {
    $label =~ s/\.fsa$/_VFDB/;
  }
  open F, $file;
  while (<F>) {
    next if /^#/;
    next if /^$/;
    chomp;
    if (/^>/) {
      @f = split / /, $_;
      if (defined $f[1]) {
        $f[1] = "DB:$label " . $f[1];
      } else {
        $f[1] = "DB:$label";
      }
      # clean up - some of the files have __ separators which make SRST2 choke
      # also if there is a / in the name, SRST2 chokes on making output files
#      @g = split /__/, $f[0];
      $f[0] =~ s/\//_/g;
      $f[0] =~ s/\'/-/g;
      $f[0] =~ s/\"/-/g;
      $f[0] =~ s/\./-/g;
      @h = split /__+/, $f[0];
      if ($#h > 3) {
        $g[0] = $h[0];
        $g[1] = join ("_", @h[1..$#h-2]);
        $g[2] = $h[$#h-1];
        $g[3] = $h[$#h];
      } else {
        @g = @h;
      }
      $g[0] =~ s/^>//;
      if (!defined $id->{$file}->{$g[0]}) {
        $id->{$file}->{$g[0]} = $i;
        $i++;
      }
      $g[0] = ">" . $id->{$file}->{$g[0]};
      print join (" ", join ("__", @g), @f[1..$#f]), "\n";
    } else {
      print uc($_), "\n";
    }
  }
  close F;
}
close DB;
