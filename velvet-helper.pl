#!/usr/bin/perl -w
#
# take a bunch of files
# do some QC on scarf files
# run velvet tests on fasta files
# put images into database
# if we're already there, maybe do the assemblies
#
use DBI;
use File::Temp;
use File::Basename;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");
my $file;
my $fasta;
my @id;
my $debug = 0;
my $tempdir = "";
my $database = "";
my $dbhost = "";
my $dbuser = "";
my $dbpass = "";
my $cwd = `pwd`;
chomp $cwd;
my $dbh;
my $primary_key;
my $person = "";
my $bestk;
my $bestcov;
my $max_cov;	# we'll filter out high coverage in the final assembly
my $min_contig;
my $assembly;
my $cluster = 0;
my $do_velvet = 1;	# whether or not to do velvet tests
my $filter_high_coverage = 1;	# whether or not to filter high coverage
				# (repeats) in final assembly, recommended for
				# SNP detection
my $command_switches = "";
my $velveth = "/home/comp/jglab/slchen/bin/velveth";
if (!-x $velveth) { $velveth = `which velveth`; chomp $velveth; }
my $velvetg = "/home/comp/jglab/slchen/bin/velvetg";
if (!-x $velvetg) { $velvetg = `which velvetg`; chomp $velvetg; }

GetOptions (
  "person=s" => \$person,
  "cluster!" => \$cluster,
  "velvet!" => \$do_velvet,
  "maxcoverage!" => \$filter_high_coverage,
  "debug!" => \$debug,
  "tempdir=s" => \$tempdir
);

if ($do_velvet) {
  $command_switches .= " -velvet";
} else {
  $command_switches .= " -novelvet";
}
if ($filter_high_coverage) {
  $command_switches .= " -maxcoverage";
} else {
  $command_switches .= " -nomaxcoverage";
}

sub connect {
  $dbh = DBI->connect('DBI:mysql:database='.$database.';host='.$dbhost, $dbuser, $dbpass);
}

sub disconnect {
  $dbh->disconnect;
}

# files are on the command line
if (!@ARGV || !length $person) {
  print "Usage: $0 [ -cluster ] [ -velvet|novelvet ] [ -nomaxcoverage ] [ -tempdir <temp directory> ] -person <your name> <scarf or fasta files>\n";
  print "Make sure you quote your name if you use a space\n";
  print "Takes new or old scarf files and does a little QC\n";
  print "If -cluster is set, print out commands ready for nq\n";
  print "If -novelvet, only do quality graphs for scarf file\n";
  print "If -nomaxcoverage, don't filter out high coverage contigs, useful if you want to get the largest assemblies including duplicated regions and don't care about SNPs, this only applies to the final assembly\n";
  print "If you specify a tempdir make sure it already exists, default is the system temp directory\n";
  print "Scarf and fasta files get run through a bunch of velvet parameters\n";
  print "Data gets put into a database where you can see them from the web\n";
  print "Web site is http://jig-slchen.wustl.edu/velvet-helper\n";
  print "Check assembly parameters and enter best K and coverage there\n";
  print "Then you can run $0 on your fasta file and it will give you a final assembly\n";
  exit;
}

if ($debug) {
  if (!length($tempdir) || !-d $tempdir) {
    $tempdir = File::Temp::tempdir( CLEANUP => 0 );
  }
  print STDERR "Tempdir: $tempdir\n";
  print STDERR "Host: ", `hostname`;
} else {
  if (!length($tempdir) || !-d $tempdir) {
    $tempdir = File::Temp::tempdir( CLEANUP => 1 );
  }
}

foreach $file (@ARGV) {
  next if !-f $file;
  if ($cluster) {
    print "$0 -nocluster $command_switches -person $person $file\n";
    next;
  }
  print "File: $file\n";
  if (test_fasta($file)) {
    $fasta = $file;
  } else {
    $fasta = do_scarf($file);
  }
  @id = fingerprint($fasta);
  ($primary_key, $bestk, $bestcov, $max_cov, $min_contig) = insert_data($person, @id);

  if ($do_velvet) {
    chdir $tempdir;
    my $fullname = $fasta;
    if (-f "$cwd/$fullname") {
      # hope we don't get in here if we're given a full path for the file
      $fullname = "$cwd/$fullname";
    }
    if ($bestk && $bestcov) {
      # just do the assembly
      $assembly = do_final_assembly($fullname, $bestk, $bestcov, $max_cov, $min_contig, $filter_high_coverage);
      $fullname =~ s/\.fasta$//;
      &insert_pngs($tempdir, $primary_key);
      if ($max_cov && $filter_high_coverage) {
        system "mv $tempdir/$assembly $fullname.assembled-$bestk-$bestcov-$max_cov";
      } else {
        system "mv $tempdir/$assembly $fullname.assembled-$bestk-$bestcov";
      }
    } else {
      &do_velvet_tests($fullname, $primary_key);
      &insert_pngs($tempdir, $primary_key);
    }
    chdir $cwd;
  } else {
    &insert_pngs($tempdir, $primary_key);
  }

}

sub fingerprint {
  my ($file) = @_;
  # give some identifying information so we can tell if we see the same file
  # assume we're getting a fasta file so ignore headers
  # use name, length, date, MD5sum, and SHAsum
  my $MD5 = `head -2000 $file | grep -v '^>' | md5sum`;
  $MD5 = substr($MD5, 0, 32);
  my $SHA = `head -2000 $file | grep -v '^>' | sha1sum`;
  $SHA = substr($SHA, 0, 40);
  my $numreads = `wc -l $file`;
  $numreads = (split /\s/, $numreads)[0] / 2;
  my @s = stat($file);
  my $length = $s[7];
  my $date = scalar(localtime($s[9]));
  my $name = File::Basename::basename($file);
  return ($name, $length, $numreads, $date, $MD5, $SHA);
}

sub test_fasta {
  my ($file) = @_;
  my $test_limit = 10;
  my $linecounter = 0;
  my $fasta_score = 0;
  open F, $file;
  while (<F>) {
    if (/^>/) { $fasta_score++; }
    $linecounter++;
    last if $linecounter >= $test_limit;
  }
  close F;
  if ($fasta_score > 3) {
    return 1;
  } else {
    return 0;
  }
}

sub do_scarf {
  my ($file) = @_;
  my $fasta = $file . ".fasta";
  # try to only run through once, so keep track of histograms we want
  # quality by nt, quality by position, number of n's per read
  # we're also going to make a fasta file
  my $nthist;	# histogram of q by nt, no N's
  my $poshist;	# histogram of q by position, no N's
  my $nhist;	# histogram of # n per read
  my @avgq;	# average qualities of all reads
  my @goodq;	# qualities of non-N nts
  my @nt;
  my @q;
  my $i;
  my $type;
  my @f;
  my $numn;
  my $numreads = 0;
  my $notify = 0;

  open F, $file;
  open FASTA, ">$fasta" || die "Can't make file $fasta for $file\n";;
  while (<F>) {
    chomp;
    next if /^#/ || /^$/;
    $numreads++;
    @nt = ();
    @q = ();
    @f = split /:/, $_;
    print FASTA ">", join (":", @f[0..4]), "\n$f[5]\n";
    $f[6] =~ s/^\s+//;
#    if (!defined $type) {
      if (length($f[6]) == length($f[5])) {
        $type = 1;	# new style
      } else {
        $type = 0;	# old style
      }
#    }
    $f[5] = uc $f[5];
    $numn = $f[5] =~ tr/N/N/;
    $nhist->{$numn}++;
    @nt = split //, $f[5];
    if ($type) {
      @q = split //, $f[6];
      foreach $i (0..$#q) {
        $q[$i] = ord($q[$i]) - 64;
      }
    } else {
      $f[6] =~ s/^\s+//;
      @q = split /\s+/, $f[6];
    }
    @goodq = ();
    foreach $i (0..$#nt) {
      next if $nt[$i] eq "N";
      push @goodq, $q[$i];
      $nthist->{$nt[$i]}->{$q[$i]}++;
      $poshist->{$i}->{$q[$i]}++;
    }
    push @avgq, arraymean(@goodq);
    $notify++;
    print "$notify\n" if $notify/100000 == int($notify/100000);
  }
  close FASTA;
  close F;

  # draw the graphs
  gri_histogram("Number_N_Per_Read_Hist", "Number N per Read", "Number Reads, total $numreads", $nhist);

  my $posbw;
  foreach $i (sort { $a <=> $b } keys %$poshist) {
    @{$posbw->{$i}} = hist2box($poshist->{$i});
  }
  gri_bw("Quality_by_Position", "Read Position", "Quality", $posbw, 1);

  my $ntbw;
  foreach $i (sort { $a cmp $b } keys %$nthist) {
    @{$ntbw->{$i}} = hist2box($nthist->{$i});
  }
  gri_bw("Quality_by_NT", "Nucleotide", "Quality", $ntbw, 0);

  return $fasta;
}

sub arraymean {
  my (@a) = @_;
  my $mean = 0;
  my $a;
  foreach $a (@a) {
    $mean += $a;
  }
  $mean /= scalar @a if scalar @a;
  return $mean;
}

sub get_primary_key {
  # return value of 0 means we didn't find anything
  my @id = @_;
  my $sql = "SELECT VelvetID FROM Velvet WHERE NumReads = ? AND MD5 = ? AND SHA = ?";
  my $sth = $dbh->prepare($sql);
  my $return = 0;
  $sth->execute($id[1], $id[3], $id[4]);
  while (@data = $sth->fetchrow_array()) {
    $return = $data[0];
  }
  return $return;
}

sub insert_data {
  my ($person, $name, $length, $numreads, $date, $MD5, $SHA) = @_;
  my $sql;
  my $sth;
  my @data;
  my $velvetID;
  my $k = '19,21,23,25,27,29,31';
  my $cov = '3..15';
  my $min_contig = 100;
  my $bestk = 0;
  my $bestcov = 0;
  my $max_cov = 0;

  &connect;

  $velvetID = get_primary_key($name, $numreads, $date, $MD5, $SHA);

  if ($velvetID) {
    print "Looks like file $name has been seen before, primary key $velvetID\n";
    # check for updated data
    $sql = "SELECT BestK, BestCoverage, TestMinContig, CoverageCutoff FROM Velvet WHERE VelvetID = ?";
    $sth = $dbh->prepare($sql);
    $sth->execute($velvetID);
    while (@data = $sth->fetchrow_array()) {
      ($bestk, $bestcov, $min_contig, $max_cov) = @data;
    }
  } else {
    # this is a new file
    $sql = "INSERT INTO Velvet (Person, OriginalFileName, OriginalSize, NumReads, OriginalDate, MD5, SHA, TestK, TestCov, TestMinContig) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
    $sth = $dbh->prepare($sql);
    $sth->execute($person, $name, $length, $numreads, $date, $MD5, $SHA, $k, $cov, $min_contig);
    $velvetID = $sth->{mysql_insertid};
  }

  &disconnect;

  return ($velvetID, $bestk, $bestcov, $max_cov, $min_contig);
}

sub insert_pngs {
  my ($dir, $primary_key) = @_;
  my $sql = "REPLACE INTO VelvetPNG (VelvetID, GraphName, PNG) VALUES (?, ?, ?)";
  my $sth;
  my @png;
  my $png;
  my $imgdata;
  my @finfo;
  my $type;

  &connect;
  $sth = $dbh->prepare($sql);

  opendir D, $dir;
  @png = readdir D;
  closedir D;

  foreach $png (@png) {
    if (-f "$dir/$png" && $png =~ /\.png$/i) {
      $type = $png; 
      $type =~ s/\.png$//;
      @finfo = stat("$dir/$png");
      open PNG, "$dir/$png";
      read PNG, $imgdata, $finfo[7];
      close PNG;

      $sth->bind_param(1, $primary_key);
      $sth->bind_param(2, $type);
      $sth->bind_param(3, $imgdata);
      $sth->execute;

      !$debug && unlink "$dir/$png";
    }
  }

  &disconnect;

}

sub hist2box {
  my ($ref) = @_;
  # take ref to discrete histogram, give back values for whiskers and box
  # we assume the histogram is on integers (i.e. quality)
  # ref->value is a count
  # order->count is a value
  my $key;
  my $total = 0;
  my $order;	# keep track of rank order of values
  my ($lw, $lq, $med, $hq, $hw);
  undef $lq;
  undef $med;
  undef $hq;

  foreach $key (sort { $a <=> $b } keys %$ref) {
    $total += $ref->{$key};
    $order->{$total} = $key;
  }
  foreach $key (sort { $a <=> $b } keys %$order) {
    if (!defined $lq && $key >= $total/4) {
      $lq = $order->{$key};
    }
    if (!defined $med && $key >= $total/2) {
      $med = $order->{$key};
    }
    if (!defined $hq && $key >= 3*$total/4) {
      $hq = $order->{$key};
    }
  }
  $lw = $lq - 1.5*($hq - $lq);
  $hw = $hq + 1.5*($hq - $lq);
  # whiskers have to land on a data point
  foreach $key (sort { $a <=> $b } keys %$ref) {
    if ($key >= $lw) {
      $lw = $key;
      last;
    }
  }
  foreach $key (reverse (sort { $a <=> $b } keys %$ref)) {
    if ($key <= $hw) {
      $hw = $key;
      last;
    }
  }
  return ($lw, $lq, $med, $hq, $hw);
}

sub do_velvet_tests {
  my ($file, $velvetID) = @_;
  # we should already be in the temp directory

  # some defaults
  my $dirname = File::Basename::basename($file);
  my $k;
  my $cov;
  my $min_contig;
  my $hist_cutoff = 50;	# just for plotting a coverage histogram
			# actually autoscale this to 2x expected coverage
  my $numreads;
  my $readlength = 36;
  my $assembled_length;
  my $expected_coverage;
  my $Rfile = "commands.R";
  my $testfile = "testdata.txt";
  my @k;
  my @cov;
  my @g;
  my @h;
  my $g;
  my @f;
  my @r;
  my $sql;
  my $sth;
  my @fasta;
  my @data;

  &connect;

  # get the real values from the database
  $sql = "SELECT NumReads, TestK, TestCov, TestMinContig FROM Velvet WHERE VelvetID = ?";
  $sth = $dbh->prepare($sql);
  $sth->execute($velvetID);
  while (@data = $sth->fetchrow_array()) {
    ($numreads, $k, $cov, $min_contig) = @data;
  }

  &disconnect;

  @k = parse_list($k);
  @cov = parse_list($cov);
  @r = ();
  push @r, "library(plotrix)";
  push @r, "library(lattice)";
  push @r, "cov <- data.frame(k = numeric(0), mids = numeric(0), counts = numeric(0))";

  # prep directories and get stats.txt files for coverage analysis
  foreach $k (@k) {
    print "Coverage tests, k = $k\n";
    @h = `$velveth $tempdir/$dirname-$k $k -fasta $file 2>&1`;
    $debug && print "$velveth $tempdir/$dirname-$k $k -fasta $file 2>&1\n";;
    $debug && print @h;
    @g = `$velvetg $tempdir/$dirname-$k 2>&1`;
    $debug && print "$velvetg $tempdir/$dirname-$k 2>&1\n";
    $debug && print @g;
    # calculate expected coverage
    if (!defined $expected_coverage) {
      $assembled_length = `grep -v '^>' $tempdir/$dirname-$k/contigs.fa | wc`;
      $assembled_length =~ s/^\s+//;
      $assembled_length = (split (/\s+/, $assembled_length))[2];
      if ($assembled_length) {
        $expected_coverage = int($numreads*$readlength/$assembled_length);
      } else {
        $expected_coverage = 25;	# default to making histogram up to 50
      }
      $hist_cutoff = int(1.5 * $expected_coverage);
    }
    push @r, "stats$k <- read.table(\"$tempdir/$dirname-$k/stats.txt\", header=TRUE)";
    push @r, "cov$k <- weighted.hist(stats$k\$short1_cov, stats$k\$lgth, breaks=0:$hist_cutoff, plot=FALSE)";
    push @r, "cov <- data.frame(k=c(cov\$k, rep($k, length(cov$k\$mids))),
                             mids=c(cov\$mids, cov$k\$mids),
                           counts=c(cov\$counts, cov$k\$counts))";
    # clean up a little to save disk space
    !$debug && unlink("$tempdir/$dirname-$k/Sequences");
    !$debug && unlink("$tempdir/$dirname-$k/Roadmaps");
    !$debug && unlink("$tempdir/$dirname-$k/PreGraph");
    !$debug && unlink("$tempdir/$dirname-$k/LastGraph");
    !$debug && unlink("$tempdir/$dirname-$k/Log");
    !$debug && unlink("$tempdir/$dirname-$k/Graph");
    !$debug && unlink("$tempdir/$dirname-$k/contigs.fa");
  }

  # draw graph
  push @r, "pdf(file=\"$tempdir/Coverage.pdf\", width=10, height=10)";
  push @r, "xyplot(counts ~ mids, data = cov, groups = k, type = 'o', auto.key = TRUE, main = \"Weighted contig coverage plot separated by k for $file\")";
  push @r, "dev.off()";

  open R, ">$tempdir/$Rfile";
  print R join ("\n", @r), "\n";
  close R;
  system "R CMD BATCH $tempdir/$Rfile > /dev/null";
  !$debug && unlink "$tempdir/$Rfile";

  # run through parameters with coverage
  @r = ();
  open O, ">$tempdir/$testfile";
  print O join ("\t", qw(k coverage nodes n50length max total)), "\n";
  foreach $k (@k) {
    print "Assembly tests, k = $k\n";
    foreach $cov (@cov) {
      @g = `$velvetg $tempdir/$dirname-$k -cov_cutoff $cov -min_contig_lgth $min_contig 2>&1 | grep Final`;
      $debug && print @g;
      foreach $g (@g) {
        chomp $g;
        $g =~ s/,//g;
        @f = split /\s+/, $g;
        print O join ("\t", $k, $cov, @f[3,8,10,12]), "\n";
      }
    }
  }
  close O;

  # make the k-coverage plots
  push @r, "library(lattice)";
  push @r, "assembly <- read.table(\"$tempdir/$testfile\", header = TRUE)";
  push @r, "assembly\$k <- as.factor(assembly\$k)";
  push @r, "assembly\$coverage <- as.factor(assembly\$coverage)";
  push @r, "pdf(file=\"$tempdir/Assembly.pdf\", width=10, height=10)";
  push @r, "xyplot(n50length + total ~ coverage, data = assembly, groups = k, type = 'o', auto.key = TRUE, scales = \"free\", main = \"Assembly data separated by k for $file\")";
  push @r, "dev.off()";

  open R, ">$tempdir/$Rfile";
  print R join ("\n", @r), "\n";
  close R;

  system "R CMD BATCH $tempdir/$Rfile > /dev/null";
  system "convert $tempdir/Coverage.pdf $tempdir/VelvetCoverage.png";
  system "convert $tempdir/Assembly.pdf $tempdir/VelvetAssembly.png";
  !$debug && unlink "$tempdir/$Rfile";

}

sub parse_list {
  my ($list) = @_;
  # assume we're only getting integers
  my @numbers = ();
  my @temp = ();
  chomp $list;
  my @atoms = split /,/, $list;
  foreach my $i (0..$#atoms) {
    if ($atoms[$i] =~ /\.\./) {
      @temp = split /\.\./, $atoms[$i];
      if ($temp[0] < $temp[1]) {
        push @numbers, ($temp[0]..$temp[1]);
      } else {
        push @numbers, reverse ($temp[1]..$temp[0]);
      }
    } elsif ($atoms[$i] =~ /-/) {
      @temp = split /-/, $atoms[$i];
      if ($temp[0] < $temp[1]) {
        push @numbers, ($temp[0]..$temp[1]);
      } else {
        push @numbers, reverse ($temp[1]..$temp[0]);
      }
    } else {
      push @numbers, $atoms[$i];
    }
  }
  return @numbers;
}

sub gri_histogram {
  my ($filename, $xname, $yname, $histref) = @_;
  # draw a simple bar graph
  my $width = 0.35;
  my $total = 0;
  my $x;
  my $xl;	# x left
  my $xr;	# x right
  my $yl;	# y low
  my $yh;	# y high
  open GRI, "| gri -nowarn_offpage -batch -output $tempdir/$filename.ps 2>&1 > /dev/null";
  print GRI "set x name \"$xname\"\n";
  print GRI "set y name \"$yname\"\n";
  print GRI "read columns x y\n";
  foreach $x (keys %$histref) {
    $total += $histref->{$x};
    next if $x == 0;
    print GRI "$x $histref->{$x}\n";
    print GRI $x-1, " $histref->{$x}\n";
    print GRI $x+1, " $histref->{$x}\n";
  }
  print GRI "\n";
  print GRI "draw axes\n";
  foreach $x (keys %$histref) {
    next if $x == 0;
    $xl = $x - $width;
    $xr = $x + $width;
    $yl = 0;
    $yh = $histref->{$x};
    print GRI "draw box $xl $yl $xr $yh\n";
  }
  close GRI;

  system "convert -density 200 $tempdir/$filename.ps $tempdir/$filename.png";
  !$debug && unlink "$tempdir/$filename.ps";
}

sub gri_bw {
  my ($filename, $xname, $yname, $bwref, $numeric) = @_;
  # use gri to draw a box and whisker plot
  # take box and whisker data as a reference
  # if $numeric, then treat bwref keys as numbers for x axis
  my $xmax = scalar keys %$bwref;
  my $ymin;
  my $ymax;
  my $x;
  my @label;
  my $i;
  if ($numeric) {
    @label = sort {$a <=> $b} keys %$bwref;
  } else {
    @label = sort {$a cmp $b} keys %$bwref;
  }
  my $width = 0.25;
  my @temp;
  my $temp;
  foreach $x (@label) {
    @temp = sort {$a <=> $b} @{$bwref->{$x}};
    $temp = pop @temp;
    $ymax = $temp if !defined $ymax || $ymax < $temp;
    $temp = shift @temp;
    $ymin = $temp if !defined $ymin || $ymin > $temp;
  }

  # the drawing
  open GRI, "| gri -nowarn_offpage -batch -output $tempdir/$filename.ps 2>&1 > /dev/null";

  # set up axes
  print GRI "set font size 5\n";
  print GRI "set x axis 0 ", scalar(@label) + 1, " 1\n";
  print GRI "set x axis labels ";
  foreach $i (0..$#label) {
    if ($numeric) {
      print GRI $i+1, " \"", $label[$i]+1, "\" ";
    } else {
      print GRI $i+1, " \"$label[$i]\" ";
    }
  }
  print GRI "\n";
  
  print GRI "set y axis unknown\n";
  print GRI "read columns x y\n";
  print GRI "0 $ymin\n";
  print GRI scalar(@label), " $ymax\n";
  print GRI "\n";
  print GRI "set x name \"$xname\"\n";
  print GRI "set y name \"$yname\"\n";
  print GRI "draw axes\n";

  # the real drawing
  my $xl;
  my $xr;
  foreach $i (0..$#label) {
    $x = $i+1;
    $xl = $x - $width/2;
    $xr = $x + $width/2;
    ($lw, $lq, $med, $hq, $hw) = @{$bwref->{$label[$i]}};
    print GRI "draw box $xl $lq $xr $hq\n";
    print GRI "draw line from $xl $med to $xr $med\n";
    print GRI "draw line from $xl $lw to $xr $lw\n";
    print GRI "draw line from $xl $hw to $xr $hw\n";
    print GRI "draw line from $x $lq to $x $lw\n";
    print GRI "draw line from $x $hq to $x $hw\n";
  }

  close GRI;
  system "convert -density 200 $tempdir/$filename.ps $tempdir/$filename.png";
  !$debug && unlink "$tempdir/$filename.ps";
}

sub do_final_assembly {
  my ($fasta, $k, $cov, $max_cov, $min_contig, $filter_high_coverage) = @_;
  my $name;
  my @f;
  my $Rfile = "commands.R";
  my $high_coverage_switch = "";
  print "Doing the assembly...\n";
  if ($max_cov && $filter_high_coverage) {
    $high_coverage_switch = "-max_coverage $max_cov";
  }
  system "$velveth velvet-temp $k -fasta $fasta";
  system "$velvetg velvet-temp -cov_cutoff $cov -min_contig_lgth $min_contig $high_coverage_switch";
#  if ($max_cov && $filter_high_coverage) {
#    # filter out high coverage
#    open F, ">velvet-temp/contigs-filtered.fa";
#    open G, "velvet-temp/contigs.fa";
#    undef $name;
#    while (<G>) {
#      if (/^>/) {
#        # check coverage
#        @f = split /_/, $_;
#        if ($f[5] > $max_cov) {
#          undef $name;
#        } else {
#          print F;
#          $name = $_;
#        }
#      } else {
#        print F if defined $name;
#      }
#    }
#    close F;
#    close G;
#    return "velvet-temp/contigs-filtered.fa";
    # make some graphs to look at coverage and assembly
    system "grep '^>' velvet-temp/contigs.fa | sed -e 's/_/	/g' > velvet-temp/contigs.fasta-headers";
    open F, ">$Rfile";
    print F <<'__END__';
library(lattice)
lines <- read.table("velvet-temp/contigs.fasta-headers")
pdf(file="ContigCoverage.pdf", width=10, height=10)
hist <- hist(lines$V6, plot=FALSE)
max <- 4*hist$mids[which.max(hist$counts)]
par(mar=c(5,4,4,4))
hist <- hist(lines$V6, breaks = max, xlim = c(0, max), xlab = "Coverage", ylab = "# of contigs", main = "Contig coverage and total assembly length")
assembly <- hist$mids
for (i in 1:length(hist$mids)) {
  assembly[i] <- sum(lines$V4[which(lines$V6 <= hist$breaks[i+1])])
}
par(new=T)
plot(hist$mids, assembly, axes=F, xlab="", ylab="", xlim = c(0,max), type='l', col=4)
axis(side=4, col=4, col.axis=4)
mtext("Total assembly length", side=4, col=4, line=2)
dev.off()

# x-axis is contig length
pdf(file="ContigLength.pdf", width=10, height=10)
hist <- hist(lines$V4, plot=FALSE)
par(mar=c(5,4,4,4))
hist <- hist(lines$V4, breaks = 50, xlab = "Length", ylab = "# of contigs", main = "Contig length and total assembly length")
assembly <- hist$mids
for (i in 1:length(hist$mids)) {
  assembly[i] <- sum(lines$V4[which(lines$V4 <= hist$breaks[i+1])])
}
par(new=T)
plot(hist$mids, assembly, axes=F, xlab="", ylab="", type='l', col=4)
axis(side=4, col=4, col.axis=4)
mtext("Total assembly length", side=4, col=4, line=2)
dev.off()
__END__

    system "R CMD BATCH $Rfile > /dev/null";
    system "convert $tempdir/ContigCoverage.pdf $tempdir/ContigCoverage.png";
    system "convert $tempdir/ContigLength.pdf $tempdir/ContigLength.png";
    return "velvet-temp/contigs.fa";
}
