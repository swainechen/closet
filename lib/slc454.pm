#!/usr/bin/perl
#
# useful procedures for read mapping/snp calling
#
use warnings;
use strict;
use File::Temp;

package slc454;
require Exporter;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw();

sub full_blast_sequences {
  # parse the blast hit if not already parsed
  # return the subject sequences hit and their E-values
  my $type = shift;	# BLASTP, BLASTX, BLASTN, TBLASTN
  my @r = @_;
  my $data;
  if (ref($type) ne "HASH") {	# not yet parsed
    $data = parse_blast_general($type, @r);
  } else {
    $data = $type;
  }
  # data structure is $data->{queryid}->{database}->{HITS}->{hitname}->{hitcounter}->{MAP}->{qcoord}->{scoord}->{QUERY} or last one can be {SUBJECT}
  my $query;
  my $hit;
  my $hitcount;
  my $database;
  my $qcoord;
  my $scoord;
  my $eval;
  my $name;
  my $sequences;
  # return data structure will be $sequences->{query}->{database}->{hitname}->{SEQUENCE} and {EVALUE}
  # and {QUERY}
  foreach $query (keys %$data) {
    foreach $database (keys %{$data->{$query}}) {
      foreach $hit (keys %{$data->{$query}->{$database}->{HITS}}) {
        foreach $hitcount (sort { $a <=> $b } keys %{$data->{$query}->{$database}->{HITS}->{$hit}}) {
          $sequences->{$query}->{$database}->{$hit}->{$hitcount}->{EVALUE} = $data->{$query}->{$database}->{HITS}->{$hit}->{$hitcount}->{EVALUE};
          $sequences->{$query}->{$database}->{$hit}->{$hitcount}->{SEQUENCE} = "";
          $sequences->{$query}->{$database}->{$hit}->{$hitcount}->{QUERY} = "";
          foreach $qcoord (sort { $a <=> $b } keys %{$data->{$query}->{$database}->{HITS}->{$hit}->{$hitcount}->{MAP}}) {
            foreach $scoord (sort { $a <=> $b } keys %{$data->{$query}->{$database}->{HITS}->{$hit}->{$hitcount}->{MAP}->{$qcoord}}) {
              $sequences->{$query}->{$database}->{$hit}->{$hitcount}->{SEQUENCE} .= $data->{$query}->{$database}->{HITS}->{$hit}->{$hitcount}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT};
              $sequences->{$query}->{$database}->{$hit}->{$hitcount}->{QUERY} .= $data->{$query}->{$database}->{HITS}->{$hit}->{$hitcount}->{MAP}->{$qcoord}->{$scoord}->{QUERY};
            }
          }
        }
      }
    }
  }
  return $sequences;
}

sub parse_blast_general {
  #
  # take blastn -m 0 output
  # give a full account of residues that align
  #
  # data structure will be:
  # $data->{queryid}->{database}
  # then these keys:
  #   QUERYLENGTH, DBLENGTH, DBSEQUENCES, HITS
  # under $data->{queryid}->{database}->{HITS}->{hitname} we have integer
  #   hitcounter keys (starting numbering from 1)
  # under each hitcounter we have
  #   HITLENGTH, HITTEXT, EVALUE, SCORE, IDENTITY, GAPS,
  #   ALIGNLENGTH, QSTART, QEND, SSTART, SEND, MAP, and maybe FRAME
  # rest of data in MAP key
  # $data->{queryid}->{database}->{HITS}->{hitname}->{hitcounter}->{MAP}->{qcoord}->{scoord}->{QUERY, SUBJECT}
  #
  my $type = shift;	# BLASTP, BLASTX, BLASTN, TBLASTN
  my @r = @_;
  my $i;
  my ($queryid, $database, $identity, $alignment, @mismatch, $mismatches, @gap, $gaps, $qstart, $qend, $sstart, $send, $eval, $score);
  my ($querylength, $dbsequences, $dblength);
  my $hit;
  my $strand;
  my $data;
  my $hitcounter;
  my $hitlength;
  my @q;
  my @s;
  my ($ql, $qr);
  my ($sl, $sr);
  my @qseq;
  my @sseq;
  my @comp;
  my $comp;
  my $seq;
  my $left;
  my $right;
  my $qcorrect;
  my $scorrect;
  my $j;
  my $qcoord;
  my $scoord;
  my $align_status;
  my ($qincrement, $sincrement, $reverse);
  my $frame;	# only for some of these
  if ($type eq 'BLASTP') {
    $qincrement = 1;
    $sincrement = 1;
    $reverse = 0;
  } elsif ($type eq 'BLASTN') {
    $qincrement = 1;
    $sincrement = 1;
    $reverse = 1;
  } elsif ($type eq 'BLASTX') {
    $qincrement = 3;
    $sincrement = 1;
    $reverse = 0;
  } elsif ($type eq 'TBLASTN') {
    $qincrement = 1;
    $sincrement = 3;
    $reverse = 0;
  } elsif ($type eq 'TBLASTX') {
    $qincrement = 3;
    $sincrement = 3;
    $reverse = 0;
  } else {
    die "Unknown blast type\n";
  }
  undef $queryid;
  undef $database;

  for ($i = 0; $i <= $#r; $i++) {
    chomp $r[$i];
  }

  for ($i = 0; $i <= $#r; $i++) {

    # short circuit if no hits
    if ($r[$i] =~ /\*+ No hits found \*+/) {
      next;
    }

    # get Query, Database
    if ($r[$i] =~ /$type/) {
      undef $queryid;
      undef $database;
      next;
    }
    if ($r[$i] =~ /^Query= (.*)$/) {
      $queryid = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /\(([0-9,]+) letters\)/) {
          $querylength = $1;
          $querylength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $queryid .= " " . $r[$i];
      }
      if (defined $queryid && defined $database) {
        $data->{$queryid}->{$database}->{QUERYLENGTH} = $querylength;
      }
      next;
    }
    if ($r[$i] =~ /^Database: (.*)$/) {
      $database = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /(\d+) sequences; ([0-9,]+) total letters/) {
          $dbsequences = $1;
          $dblength = $2;
          $dblength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $r[$i] =~ s/\s+$//;
        $database .= " " . $r[$i];
      }
      if (defined $queryid && defined $database) {
        $data->{$queryid}->{$database}->{QUERYLENGTH} = $querylength;
        $data->{$queryid}->{$database}->{DBLENGTH} = $dblength;
        $data->{$queryid}->{$database}->{DBSEQUENCES} = $dbsequences;
      }
      next;
    }

    # get hit
    # if there are multiple names for the hit (such as blasting against nr) we just pick up the first one
    if ($r[$i] =~ /^>(.*)/) {
      $hit = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /^\s*Length\s=\s(\d+)/) {
          $hitlength = $1;
          last;
        } else {
          $hit .= $r[$i];
          $hit =~ s/\s+/ /g;
        }
      }
      $hitcounter = 0;
      $identity = 0;
      $alignment = 0;
      $gaps = 0;
      $mismatches = 0;
      $score = 0;
      $eval = -1;
      $qstart = 0;
      $qend = 0;
      $sstart = 0;
      $send = 0;
      @q = ();
      @s = ();
      $ql = 0; $qr = 0;
      @qseq = ();
      $sl = 0; $sr = 0;
      @sseq = ();
      @comp = ();
      $align_status = 0;	# 0 = haven't found one yet
				# 1 = in the middle of parsing
				# 2 = have one and it's parsed fully

      while ($i < $#r) {
        last if $align_status == 2;
        $i++;

        # another hit
        if ($r[$i] =~ /^>/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # another whole record
        if ($r[$i] =~ /^$type/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # all the header stuff
        # sometimes we get Expect(2) = xxx or something
        if ($r[$i] =~ /^\s*Score =\s+(\d+\.?\d*[Ee]?[+-]?\d*) bits.*Expect(?:\(\d+\))? =\s+(.*?)(,\s+Method.*)?$/) {
          $hitcounter++;
          $score = $1;
          $eval = $2;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{HITLENGTH} = $hitlength;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{SCORE} = $score;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{EVALUE} = $eval;
          next;
        }
        if ($r[$i] =~ /^\s*Identities = (\d+)\/(\d+)/) {
          $identity = $1;
          $alignment = $2;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{IDENTITY} = $identity;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{ALIGNLENGTH} = $alignment;
          if ($r[$i] =~ /Gaps = (\d+)\/(\d+)/) {
            $gaps = $1;
            if ($alignment != $2) {
              die "Parse error, query $queryid align length $alignment but parsed $2 in hit $hit\n";
            }
          } else {
            $gaps = 0;
          }
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{GAPS} = $gaps;
          next;
        }
        if ($r[$i] =~ /^\s*Strand = (.*) \/ (.*)/) {
          if ($1 eq $2) {
            $strand = 1;	# whether on the same strand
          } else {
            $strand = 0;
          }
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{STRAND} = $strand;
          next;
        }
        if ($r[$i] =~ /^\s*Frame = ([-+]\d)/) {
          $frame = $1;
          if ($frame < 0) {
            $strand = 0;
          } else {
            $strand = 1;
          }
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{STRAND} = $strand;
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{FRAME} = $frame;
          next;
        }

        # the alignment
        if ($r[$i] =~ /^Query:\s+(\d+)\s*([a-zA-Z*-]+)\s+(\d+)\s*$/) {
          if (!$align_status) {
            $align_status = 1;
          }
          push @{$data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($strand || !$reverse) {
            $ql = $1;
            $qr = $3;
            push @q, $1, $3;
            $seq = $2;
            @qseq = split //, $2;
#            $r[$i] =~ /$seq/;	# this gives an error if one residue on line, and it happens to be Q (matches Query)
#            $left = length $`;
            $r[$i] =~ /^(Query:\s+\d+\s*)/;
            $left = length $1;
            $right = $left + length($seq) - 1;
          } else {
            $ql = $3;
            $qr = $1;
            push @q, $1, $3;
            $seq = $2;
            @qseq = reverse (split //, $2);
            $r[$i] =~ /$seq/;
            $left = length $`;
            $right = $left + length($seq) - 1;
          }
          # next line is the comparison line
          $i++;
          push @{$data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{HITTEXT}}, $r[$i];
          $comp = substr($r[$i], $left, $right - $left + 1);
          while (length $comp < $right - $left + 1) {
            $comp .= " ";
          }
          if ($strand || !$reverse) {
            @comp = split //, $comp;
          } else {
            @comp = reverse (split //, $comp);
          }

          # next line is the subject
          $i++;
          push @{$data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($r[$i] !~ /^Sbjct/) {
            die "Error parsing alignment for $hit\n";
          }
          # sometimes it seems the position runs right into the subject sequence
          # with no space, also sometimes for assembled stuff you might get a
          # stop codon which comes out as '*'
          $r[$i] =~ /^Sbjct:\s+(\d+)\s*([a-zA-Z*-]+)\s+(\d+)\s*$/;
          if ($strand || !$reverse) {
            $sl = $1;
            $sr = $3;
            push @s, $1, $3;
            @sseq = split //, $2;
          } else {
            $sl = $3;
            $sr = $1;
            push @s, $1, $3;
            @sseq = reverse (split //, $2);
          }

          # do the comparisons
          $qcorrect = 0;	# gaps don't count for position
          $scorrect = 0;
          foreach $j (0..$#qseq) {
            # error checking first
            if ($comp[$j] eq '|' && !is_same($qseq[$j], $sseq[$j], $type)) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j]\n";
            }
            if ($comp[$j] eq ' ' && is_same($qseq[$j], $sseq[$j], $type)) {
              die "Error nonidentity parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j]\n";
            }
            if (is_same($comp[$j], $qseq[$j], $type) && !is_same($qseq[$j], $sseq[$j], $type)) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if (is_same($comp[$j], $sseq[$j], $type) && !is_same($qseq[$j], $sseq[$j], $type)) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if ($comp[$j] eq '+' && is_same($qseq[$j], $sseq[$j], $type)) {
              die "Error nonidentity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j]\n";
            }

            if ($qseq[$j] eq '-') {
              $qcorrect += $qincrement;
            }
            if ($sseq[$j] eq '-') {
              $scorrect += $sincrement;
            }
            $qcoord = coordinate($ql, $qr, $j*$qincrement, $qcorrect);
            $scoord = coordinate($sl, $sr, $j*$sincrement, $scorrect);
            $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{QUERY} = $qseq[$j];
            $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT} = $sseq[$j];
          }
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{QSTART} = $q[0];
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{QEND} = $q[$#q];
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{SSTART} = $s[0];
          $data->{$queryid}->{$database}->{HITS}->{$hit}->{$hitcounter}->{SEND} = $s[$#s];
        }
      }

    }
  }

  return $data;

  sub is_same {
    my ($a, $b, $type) = @_;
    $a = uc($a);
    $b = uc($b);
    return 1 if $a eq $b;
    return 0 if $type eq 'BLASTN';	# these special cases are for aa only
    return 1 if $a eq 'B' && ($b eq 'D' || $b eq 'N');
    return 1 if $a eq 'Z' && ($b eq 'E' || $b eq 'Q');
    return 1 if $b eq 'B' && ($a eq 'D' || $a eq 'N');
    return 1 if $b eq 'Z' && ($a eq 'E' || $a eq 'Q');
    return 0;
  }

}

sub full_blastn_alignment_multihit {
  #
  # take blastn -m 0 output
  # give a full account of nts that align
  #
  my @r = @_;
  my $i;
  my ($queryid, $database, $identity, $alignment, @mismatch, $mismatches, @gap, $gaps, $qstart, $qend, $sstart, $send, $eval, $score);
  my ($querylength, $dbsequences, $dblength);
  my $hit;
  my $strand;
  my $data;
  my $hitcounter = 0;
  my @q;
  my @s;
  my ($ql, $qr);
  my ($sl, $sr);
  my @qseq;
  my @sseq;
  my @comp;
  my $comp;
  my $seq;
  my $left;
  my $right;
  my $qcorrect;
  my $scorrect;
  my $j;
  my $qcoord;
  my $scoord;
  my $align_status;
  undef $queryid;
  undef $database;

  for ($i = 0; $i <= $#r; $i++) {

    # short circuit if no hits
    if ($r[$i] =~ /\*+ No hits found \*+/) {
#      return undef;
      undef $queryid;
      next;
    }

    # get Query, Database
#    if (!defined $queryid && $r[$i] =~ /^Query= (.*)$/) {
    if ($r[$i] =~ /^Query= (.*)$/) {
      $queryid = $1;
      $hitcounter = 0;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /\(([0-9,]+) letters\)/) {
          $querylength = $1;
          $querylength =~ s/,//g;
          last;
        }
        if ($r[$i] =~ /^Length=([0-9,]+)/) {
          $querylength = $1;
          $querylength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $queryid .= " " . $r[$i] if length $r[$i];
      }
      $data->{$queryid}->{$hitcounter}->{LENGTH} = $querylength;
      if (defined $database && !defined $data->{$queryid}->{$hitcounter}->{DATABASE}) {
        $data->{$queryid}->{$hitcounter}->{DATABASE} = $database;
      }
      next;
    }
#    if (!defined $database && $r[$i] =~ /^Database: (.*)$/) {
    if ($r[$i] =~ /^Database: (.*)$/) {
      $database = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /(\d+) sequences; ([0-9,]+) total letters/) {
          $dbsequences = $1;
          $dblength = $2;
          $dblength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $database .= $r[$i];
      }
      if (defined $queryid) {
        $data->{$queryid}->{$hitcounter}->{DATABASE} = $database;
      }
      next;
    }

    next if !defined $queryid;

    # get hit
    # if there are multiple names for the hit (such as blasting against nr) we just pick up the first one
    if ($r[$i] =~ /^>(.*)/) {
      $hit = $1;
      $data->{$queryid}->{$hitcounter}->{HIT} = $hit;
      $hitcounter++;
      $identity = 0;
      $alignment = 0;
      $gaps = 0;
      $mismatches = 0;
      $score = 0;
      $eval = -1;
      $qstart = 0;
      $qend = 0;
      $sstart = 0;
      $send = 0;
      @q = ();
      @s = ();
      $ql = 0; $qr = 0;
      @qseq = ();
      $sl = 0; $sr = 0;
      @sseq = ();
      @comp = ();
      $align_status = 0;	# 0 = haven't found one yet
				# 1 = in the middle of parsing
				# 2 = have one and it's parsed fully

      while ($i < $#r) {
        last if $align_status == 2;
        $i++;
        last if $r[$i] =~ /^BLAST/;
        if ($r[$i] =~ /^Query=/) {
          $i--;
          last;
        }

        # another hit
        if ($r[$i] =~ /^>/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # all the header stuff
        if ($r[$i] =~ /^\s*Score\s*=\s*(\d+\.?\d*[Ee]?[+-]?\d*)\s+bits.*Expect\s*=\s+(.*)$/) {
          if ($score) {		# another alignment, same hit
            $data->{$queryid}->{$hitcounter}->{QSTART} = $q[0];
            $data->{$queryid}->{$hitcounter}->{QEND} = $q[$#q];
            $data->{$queryid}->{$hitcounter}->{SSTART} = $s[0];
            $data->{$queryid}->{$hitcounter}->{SEND} = $s[$#s];
            $identity = 0;
            $alignment = 0;
            $gaps = 0;
            $mismatches = 0;
            $score = 0;
            $eval = -1;
            $qstart = 0;
            $qend = 0;
            $sstart = 0;
            $send = 0;
            @q = ();
            @s = ();
            $ql = 0; $qr = 0;
            @qseq = ();
            $sl = 0; $sr = 0;
            @sseq = ();
            @comp = ();
            $align_status = 0;
#            next;
          }
          $score = $1;
          $eval = $2;
          $data->{$queryid}->{$hitcounter}->{SCORE} = $score;
          $data->{$queryid}->{$hitcounter}->{EVALUE} = $eval;
          next;
        }
        if ($r[$i] =~ /^\s*Identities\s*=\s*(\d+)\/(\d+)/) {
          $identity = $1;
          $alignment = $2;
          $data->{$queryid}->{$hitcounter}->{IDENTITY} = $identity;
          $data->{$queryid}->{$hitcounter}->{ALIGNLENGTH} = $alignment;
          if ($r[$i] =~ /Gaps = (\d+)\/(\d+)/) {
            $gaps = $1;
            if ($alignment != $2) {
              die "Parse error, query $queryid align length $alignment but parsed $2 in hit $hit\n";
            }
          } else {
            $gaps = 0;
          }
          $data->{$queryid}->{$hitcounter}->{GAPS} = $gaps;
          next;
        }
        if ($r[$i] =~ /^\s*Strand\s*=\s*(\w+)\s*\/\s*(\w+)/) {
          if ($1 eq $2) {
            $strand = 1;	# whether on the same strand
          } else {
            $strand = 0;
          }
          $data->{$queryid}->{$hitcounter}->{STRAND} = $strand;
          next;
        }

        # the alignment
        if ($r[$i] =~ /^Query\s+(\d+)\s+([a-zA-Z*-]+)\s+(\d+)\s*$/) {
          if (!$align_status) {
            $align_status = 1;
          }
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($strand) {
            $ql = $1;
            $qr = $3;
            push @q, $1, $3;
            $seq = $2;
            @qseq = split //, $2;
            $r[$i] =~ /$seq/;
            $left = length $`;
            $right = $left + length($seq) - 1;
          } else {
            $ql = $3;
            $qr = $1;
            push @q, $1, $3;
            $seq = $2;
            @qseq = reverse (split //, $2);
            $r[$i] =~ /$seq/;
            $left = length $`;
            $right = $left + length($seq) - 1;
          }

          # next line is the comparison line
          $i++;
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          $comp = substr($r[$i], $left, $right - $left + 1);
          while (length $comp < $right - $left + 1) {
            $comp .= " ";
          }
          if ($strand) {
            @comp = split //, $comp;
          } else {
            @comp = reverse (split //, $comp);
          }

          # next line is the subject
          $i++;
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($r[$i] !~ /^Sbjct/) {
            die "Error parsing alignment for $hit\n";
          }
          $r[$i] =~ /^Sbjct\s+(\d+)\s+([a-zA-Z*-]+)\s+(\d+)\s*$/;
          if ($strand) {
            $sl = $1;
            $sr = $3;
            push @s, $1, $3;
            @sseq = split //, $2;
          } else {
            $sl = $3;
            $sr = $1;
            push @s, $1, $3;
            @sseq = reverse (split //, $2);
          }

          # do the comparisons
          $qcorrect = 0;	# gaps don't count for position
          $scorrect = 0;
          foreach $j (0..$#qseq) {
            # error checking first
            if ($comp[$j] eq '|' && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in $queryid, $hit\n";
            }
            if ($comp[$j] eq ' ' && $qseq[$j] eq $sseq[$j]) {
              die "Error nonidentity parsed as identity in $queryid, $hit\n";
            }

#            next if $comp[$j] eq '|';

            if ($qseq[$j] eq '-') {
              $qcorrect++;
            }
            if ($sseq[$j] eq '-') {
              $scorrect++;
            }
            $qcoord = coordinate($ql, $qr, $j, $qcorrect);
            $scoord = coordinate($sl, $sr, $j, $scorrect);
            $data->{$queryid}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{QUERY} = $qseq[$j];
            $data->{$queryid}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT} = $sseq[$j];
          }
        }
      }

      $data->{$queryid}->{$hitcounter}->{QSTART} = $q[0];
      $data->{$queryid}->{$hitcounter}->{QEND} = $q[$#q];
      $data->{$queryid}->{$hitcounter}->{SSTART} = $s[0];
      $data->{$queryid}->{$hitcounter}->{SEND} = $s[$#s];

    }
  }

  return $data;

}

sub full_blastn_alignment {
  #
  # take blastn -m 0 output
  # give a full account of nts that align
  #
  my @r = @_;
  my $i;
  my ($queryid, $database, $identity, $alignment, @mismatch, $mismatches, @gap, $gaps, $qstart, $qend, $sstart, $send, $eval, $score);
  my ($querylength, $dbsequences, $dblength);
  my $hit;
  my $strand;
  my $data;
  my @q;
  my @s;
  my ($ql, $qr);
  my ($sl, $sr);
  my @qseq;
  my @sseq;
  my @comp;
  my $comp;
  my $seq;
  my $left;
  my $right;
  my $qcorrect;
  my $scorrect;
  my $j;
  my $qcoord;
  my $scoord;
  my $align_status;
  undef $queryid;
  undef $database;

  for ($i = 0; $i <= $#r; $i++) {

    # short circuit if no hits
    if ($r[$i] =~ /\*+ No hits found \*+/) {
      return undef;
      next;
    }

    # get Query, Database
    if (!defined $queryid && $r[$i] =~ /^Query= (.*)$/) {
      $queryid = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /\(([0-9,]+) letters\)/) {
          $querylength = $1;
          $querylength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $queryid .= " " . $r[$i];
      }
      $data->{$queryid}->{LENGTH} = $querylength;
      next;
    }
    if (!defined $database && $r[$i] =~ /^Database: (.*)$/) {
      $database = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /(\d+) sequences; ([0-9,]+) total letters/) {
          $dbsequences = $1;
          $dblength = $2;
          $dblength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $database .= $r[$i];
      }
      $data->{$queryid}->{DATABASE} = $database;
      next;
    }

    # get hit
    if ($r[$i] =~ /^>(.*)/) {
      $hit = $1;
      $data->{$queryid}->{HIT} = $hit;
      $identity = 0;
      $alignment = 0;
      $gaps = 0;
      $mismatches = 0;
      $score = 0;
      $eval = -1;
      $qstart = 0;
      $qend = 0;
      $sstart = 0;
      $send = 0;
      @q = ();
      @s = ();
      $ql = 0; $qr = 0;
      @qseq = ();
      $sl = 0; $sr = 0;
      @sseq = ();
      @comp = ();
      $align_status = 0;	# 0 = haven't found one yet
				# 1 = in the middle of parsing
				# 2 = have one and it's parsed fully

      while ($i < $#r) {
        last if $align_status == 2;
        $i++;

        # another hit
        if ($r[$i] =~ /^>/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # all the header stuff
#        if ($r[$i] =~ /^\s*Score =\s+(\d+\.?\d*) bits.*Expect =\s+(.*)$/) {
        if ($r[$i] =~ /^\s*Score =\s+(\d+\.?\d*[Ee]?[+-]?\d*) bits.*Expect =\s+(.*)$/) {
          if ($score) {		# another alignment, same hit, skip this for now
            $align_status = 2 if $align_status;
            next;
          }
          $score = $1;
          $eval = $2;
          $data->{$queryid}->{SCORE} = $score;
          $data->{$queryid}->{EVALUE} = $eval;
          next;
        }
        if ($r[$i] =~ /^\s*Identities = (\d+)\/(\d+)/) {
          $identity = $1;
          $alignment = $2;
          $data->{$queryid}->{IDENTITY} = $identity;
          $data->{$queryid}->{ALIGNLENGTH} = $alignment;
          if ($r[$i] =~ /Gaps = (\d+)\/(\d+)/) {
            $gaps = $1;
            if ($alignment != $2) {
              die "Parse error, query $queryid align length $alignment but parsed $2 in hit $hit\n";
            }
          } else {
            $gaps = 0;
          }
          $data->{$queryid}->{GAPS} = $gaps;
          next;
        }
        if ($r[$i] =~ /^\s*Strand = (.*) \/ (.*)/) {
          if ($1 eq $2) {
            $strand = 1;	# whether on the same strand
          } else {
            $strand = 0;
          }
          $data->{$queryid}->{STRAND} = $strand;
          next;
        }

        # the alignment
        if ($r[$i] =~ /^Query:\s+(\d+)\s+([a-z-]+)\s+(\d+)\s*$/) {
          if (!$align_status) {
            $align_status = 1;
          }
          push @{$data->{$queryid}->{HITTEXT}}, $r[$i];
          if ($strand) {
            $ql = $1;
            $qr = $3;
            push @q, $1, $3;
            $seq = $2;
            @qseq = split //, $2;
            $r[$i] =~ /$seq/;
            $left = length $`;
            $right = $left + length($seq) - 1;
          } else {
            $ql = $3;
            $qr = $1;
            push @q, $1, $3;
            $seq = $2;
            @qseq = reverse (split //, $2);
            $r[$i] =~ /$seq/;
            $left = length $`;
            $right = $left + length($seq) - 1;
          }

          # next line is the comparison line
          $i++;
          push @{$data->{$queryid}->{HITTEXT}}, $r[$i];
          $comp = substr($r[$i], $left, $right - $left + 1);
          while (length $comp < $right - $left + 1) {
            $comp .= " ";
          }
          if ($strand) {
            @comp = split //, $comp;
          } else {
            @comp = reverse (split //, $comp);
          }

          # next line is the subject
          $i++;
          push @{$data->{$queryid}->{HITTEXT}}, $r[$i];
          if ($r[$i] !~ /^Sbjct/) {
            die "Error parsing alignment for $hit\n";
          }
          $r[$i] =~ /^Sbjct:\s+(\d+)\s+([a-z-]+)\s+(\d+)\s*$/;
          if ($strand) {
            $sl = $1;
            $sr = $3;
            push @s, $1, $3;
            @sseq = split //, $2;
          } else {
            $sl = $3;
            $sr = $1;
            push @s, $1, $3;
            @sseq = reverse (split //, $2);
          }

          # do the comparisons
          $qcorrect = 0;	# gaps don't count for position
          $scorrect = 0;
          foreach $j (0..$#qseq) {
            # error checking first
            if ($comp[$j] eq '|' && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in $queryid, $hit\n";
            }
            if ($comp[$j] eq ' ' && $qseq[$j] eq $sseq[$j]) {
              die "Error nonidentity parsed as identity in $queryid, $hit\n";
            }

#            next if $comp[$j] eq '|';

            if ($qseq[$j] eq '-') {
              $qcorrect++;
            }
            if ($sseq[$j] eq '-') {
              $scorrect++;
            }
            $qcoord = coordinate($ql, $qr, $j, $qcorrect);
            $scoord = coordinate($sl, $sr, $j, $scorrect);
            $data->{$queryid}->{MAP}->{$qcoord}->{$scoord}->{QUERY} = $qseq[$j];
            $data->{$queryid}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT} = $sseq[$j];
          }
        }
      }

      $data->{$queryid}->{QSTART} = $q[0];
      $data->{$queryid}->{QEND} = $q[$#q];
      $data->{$queryid}->{SSTART} = $s[0];
      $data->{$queryid}->{SEND} = $s[$#s];

    }
  }

  return $data;

}


sub parse_blastn {
  my @r = @_;
  my $i;
  my $data;
  my ($queryid, $database, $identity, $alignment, @mismatch, $mismatches, @gap, $gaps, $qstart, $qend, $sstart, $send, $eval, $score);
  my ($querylength, $dbsequences, $dblength);
  my $hit;
  my $strand;
  my $hitcounter = 0;
  my @q;
  my @s;
  my ($ql, $qr);
  my ($sl, $sr);
  my @qseq;
  my @sseq;
  my @comp;
  my $comp;
  my $seq;
  my $left;
  my $right;
  my $qcorrect;
  my $scorrect;
  my $j;
  my $qcoord;
  my $scoord;
  my $align_status;
  undef $queryid;
  undef $database;

  for ($i = 0; $i <= $#r; $i++) {

    # short circuit if no hits
    if ($r[$i] =~ /\*+ No hits found \*+/) {
      return undef;
    }

    # get Query, Database
    if (!defined $queryid && $r[$i] =~ /^Query= (.*)$/) {
      $queryid = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /\((\d+) letters\)/) {
          $querylength = $1;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $queryid .= $r[$i];
      }
      $data->{$queryid}->{LENGTH} = $querylength;
      next;
    }
    if (!defined $database && $r[$i] =~ /^Database: (.*)$/) {
      $database = $1;
      while ($i < $#r) {
        $i++;
        if ($r[$i] =~ /(\d+) sequences; ([0-9,]+) total letters/) {
          $dbsequences = $1;
          $dblength = $2;
          $dblength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $database .= $r[$i];
      }
      $data->{$queryid}->{DATABASE} = $database;
      next;
    }

    # get hit
    if ($r[$i] =~ /^>(.*)/) {
      $hit = $1;
      $data->{$queryid}->{HIT} = $hit;
      $identity = 0;
      $alignment = 0;
      $gaps = 0;
      $mismatches = 0;
      $score = 0;
      $eval = -1;
      $qstart = 0;
      $qend = 0;
      $sstart = 0;
      $send = 0;
      @q = ();
      @s = ();
      $ql = 0; $qr = 0;
      @qseq = ();
      $sl = 0; $sr = 0;
      @sseq = ();
      @comp = ();
      $align_status = 0;	# 0 = haven't found one yet
				# 1 = in the middle of parsing
				# 2 = have one and it's parsed fully

      while ($i < $#r) {
        last if $align_status == 2;
        $i++;

        # another hit
        if ($r[$i] =~ /^>/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # all the header stuff
#        if ($r[$i] =~ /^\s*Score =\s+(\d+\.?\d*) bits.*Expect =\s+(.*)$/) {
        if ($r[$i] =~ /^\s*Score\s*=\s+(\d+\.?\d*[Ee]?[+-]?\d*) bits.*Expect =\s+(.*)$/) {
          if ($score) {		# another alignment, same hit, skip this for now
            $align_status = 2 if $align_status;
            next;
          }
          $score = $1;
          $eval = $2;
          $data->{$queryid}->{SCORE} = $score;
          $data->{$queryid}->{EVALUE} = $eval;
          next;
        }
        if ($r[$i] =~ /^\s*Identities = (\d+)\/(\d+)/) {
          $identity = $1;
          $alignment = $2;
          $data->{$queryid}->{IDENTITY} = $identity;
          $data->{$queryid}->{ALIGNLENGTH} = $alignment;
          if ($r[$i] =~ /Gaps = (\d+)\/(\d+)/) {
            $gaps = $1;
            if ($alignment != $2) {
              die "Parse error, query $queryid align length $alignment but parsed $2 in hit $hit\n";
            }
          } else {
            $gaps = 0;
          }
          $data->{$queryid}->{GAPS} = $gaps;
          next;
        }
        if ($r[$i] =~ /^\s*Strand = (.*) \/ (.*)/) {
          if ($1 eq $2) {
            $strand = 1;	# whether on the same strand
          } else {
            $strand = 0;
          }
          $data->{$queryid}->{STRAND} = $strand;
          next;
        }

        # the alignment
        if ($r[$i] =~ /^Query:\s+(\d+)\s+([a-z-]+)\s+(\d+)\s*$/) {
          if (!$align_status) {
            $align_status = 1;
            $data->{$queryid}->{SNP} = {};
            $data->{$queryid}->{INS} = {};
            $data->{$queryid}->{DEL} = {};
          }
          $ql = $1;
          $qr = $3;
          push @q, $1, $3;
          $seq = $2;
          @qseq = split //, $2;
          $r[$i] =~ /$seq/;
          $left = length $`;
          $right = $left + length($seq) - 1;

          # next line is the comparison line
          $i++;
          $comp = substr($r[$i], $left, $right - $left + 1);
          while (length $comp < $right - $left + 1) {
            $comp .= " ";
          }
          @comp = split //, $comp;

          # next line is the subject
          $i++;
          if ($r[$i] !~ /^Sbjct/) {
            die "Error parsing alignment for $hit\n";
          }
          $r[$i] =~ /^Sbjct:\s+(\d+)\s+([a-z-]+)\s+(\d+)\s*$/;
          $sl = $1;
          $sr = $3;
          push @s, $1, $3;
          @sseq = split //, $2;

          # do the comparisons
          $qcorrect = 0;	# gaps don't count for position
          $scorrect = 0;
          foreach $j (0..$#qseq) {
            # error checking first
            if ($comp[$j] eq '|' && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in $queryid, $hit\n";
            }
            if ($comp[$j] eq ' ' && $qseq[$j] eq $sseq[$j]) {
              die "Error nonidentity parsed as identity in $queryid, $hit\n";
            }
            next if $comp[$j] eq '|';
            if ($qseq[$j] eq '-') {
              $qcorrect++;
              $qcoord = coordinate($ql, $qr, $j, $qcorrect);
              $data->{$queryid}->{DEL}->{$qcoord}->{COORD} = coordinate($sl, $sr, $j, $scorrect);
              $data->{$queryid}->{DEL}->{$qcoord}->{NUCLEOTIDE} = $sseq[$j];
              next;
            }
            if ($sseq[$j] eq '-') {
              $scorrect++;
              $qcoord = coordinate($ql, $qr, $j, $qcorrect);
              $data->{$queryid}->{INS}->{$qcoord}->{COORD} = coordinate($sl, $sr, $j, $scorrect);
              $data->{$queryid}->{INS}->{$qcoord}->{NUCLEOTIDE} = $qseq[$j];
              next;
            }
            $qcoord = coordinate($ql, $qr, $j, $qcorrect);
            $data->{$queryid}->{SNP}->{$qcoord}->{COORD} = coordinate($sl, $sr, $j, $scorrect);
            $data->{$queryid}->{SNP}->{$qcoord}->{NUCLEOTIDE} = $sseq[$j];
            $data->{$queryid}->{SNP}->{$qcoord}->{MUTATION} = $qseq[$j];
          }
        }
      }

      $data->{$queryid}->{QSTART} = $q[0];
      $data->{$queryid}->{QEND} = $q[$#q];
      $data->{$queryid}->{SSTART} = $s[0];
      $data->{$queryid}->{SEND} = $s[$#s];

    }

  }

  return $data;

}

sub coordinate {
  my ($l, $r, $i, $c) = @_;
  # take left and right coordinate, position in string (0-based), correction
  # (i.e. # of gaps)
  # return corresponding coordinate of that position
  if ($l < $r) {
    return ($l+$i-$c);
  } else {
    return ($l-$i+$c);
  }
}

sub parse_repeat {
  # this takes name of a data file generated by icasss
  # returns a hash ref
  my ($file) = @_;
  my $data;
  my $score;
  my $match;
  my ($s1, $e1, $s2, $e2);
  my $t;

  open F, $file;
  while (<F>) {
    if (/SCORE = (\d+) MATCH = (.*) \% .* PARENT Length = (\d+)/) {
      if ($1 < 0.9*$3) {        # don't want whole genome match
        $score = $1;
        $match = $2;
      }
    }
    if (/Match region is (\d+) to (\d+) on .* and (\d+) to (\d+) on /) {
      ($s1, $e1, $s2, $e2) = ($1, $2, $3, $4);
      if ($s1 > $e1) {
        $t = $s1;
        $s1 = $e1;
        $e1 = $t;
      }
      if ($s2 > $e2) {
        $t = $s2;
        $s2 = $e2;
        $e2 = $t;
      }
      next if !$score || !$match;
      $data->{$s1}->{$e1}->{SCORE} = $score;
      $data->{$s1}->{$e1}->{MATCH} = $match;
      $data->{$s2}->{$e2}->{SCORE} = $score;
      $data->{$s2}->{$e2}->{MATCH} = $match;
      undef $score;
      undef $match;
      undef $s1;
      undef $e1;
      undef $s2;
      undef $e2;
    }
  }
  close F;
  return $data;
}

sub parse_injection {
  # this takes name of a data file generated by injection.pl
  # returns a hash ref
  my ($file) = @_;
  my $data;
  my @f;
  my $i;
  my @g;
  open F, $file;
  while (<F>) {
    chomp;
    next if /^$/ || /^#/;
    @f = split /\t/, $_; 
    if ($#f > 1) {
      foreach $i (1..$#f) {
        @g = split /\.\./, $f[$i];
        if ($g[0] <= $g[1]) {
          $data->{$g[0]}->{$g[1]} = 1;
        } else {
          $data->{$g[1]}->{$g[0]} = 1;
        }
      }
    }
  }
  close F;
  return $data;
}

sub inrepeat {
  # this takes a hash ref from parse_repeat
  # and a position, returns true/false if position is within a repeat
  my ($data, $nt) = @_;
  my $i;
  my $j;
  foreach $i (keys %$data) {
    next if $nt < $i;
    foreach $j (keys %{$data->{$i}}) {
      return (1) if $nt <= $j;
    }
  }
  return 0;
}

sub print_known_readsets {
  my ($dbh) = @_;
  my $sth = $dbh->prepare("SELECT ReadSetID, SetName FROM ReadSets");
  $sth->execute;
  my @data;
  print "Known ReadSetIDs:\n";
  while (@data = $sth->fetchrow_array) {
    print join ("\t", @data), "\n";
  }
}

sub print_known_mappings {
  my ($dbh) = @_;
  my $sth = $dbh->prepare("SHOW TABLES");
  $sth->execute;
  my @data;
  while (@data = $sth->fetchrow_array) {
    print "$data[0]\n" if $data[0] =~ /^Mapping_/;
  }
}

sub poisson_cumulative {
  my ($total_sequence, $chromosome_length, $prob) = @_;

  sub fact {
    my ($n) = @_;
    my $i;
    my $f = 1;
    foreach $i (2..$n) {
      $f *= $i;
    }
    return $f;
  }

  # calculate Poisson parameters
  my $k;	# k, number of occurrences
  my $l;	# lambda
  my $p;	# cumulative probability
  $l = $total_sequence/$chromosome_length;
  $p = 0;
  $k = -1;
  while ($p < $prob) {
    $k++;
    $p += exp(-$l) * $l**$k / fact($k);
  }
  return $k;
}

sub full_blast_alignment_multihit {
  #
  # take any blast output
  # give a full account of nts that align
  # data structure
  # $data is a ref to a hash with keys that are all the queryids
  # $data->{$queryid} has keys that are numbers that count the hits, 0-based
  # $data->{$queryid}->{$hitcounter} has keys:
  #   LENGTH, DATABASE, HIT, HITTEXT
  #   SCORE, EVAL, IDENTITY, ALIGNLENGTH, GAPS
  #   STRAND, maybe FRAME if relevant
  #   QSTART, QEND, SSTART, SEND
  #   and MAP
  # $data->{$queryid}->{$hitcounter}->{MAP} has query coord, then subject coord
  # then QUERY and SUBJECT, which are both sequences
  # i.e. $data->{$queryid}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT}
  #
  my $type = shift;	# BLASTP, BLASTX, BLASTN, TBLASTN, TBLASTX
  my @r = @_;
  my $i;
  my ($queryid, $database, $identity, $alignment, @mismatch, $mismatches, @gap, $gaps, $qstart, $qend, $sstart, $send, $eval, $score);
  my ($querylength, $dbsequences, $dblength);
  my $hit;
  my $strand;
  my $data;
  my $hitcounter = 0;
  my @q;
  my @s;
  my ($ql, $qr);
  my ($sl, $sr);
  my @qseq;
  my @sseq;
  my @comp;
  my $comp;
  my $seq;
  my $left;
  my $right;
  my $qcorrect;
  my $scorrect;
  my $j;
  my $qcoord;
  my $scoord;
  my $align_status;
  my ($qincrement, $sincrement, $reverse);
  my $frame;
  undef $queryid;
  undef $database;

  if ($type eq 'BLASTP') {
    $qincrement = 1;
    $sincrement = 1;
    $reverse = 0;
  } elsif ($type eq 'BLASTN') {
    $qincrement = 1;
    $sincrement = 1;
    $reverse = 1;
  } elsif ($type eq 'BLASTX') {
    $qincrement = 3;
    $sincrement = 1;
    $reverse = 0;
  } elsif ($type eq 'TBLASTN') {
    $qincrement = 1;
    $sincrement = 3;
    $reverse = 0;
  } elsif ($type eq 'TBLASTX') {
    $qincrement = 3;
    $sincrement = 3;
    $reverse = 0;
  } else {
    die "Unknown blast type\n";
  }

  for ($i = 0; $i <= $#r; $i++) {
    chomp $r[$i];

    # short circuit if no hits
    if ($r[$i] =~ /\*+ No hits found \*+/) {
      undef $queryid;
      next;
    }

    # get Query, Database
    if ($r[$i] =~ /^Query= (.*)$/) {
      $queryid = $1;
      $hitcounter = 0;
      while ($i < $#r) {
        $i++;
        chomp $r[$i];
        if ($r[$i] =~ /\(([0-9,]+) letters\)/) {
          $querylength = $1;
          $querylength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $queryid .= " " . $r[$i];
      }
      $data->{$queryid}->{$hitcounter}->{LENGTH} = $querylength;
      next;
    }
    if ($r[$i] =~ /^Database: (.*)$/) {
      $database = $1;
      while ($i < $#r) {
        $i++;
        chomp $r[$i];
        if ($r[$i] =~ /(\d+) sequences; ([0-9,]+) total letters/) {
          $dbsequences = $1;
          $dblength = $2;
          $dblength =~ s/,//g;
          last;
        }
        $r[$i] =~ s/^\s+//;
        $r[$i] =~ s/\s+$//;
        $database .= " " . $r[$i];
      }
      $data->{$queryid}->{$hitcounter}->{DATABASE} = $database;
      next;
    }

    next if !defined $queryid;

    # get hit
    # if there are multiple names for the hit (such as blasting against nr) we just pick up the first one
    if ($r[$i] =~ /^>(.*)/) {
      $hit = $1;
      $data->{$queryid}->{$hitcounter}->{HIT} = $hit;
      $identity = 0;
      $alignment = 0;
      $gaps = 0;
      $mismatches = 0;
      $score = 0;
      $eval = -1;
      $qstart = 0;
      $qend = 0;
      $sstart = 0;
      $send = 0;
      @q = ();
      @s = ();
      $ql = 0; $qr = 0;
      @qseq = ();
      $sl = 0; $sr = 0;
      @sseq = ();
      @comp = ();
      $align_status = 0;	# 0 = haven't found one yet
				# 1 = in the middle of parsing
				# 2 = have one and it's parsed fully

      while ($i < $#r) {
        last if $align_status == 2;
        $i++;
        last if $r[$i] =~ /^BLAST/;

        # another hit
        if ($r[$i] =~ /^>/) {
          $align_status = 2 if $align_status;
          $i--;
          last;
        }

        # all the header stuff
        if ($r[$i] =~ /^\s*Score =\s+(\d+\.?\d*[Ee]?[+-]?\d*) bits.*Expect =\s+(.*)$/) {
          if ($score) {		# another alignment, same hit
            $data->{$queryid}->{$hitcounter}->{QSTART} = $q[0];
            $data->{$queryid}->{$hitcounter}->{QEND} = $q[$#q];
            $data->{$queryid}->{$hitcounter}->{SSTART} = $s[0];
            $data->{$queryid}->{$hitcounter}->{SEND} = $s[$#s];
            $hitcounter++;
            $identity = 0;
            $alignment = 0;
            $gaps = 0;
            $mismatches = 0;
            $score = 0;
            $eval = -1;
            $qstart = 0;
            $qend = 0;
            $sstart = 0;
            $send = 0;
            @q = ();
            @s = ();
            $ql = 0; $qr = 0;
            @qseq = ();
            $sl = 0; $sr = 0;
            @sseq = ();
            @comp = ();
            $align_status = 0;
          }
          $score = $1;
          $eval = $2;
          $data->{$queryid}->{$hitcounter}->{SCORE} = $score;
          $data->{$queryid}->{$hitcounter}->{EVALUE} = $eval;
          next;
        }
        if ($r[$i] =~ /^\s*Identities = (\d+)\/(\d+)/) {
          $identity = $1;
          $alignment = $2;
          $data->{$queryid}->{$hitcounter}->{IDENTITY} = $identity;
          $data->{$queryid}->{$hitcounter}->{ALIGNLENGTH} = $alignment;
          if ($r[$i] =~ /Gaps = (\d+)\/(\d+)/) {
            $gaps = $1;
            if ($alignment != $2) {
              die "Parse error, query $queryid align length $alignment but parsed $2 in hit $hit\n";
            }
          } else {
            $gaps = 0;
          }
          $data->{$queryid}->{$hitcounter}->{GAPS} = $gaps;
          next;
        }
        if ($r[$i] =~ /^\s*Strand =\s+(\w+)\s+\/\s+(\w+)/) {
          if ($1 eq $2) {
            $strand = 1;	# whether on the same strand
          } else {
            $strand = 0;
          }
          $data->{$queryid}->{$hitcounter}->{STRAND} = $strand;
          next;
        }
        if ($r[$i] =~ /^\s*Frame = ([-+]\d)/) {
          $frame = $1;
          if ($frame < 0) {
            $strand = 0;
          } else {
            $strand = 1;
          }
          $data->{$queryid}->{$hitcounter}->{STRAND} = $strand;
          $data->{$queryid}->{$hitcounter}->{FRAME} = $frame;
          next;
        }

        # the alignment
        if ($r[$i] =~ /^Query:\s+(\d+)\s*([a-zA-Z-*]+)\s+(\d+)\s*$/) {
          if (!$align_status) {
            $align_status = 1;
          }
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($strand) {
            $ql = $1;
            $qr = $3;
            push @q, $1, $3;
            $seq = $2;
            @qseq = split //, $2;
            $r[$i] =~ /^(Query:\s+\d+\s*)/;
            $left = length $1;
            $right = $left + length($seq) - 1;
          } else {
            $ql = $3;
            $qr = $1;
            push @q, $1, $3;
            $seq = $2;
            @qseq = reverse (split //, $2);
            $r[$i] =~ /^(Query:\s+\d+\s*)/;
            $left = length $1;
            $right = $left + length($seq) - 1;
          }

          # next line is the comparison line
          $i++;
          chomp $r[$i];
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          $comp = substr($r[$i], $left, $right - $left + 1);
          while (length $comp < $right - $left + 1) {
            $comp .= " ";
          }
          if ($strand) {
            @comp = split //, $comp;
          } else {
            @comp = reverse (split //, $comp);
          }

          # next line is the subject
          $i++;
          chomp $r[$i];
          push @{$data->{$queryid}->{$hitcounter}->{HITTEXT}}, $r[$i];
          if ($r[$i] !~ /^Sbjct/) {
            die "Error parsing alignment for $hit\n";
          }
          $r[$i] =~ /^Sbjct:\s+(\d+)\s*([a-zA-Z-*]+)\s+(\d+)\s*$/;
          if ($strand) {
            $sl = $1;
            $sr = $3;
            push @s, $1, $3;
            @sseq = split //, $2;
          } else {
            $sl = $3;
            $sr = $1;
            push @s, $1, $3;
            @sseq = reverse (split //, $2);
          }

          # do the comparisons
          $qcorrect = 0;	# gaps don't count for position
          $scorrect = 0;
          foreach $j (0..$#qseq) {
            # error checking first
            if ($comp[$j] eq '|' && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if ($comp[$j] eq ' ' && $qseq[$j] eq $sseq[$j]) {
              die "Error nonidentity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if ($comp[$j] eq $qseq[$j] && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if ($comp[$j] eq $sseq[$j] && $qseq[$j] ne $sseq[$j]) {
              die "Error identity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }
            if ($comp[$j] eq '+' && $qseq[$j] eq $sseq[$j]) {
              die "Error nonidentity not parsed as identity in Query $queryid, Hit $hit, Column $j, Query $qseq[$j], Subject $sseq[$j], Comp $comp[$j], line $r[$i]\n";
            }

            if ($qseq[$j] eq '-') {
              $qcorrect += $qincrement;
            }
            if ($sseq[$j] eq '-') {
              $scorrect += $sincrement;
            }
            $qcoord = coordinate($ql, $qr, $j*$qincrement, $qcorrect);
            $scoord = coordinate($sl, $sr, $j*$sincrement, $scorrect);
            $data->{$queryid}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{QUERY} = $qseq[$j];
            $data->{$queryid}->{$hitcounter}->{MAP}->{$qcoord}->{$scoord}->{SUBJECT} = $sseq[$j];
          }
        }
      }

      $data->{$queryid}->{$hitcounter}->{QSTART} = $q[0];
      $data->{$queryid}->{$hitcounter}->{QEND} = $q[$#q];
      $data->{$queryid}->{$hitcounter}->{SSTART} = $s[0];
      $data->{$queryid}->{$hitcounter}->{SEND} = $s[$#s];
    }
  }
  return $data;
}

# End with 1 to make Perl happy
1;
