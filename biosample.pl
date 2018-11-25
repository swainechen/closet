#!/usr/bin/perl -w
#
# parse biosample xml
#
use XML::Simple;
use Data::Dumper;
my @record = ();
my $parse = {};
while (<>) {
  if (/^<\?xml version="1.0"/) {
    # new record
    if (scalar @record) {
      parse_biosample($parse, @record);
    }
    @record = ($_);
  } else {
    push @record, $_;
  }
}
if (scalar @record) {
  parse_biosample($parse, @record);
}
my $fields = ();
foreach $key (keys %$parse) {
  foreach $i (keys %{$parse->{$key}}) {
    $fields->{$i} = 1;
  }
}
print join ("\t", "# SRA Sample", sort keys %$fields), "\n";
foreach $key (keys %$parse) {
  @out = ($key);
  foreach $i (sort keys %$fields) {
    if (!defined $parse->{$key}->{$i}) {
      push @out, "";
    } else {
      push @out, $parse->{$key}->{$i};
    }
  }
  print join ("\t", @out), "\n";
}

sub parse_biosample {
  my ($r, @x) = @_;
  my $d = XMLin(join ("\n", @x));
  my $t;
  my $hold;
  my $key = "";
  # this really should be a biosample at the top
  if (defined $d->{BioSample}) {
    $t = $d->{BioSample};
    if (defined $t->{Ids}) {
      # from here should have biosample ID and SRA [EDS]RSxxx ID
      if (defined $t->{Ids}->{Id} && ref($t->{Ids}->{Id}) eq "ARRAY") {
        $t = $t->{Ids}->{Id};
        foreach $i (0..$#$t) {
          if (defined $t->[$i]->{db}) {
            if ($t->[$i]->{db} eq "BioSample" && defined $t->[$i]->{content}) {
              $hold->{BioSample} = $t->[$i]->{content};
            } elsif ($t->[$i]->{db} eq "SRA" && defined $t->[$i]->{content}) {
              $hold->{SRA} = $t->[$i]->{content};
              $key = $t->[$i]->{content};
            } elsif ( defined $t->[$i]->{db_label} &&
                      $t->[$i]->{db_label} eq "Sample name" &&
                      defined $t->[$i]->{content} ) {
              $hold->{Ids_SampleName} = $t->[$i]->{content};
            }
          }
        }
        $t = $d->{BioSample};
      }
    }
    if (defined $t->{Description}) {
      if (defined $t->{Description}->{Title}) {
        $hold->{Description_Title} = $t->{Description}->{Title};
      }
    }
    if (defined $t->{Attributes}) {
      if (defined $t->{Attributes}->{Attribute} &&
          ref($t->{Attributes}->{Attribute}) eq "ARRAY") {
        $t = $t->{Attributes}->{Attribute};
        foreach $i (0..$#$t) {
          if (defined $t->[$i]->{harmonized_name}) {
            $hold->{$t->[$i]->{harmonized_name}} = $t->[$i]->{content}
          }
        }
        $t = $d->{BioSample};
      }
    }
    if ($key ne "") {
      $r->{$key} = $hold;
    }
  }
#  foreach $i (sort keys %$r) {
#    print join ("\t", $i, $r->{$i}), "\n";
#  }
#  print Dumper($d);
}
