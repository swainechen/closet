#!/usr/bin/perl -w
#
# monitor machines and launch when there is capacity
#
use File::Temp;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

###############################
# These need to be configured #
###############################
my @zonelist = ("ap-southeast-1a", "ap-southeast-1b", "ap-southeast-1c");
my $TAG = "chenlab_auto";
my $AMI = "ami-xxxxxxxxxxxxxxxxx";
my $KEY = "key_name_here";
my $PROFILE = "aws_configure_profile";
my $query_wait = 180;	# in seconds
my $submit_wait = 0;	#in seconds
my $MAX_SIMULTANEOUS = 100;

# set up limits - hash ref keyed on Instance Type then Zone
my $limit;
$limit->{__TOTAL__} = 4510;
foreach $i (@zonelist) {
  $limit->{"r5a.large"}->{$i} = 500;
  $limit->{"r5.large"}->{$i} = 500;
  $limit->{"r4.large"}->{$i} = 500;
}

$USERDATA_TEMPLATE='#!/bin/bash
# some initial info, tagging
RUN=__RUN__
TAG=__TAG__
ID=$(GET http://169.254.169.254/latest/meta-data/instance-id)
su - ubuntu -c "aws ec2 create-tags --resources $ID --tags \"Key=Name,Value=Auto-$RUN\""
su - ubuntu -c "aws ec2 create-tags --resources $ID --tags \"Key=Project,Value=$TAG\""

# the main commands - note these are all run as root, and the working directory is / (not /root)
# wget https://some_reference_sequence
# bwa index some_reference_sequence
# aws s3 cp s3://some_bucket/${RUN}_R1.fastq.gz /home/ubuntu
# aws s3 cp s3://some_bucket/${RUN}_R1.fastq.gz /home/ubuntu
# bwa mem some_reference_sequence ${RUN}_R1.fastq.gz ${RUN}_R2.fastq.gz > ${RUN}.bam

# make sure to copy things back to some S3 bucket or persistent storage
# aws s3 cp ${RUN}.bam s3://some_bucket/output/

# shut down the instance
init 0
';

# this first command is just for on demand
$CLI_EC2_TEMPLATE="aws --profile __PROFILE__ ec2 run-instances --image-id __AMI__ --count 1 --instance-type __TYPE__ --key-name __KEY__ --user-data __USERDATA__ --instance-initiated-shutdown-behavior terminate --placement AvailabilityZone=__ZONE__";
# this version is for spot instances
$CLI_EC2_TEMPLATE="aws --profile __PROFILE__ ec2 run-instances --image-id __AMI__ --count 1 --instance-type __TYPE__ --count 1 --key-name __KEY__ --user-data __USERDATA__ --instance-initiated-shutdown-behavior terminate --placement AvailabilityZone=__ZONE__ --iam-instance-profile Name=\"IAM-slchen-lab-spot\" --instance-market-options \"MarketType=spot,SpotOptions={SpotInstanceType=one-time}\"";

################################
# End of typical configuration #
################################

my $tempdir = File::Temp::tempdir( CLEANUP => 1 );
my $TMPFILE = "$tempdir/temp.userdata";
my $STATUS = "/usr/local/bin/get_ec2_status.pl";
my $i;
my @input;
my $summary;
my $available_total;
my $available;
my $start = 0;
my $print_help = 0;
my $type;
my $zone;
my $this_count;
my $hit_limit_instance;
my $hit_limit_type;
my $t1;
my $t2;

GetOptions (
  'start=i' => \$start,
  'profile=s' => \$PROFILE,
  'ami=s' => \$AMI,
  'wait=i' => \$query_wait,
  'tag=s' => \$TAG,
  'help' => \$print_help
);

if ($print_help) {
  print "Usage: $0 <RunIDs> [ -start <int> ] [ -profile <AWS profile> ] [ -ami <AMI ID> ] [ -wait <seconds> ] [ -tag <custom tag> ]\n";
  print "       $0 -help\n";
  print "  start is a 0-based index for which RunID to start with, made for rerunning if there's an error\n";
  print "  profile should be something in the ~/.aws/credentials file\n";
  print "  ami should be like ami-xxxxxxxxxx\n";
  print "  tag is any string, probably need to escape it if has special characters\n";
  print "  help prints this help\n";
  print "Will submit jobs, will quit as soon as there is any error\n";
  print "Current defaults:\n";
  print "  profile = $PROFILE\n";
  print "  ami = $AMI\n";
  print "  tag = $TAG\n";
  print "  wait = $query_wait\n";
  exit;
}

@input = ();
while (<>) {
  next if /^#/;
  next if /^$/;
  chomp;
  push @input, $_;
}

print "Found ", scalar(@input) - $start, " jobs to submit with AMI $AMI, TAG $TAG using types:\n";
foreach $i (sort keys %$limit) {
  next if $i =~ /^__/;
  print "  $i\n";
}
print "Will submit in the following zones:\n";
foreach $i (@zonelist) {
  print "  $i\n";
}
print "Starting from $start (0-based) - note Slots below are 1-based counting:\n";

# main watching and submitting loop
$i = $start;
$this_count = 0;
while ($i < scalar(@input)) {
  $hit_limit_instance = 0;
  $hit_limit_type = 0;
  $summary = get_status();
  while ($summary->{__TOTAL__} >= $limit->{__TOTAL__}) {
    system("date");
    print "Waiting $query_wait seconds for available instance quota...\n";
    sleep $query_wait;
    $summary = get_status();
  }
  $available_total = $limit->{__TOTAL__} - $summary->{__TOTAL__};
  # should be hash ref keyed on State, then InstanceType, then Zone
  $this_count = 0;
  foreach $type (keys %$limit) {
    next if $type =~ /^__/;
    foreach $zone (keys %{$limit->{$type}}) {
      if (!defined $summary->{"running"}->{$type}->{$zone}) {
        $available = $limit->{$type}->{$zone}
      } else {
        $available = $limit->{$type}->{$zone} - $summary->{"running"}->{$type}->{$zone}
      }
      foreach $j (1..$available) {
        last if !$available_total;
        unlink $TMPFILE;
        open F, ">", $TMPFILE;
        $userdata_out = $USERDATA_TEMPLATE;
        $userdata_out =~ s/__RUN__/$input[$i]/g;
        $userdata_out =~ s/__TAG__/$TAG/g;
        print F $userdata_out;
        close F;
        $launch_command = $CLI_EC2_TEMPLATE;
        $launch_command =~ s/__PROFILE__/$PROFILE/g;
        $launch_command =~ s/__AMI__/$AMI/g;
        $launch_command =~ s/__TYPE__/$type/g;
        $launch_command =~ s/__KEY__/$KEY/g;
        $launch_command =~ s|__USERDATA__|file://$TMPFILE|g;
        $launch_command =~ s/__ZONE__/$zone/g;
        print "Slot ", $i+1, ", Run $input[$i] (using $j/$available for $type in $zone):";
        $launch_result = `$launch_command 2>&1`;
        if ($? == 0) {
          print " OK!\n";
          $i++;
          $this_count++;
          $available_total--;
          exit if $i >= scalar(@input);
          if ($this_count >= $MAX_SIMULTANEOUS) {
            print "Mark at $this_count submissions:\n";
            system("date");
            sleep $submit_wait;
            $this_count = 0;
          }
        } else {
          if ($launch_result =~ /An error occurred \(InstanceLimitExceeded\)/) {
            # we can recover if we just wait hopefully
            if ($launch_result =~ /You have requested more instances \(\d+\) than your current instance limit of (\d+) allows for the specified instance type/) {
              print " Instance Type Limit - skipping\n";
              $t1 = $1;
              if ($limit->{$type}->{$zone} != $t1) {
                print "Limit mismatch, specified ", $limit->{$type}->{$zone}, " for $type but AWS reports the limit is $t1; updating\n";
                foreach $t2 (keys %{$limit->{$type}}) {
                  $limit->{$type}->{$zone} = $t1;
                }
              }
              $hit_limit_type = 1;
            } elsif ($launch_result =~ /Your quota allows for 0 more running instance/) {
              print " Total Instance Limit - waiting $query_wait seconds (set with -wait option)\n";
              $hit_limit_instance = 1;
            }
            last;
          } elsif ($launch_result =~ /An error occurred \(InsufficientInstanceCapacity\)/) {
              # we can hopefully just shift to another zone
            print " No more instances of $type in $zone\n";
            # when both are 0, should just keep going and go to the next zone
            $hit_limit_type = 0;
            $hit_limit_instance = 0;
            last;
          } elsif ($launch_result =~ /An error occurred \(VolumeLimitExceeded\)/) {
            print " Storage limit exceeded - waiting 20 min to try again\n";
            sleep 20*60;
            last;
          } else {
            print " ERROR!\n";
            print STDERR $launch_result;
            print STDERR "\n\nAborting. Ended with error on index $i (can rerun with -start $i)\n";
            exit;
          }
        }
      }
      last if $hit_limit_instance || $hit_limit_type;
    }
    # this loop will iterate through types
    # so if it's an instance type limit we can still try the next one
    last if $hit_limit_instance;
  }
  system("date");
  sleep $query_wait;
}

sub get_status {
  my $list = `$STATUS`;
  # should get
  # 0 Name
  # 1 InstanceID
  # 2 InstanceType
  # 3 Zone
  # 4 State
  # 5 ImageID
  # 6 PublicIP
  # 7 LaunchTime
  my @l = split /\n/, $list;
  my $r;
  my @f;
  my $i;
  my $running = 0;
  foreach $i (@l) {
    @f = split /\t/, $i;
    $r->{$f[4]}->{$f[2]}->{$f[3]} = 0 if !defined $r->{$f[4]}->{$f[2]}->{$f[3]};
    $r->{$f[4]}->{$f[2]}->{$f[3]}++;
    $running++ if $f[4] eq "running";
  }
  $r->{__TOTAL__} = $running;
  return $r;
}
