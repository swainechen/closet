#!/usr/bin/perl
#
use warnings;
use strict;
use JSON;
use Data::Dumper;
use Getopt::Long;
&Getopt::Long::Configure("pass_through");

my $profile = "default";
GetOptions (
  'profile=s' => \$profile
);

my $json = `aws --profile $profile ec2 describe-instances`;
my $parsed = JSON::decode_json($json);
my $i;
my $j;
my $k;
my $r;
my $AMI;
my $instance_ID;
my $instance_type;
my $launchtime;
my $name;
my $pub_IP;
my $vpc_ID;
my $state;
my $zone;
my $placement_ref;
my $state_ref;
my $tag_ref;
my $inst_array;
my $res_array = $parsed->{Reservations};
my @out = ();

# output everything as on EC2 console
foreach $i (0..$#{$res_array}) {
  # now drop to instances, index 0 to get data
  $inst_array = $parsed->{Reservations}->[$i]->{Instances};
  foreach $j (0..$#{$inst_array}) {
    $r = $inst_array->[$j];
    # now collect info
    $vpc_ID = $r->{VpcId};
    if (defined $r->{PublicIpAddress}) {
      $pub_IP = $r->{PublicIpAddress};
    } else {
      $pub_IP = "";
    }
    $AMI = $r->{ImageId};
    $instance_type = $r->{InstanceType};
    $instance_ID = $r->{InstanceId};
    $launchtime = $r->{LaunchTime};
    # deal with refs
    # get name out of tag - this is an array ref, each element a hash ref
    $tag_ref = $r->{Tags};
    $name = "";
    foreach $k (0..$#{$tag_ref}) {
      if ($tag_ref->[$k]->{Key} eq "Name") {
        $name = $tag_ref->[$k]->{Value};
        last;
      }
    }
    # get running/stopped etc, this is a hash with state in the Name key 
    $state_ref = $r->{State};
    if (defined $state_ref->{Name}) {
      $state = $state_ref->{Name};
    } else {
      $state = "";
    }
    # get availability zone, this is a hash and should have AvailabilityZone key
    $placement_ref = $r->{Placement};
    if (defined $placement_ref->{AvailabilityZone}) {
      $zone = $placement_ref->{AvailabilityZone};
    } else {
      $zone = "";
    }
  }
  push @out, join ("\t", $name, $instance_ID, $instance_type, $zone, $state, $AMI, $pub_IP, $launchtime);
}
print "# ", join ("\t", qw(Name InstanceID InstanceType Zone State ImageID PublicIP LaunchTime)), "\n";
print join ("\n", @out), "\n";

# print "for hash top->{Reservations}->[1]->{Instances}->[0]\n";
# $ref = $parsed->{Reservations}->[1]->{Instances}->[0];
# print join (", ", keys %$ref), "\n";
# print Dumper $ref;

# seem to get this structure:
#   $parsed->{Reservations}->[0]->{Instances}->[0]
# each machine gets a new array element under Reservations (Instances usually just one element)
# that is a hash ref for each element of the Instances array ref
# Then get these keys:
#   AmiLaunchIndex - scalar (0)
#   Architecture - scalar (x86_64)
#   BlockDeviceMappings - array ref
#   CapacityReservationSpecification - hash key CapacityReservationPreference (open)
#   ClientToken - scalar
#   CpuOptions - hash, keys ThreadsPerCore, CoreCount
#   EbsOptimized - complex bless
#   EnaSupport - maps back to BlockDeviceMappings
#   Hypervisor - scalar (xen)
#   ImageId - scalar
#   InstanceId - scalar
#   InstanceType - scalar
#   KeyName - scalar
#   LaunchTime - scalar
#   Monitoring - hash (State => disabled)
#   NetworkInterfaces - array, each element then a hash
#   Placement - hash, keys AvailabilityZone, GroupName, Tenancy
#   PrivateDnsName - scalar
#   PrivateIpAddress - scalar (dotted quad)
#   ProductCodes - array
#   PublicDnsName - scalar
#   PublicIpAddress - scalar (dotted quad)
#   RootDeviceName - scalar
#   RootDeviceType - scalar (ebs)
#   SecurityGroups - array, then each element a hash with keys GroupId and GroupName
#   SourceDestCheck - maps back to BlockDeviceMappings
#   State - hash, keys Name (running), Code (16)
#   StateTransitionReason - scalar
#   SubnetId - scalar
#   Tags - array then hash with keys Key, Value
#   VirtualizationType - scalar (hvm)
#   VpcId - scalar
