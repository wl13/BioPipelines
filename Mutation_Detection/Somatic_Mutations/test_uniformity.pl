#!/usr/bin/perl -w
#
#   test_uniformity.pl -- Simulate D = Variance / mean
#
#
#   Author: Nowind
#   Created: 2011-10-19
#   Updated: 2018-10-30
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 18/10/27: The initial version.
#   Version 1.0.1 18/10/30: Updated: Ouput both smaller and greater p-values.




use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use Statistics::Descriptive;
use List::Util qw( sum );


use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);
use MyPerl::Statistics;

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
   $options{times}   = 1;
GetOptions(
            "query-file=s"              => \$options{query_file},
            "output=s"                  => \$options{output},
            "times=i"                   => \$options{times}
           );

unless( $options{query_file} ) {
    print <<EOF;

$0  -- Simulate D = Variance / mean.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -q, --query-file  <string>
        file contain list of numbers, required
    -o, --output <filename>
        output file, default to STDOUT

    -t, --times
        random times [default: 1]

EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($options{output}) {
    open (STDOUT, ">", "$options{output}") or die $!;
}

my $fh   = getInputFilehandle($options{query_file});
my @obs_nums = ();
while (<$fh>)
{
    next if (/^\#/ || /^\s+$/);
    
    chomp;
    
    my $number = (split /\s+/)[0];
    
    push @obs_nums, $number;
}

###print STDERR "obs: @obs_nums\n";

$options{total}        = sum(@obs_nums);
$options{no_of_groups} = scalar @obs_nums;

my $obs_stat = Statistics::Descriptive::Full->new();
   $obs_stat->add_data(@obs_nums);

$options{obs_D} = $obs_stat->variance() / $obs_stat->mean();

test_uniformity(\%options);


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


sub test_uniformity
{
    my ($opts) = shift;
    
    ##
    ## random simulation
    ##
    print "$HEADER##" . (scalar localtime()) . "\n";
    print "##time simulated: $opts->{times}\n";
    
    my @rand_values = ();
    my $i = 0;
    while (++$i <= $opts->{times})
    {
        print STDERR "\r>> Start random process ... duplicate $i\/$opts->{times}";        
        
        ## https://stackoverflow.com/questions/22380890/generate-n-random-numbers-whose-sum-is-m-and-all-numbers-should-be-greater-than
        ## Generate $n numbers between 0 (incl) and 1 (excl).
        my @groups_rand = ();
        for (1..$opts->{no_of_groups})
        {
           push @groups_rand, rand();
        }
        
        ## Calculate factor.
        my $factor = ($opts->{total}) / sum(@groups_rand);
        
        for my $j (0..$#groups_rand)
        {
           $groups_rand[$j] = int($groups_rand[$j] * $factor);
        }
        
        ## Handle loss of fractional component.
        my $fudge = $opts->{total} - sum(@groups_rand);
        for (1..$fudge) {
           ## Adds one to a random number.
           ++$groups_rand[rand(@groups_rand)];
        }
        
        my $groups_rand = join "\t", @groups_rand;
        print STDOUT "Rand$i\t$groups_rand\n";
        
        my $stat = Statistics::Descriptive::Full->new();
           $stat->add_data(@groups_rand);
       
        push @rand_values, $stat->variance() / $stat->mean();
    }
    print STDERR "\tdone!\n";
    
    my $rand_values = join "\t", @rand_values;
    print STDOUT "Rand:$rand_values\n";
    
    my $above_observed = grep { $_ >= $opts->{obs_D} } @rand_values;
    my $below_observed = grep { $_ <= $opts->{obs_D} } @rand_values;
    
    #my $rand_average   = sum(@rand_values) / $opts->{times};
    my $rand_stat = Statistics::Descriptive::Full->new();
       $rand_stat->add_data(@rand_values);
           
    my $rand_average = $rand_stat->mean();
    my $rand_stdev   = $rand_stat->standard_deviation();
           
    # Unbiased estimation of empirical P (expected type I error rate) = (n + 1)/(m + 1)
    # n: the number of observations as or more extreme than that observed in the real test reporting statistic
    # m: the number of randomization (A note on the calculation of empirical P values from Monte Carlo procedures)
    my $p_value_greater  = ($above_observed + 1) / ($opts->{times} + 1);  ## significance level of observed < expected
    my $p_value_smaller  = ($below_observed + 1) / ($opts->{times} + 1);  ## significance level of observed > expected
    
    print "Observed: $opts->{obs_D}\tRand(mean): $rand_average($rand_stdev)\tP-value(smaller): $p_value_smaller\tP-value(greater): $p_value_greater\n";
}


