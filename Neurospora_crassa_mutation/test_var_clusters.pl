#!/usr/bin/perl -w
#
#   test_var_clusters.pl -- Test how many variants are in a defined cluster
#
#
#   Author: Nowind
#   Created: 2011-10-19
#   Updated: 2019-07-02
#   Version: 2.0.1
#
#   Change logs:
#   Version 1.0.0 18/11/19: The initial version.
#   Version 1.0.1 18/12/18: Update: Use groups as group id.
#   Version 1.0.2 18/12/19: Update: Add stats for mutation spectra.
#   Version 2.0.0 19/01/01: Update: Add simulation tests.
#   Version 2.0.1 19/07/02: Update: Remove unused codes; update descriptions of options.




use strict;

use Data::Dumper;
use Data::Random qw(:all);
use Getopt::Long;
use File::Find::Rule;
use Statistics::Descriptive;
use List::Util qw( sum max );


use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);
use MyPerl::Statistics;

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
   $options{max_distance} = 1000;
   $options{min_size}     = 1;
   $options{rand_times}   = 1;
GetOptions(
            "query-file=s"              => \$options{query_file},
            "context-freq=s"            => \$options{context_freq},
            "output=s"                  => \$options{output},
            "distance=i"                => \$options{max_distance},
            "size=i"                    => \$options{min_size},
            "fasta-file=s"              => \$options{fasta},
            "times=i"                   => \$options{rand_times},
           );

unless( $options{query_file} ) {
    print <<EOF;

$0  -- Test how many variants are in a defined cluster

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -q, --query-file  <filename>
        file contain observed mutation sites in VCF format, required
    -o, --output      <filename>
        output file, default to STDOUT

    -d, --distance <int>
        maximum distance of two variants to be considered as "clustered"
        [default: 1000 (bp)]
    -s, --size     <int>
        minimum cluster size (number of variants within a cluster) for
        output [default: 1]
    
    
    -r, --ref-seq  <filename>
        give a reference genome file with unwanted regions masked with "N"
        to simulate same amount of mutations with same base substitution
        patterns, required for simulation analysis
    
    -t, --times    <int>
        random times [default: 1]

EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($options{output}) {
    open (STDOUT, ">", "$options{output}") or die $!;
}

print STDERR "\r>> Start processing ... ";

my $fh       = getInputFilehandle($options{query_file});
my %var_positions = ();
my %var_sources   = ();
my %spectra       = ();
my $total_mut     = ();
my %group_counts  = ();
my %ref_alleles   = ();
while (<$fh>)
{
    next if (/^\#/ || /^\s+$/);
    
    my ($chr, $pos, $groups, $ref, $alt, $qual, $filter, $info, $format, $sample) = (split /\s+/);
    
    my $group_id = $groups;  ## obtain the tetrad/ascus id
    my ($source) = ($info =~ /Parental_Source=(\w+)/);  ## obtain the parental source
    
    push @{$var_positions{$group_id}->{$chr}}, $pos; ## save positions

    $group_counts{$group_id} ++;
    $total_mut               ++;
    $ref_alleles{$ref}       ++; ## count pre-mutated alleles
    
    unless($options{fasta}) {
        $var_sources{$group_id}->{$chr}->{$pos} = $source;
        $spectra{$group_id}->{$chr}->{$pos}    = "$ref\>$alt";  ## save mutation orientations
    }
}

print "$HEADER##" . (scalar localtime()) . "\n";

if ($options{fasta}) { ## make simulation on specified regions
    ###print "#source\tgroup_id\tcluster_id\tchrom\tstart\tend\tlength\tsize\n";
    print "#source\tno_of_clusters\tlength(mean)\tlength(stdev)\tvars_in_clusters\tvars(mean)\tvars(stdev)\n";
}
else {  ## only take statistics from observed results
    print "#group_id\tcluster_id\tchrom\tstart\tend\tlength\tsize" .
          "\tvar_sources:P1\tP2\tBP_Region\tUncertain\tspectra\tmajor_change\tmajor_perc\n";
}

if ($options{fasta}) {
    cluster_positions(\%var_positions, \%options, "Observed");
}
else {
    cluster_positions(\%var_positions, \%options, \%var_sources);
}

print STDERR "\tdone!\n";


###print STDERR Dumper(%ref_alleles); exit;
###print STDERR Dumper(%group_ids); exit;

##
## simulation analysisread into all sequences
##
my %input_seqs = ();
if ($options{fasta}) {
    ##
    ## read into sequences with unwanted regions masked by "Ns"
    ##
    print STDERR ">> Start reading $options{fasta} ... ";
    parse_fasta_SEQs(\%input_seqs, $options{fasta});
    print STDERR "done!\n";
    
    ##
    ## sampling sequences to array of positions, only nucleotieds
    ## from the pre-mutated allels would be used, for instance,
    ## if only G->A mutations are given in query file, only G sites
    ## will be sampled
    ##
    print STDERR ">> Start preparing simulation sites ... ";
    my $ra_nt_sets = sampling_seqs(\%input_seqs, \%ref_alleles);
    print STDERR "done!\n";
    
    my @group_ids  = sort keys %group_counts;
    
    ##
    ## make simulations
    ##
    my $i=0;
    while (++$i <= $options{rand_times})
    {
        print STDERR "\r>> Start simulating ... $i/$options{rand_times}";
        
        my %rand_vars = ();
        for my $gid (@group_ids)
        {
            ## randomly pick genomic sites with a same amount of the real
            ## mutations in a single genome
            my @random_set = rand_set( set => $ra_nt_sets, size => $group_counts{$gid} );
            
            for (my $j=0; $j<@random_set; $j++)
            {
                my ($chr, $pos) = (split /\:/, $random_set[$j]);
                
                push @{$rand_vars{$gid}->{$chr}}, $pos;
            }
        }
        
        cluster_positions(\%rand_vars, \%options, "Rand$i");
    }
    
    print STDERR "\tdone!\n";
}


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

=head2 sampling_seqs

    About   : Query sequence data.
    Usage   : sampling_seqs($rh_seqs, $ra_group_ids, $rh_alleles);
    Args    : Hash reference to fasta sequences;
              Array reference to group ids;
              Hash reference to pre-mutated alleles
    Returns : Array reference to all sampled positions

=cut
sub sampling_seqs
{
    my ($rh_seqs, $rh_alleles) = @_;
    
    my @nts_all = ();
    for my $chrom (sort keys %{$rh_seqs})
    {
        my $seq = uc $rh_seqs->{$chrom};
        
        my @nts = split //, $seq;
        
        for (my $i=0; $i<@nts; $i++)
        {
            next unless ( $rh_alleles->{$nts[$i]} && $rh_alleles->{$nts[$i]} > 0 );
            
            next if (($nts[$i] eq "N") || ($nts[$i] eq "n"));
            
            my $pos = $i + 1;
            push @nts_all, "$chrom:$pos";
        }
    }
    
    return \@nts_all;
}




=head2 cluster_positions

    About   : Query sequence data.
    Usage   : sampling_seqs($rh_seqs, $ra_group_ids, $rh_alleles);
    Args    : Hash reference to variant sites;
              Hash reference to options;
              Additional infos
    Returns : Null

=cut
sub cluster_positions
{
    my ($rh_variants, $opts, $addition) = @_;
    
    my @clust_lengths = ();
    my @clust_sizes   = ();
    for my $gid (sort keys %{$rh_variants})
    {
        my %clusters = ();
        
        for my $chrom (sort keys %{$rh_variants->{$gid}})
        {
            my @var_pos  = sort {$a <=> $b} @{$rh_variants->{$gid}->{$chrom}};
            my $clust_id = '00000';
            
            push @{$clusters{$chrom}->{$clust_id}->{pos}}, $var_pos[0];
            
            unless ($opts->{fasta}) {
                my $var_src = $addition->{$gid}->{$chrom}->{$var_pos[0]};
                
                $clusters{$chrom}->{$clust_id}->{src}->{$var_src} ++;
            }
            
            for (my $i=1; $i<@var_pos; $i++)
            {
                if ($var_pos[$i] - $var_pos[$i-1] > $opts->{max_distance}) {
                    ++ $clust_id;
                }
                
                push @{$clusters{$chrom}->{$clust_id}->{pos}}, $var_pos[$i];
                
                unless ($opts->{fasta}) {
                    my $var_src = $addition->{$gid}->{$chrom}->{$var_pos[$i]};
                    
                    $clusters{$chrom}->{$clust_id}->{src}->{$var_src} ++;
                }
            }
        }
        
        my $out_id      = '00000';
        for my $chrom (sort keys %clusters)
        {
            for my $clust_id (sort keys %{$clusters{$chrom}})
            {
                my @clust_sites = @{$clusters{$chrom}->{$clust_id}->{pos}};
                my $clust_size  = scalar @clust_sites;
                
                next unless ($clust_size >= $opts->{min_size});
                
                my $clust_start = $clust_sites[0];
                my $clust_end   = $clust_sites[-1];
                my $clust_len   = $clust_end - $clust_start + 1;
                
                ++ $out_id;
                
                if ($clust_size >= $opts->{min_size} && $clust_size >= 2) {
                    push @clust_lengths, $clust_len;
                    push @clust_sizes,   $clust_size;
                }
                
                if ($opts->{fasta}) {
                    ###print STDOUT "$addition\t$gid\tCL$out_id\t$chrom\t$clust_start\t$clust_end\t$clust_len\t$clust_size\n";
                }
                else {
                    ## sources of variants
                    my @var_src = ();
                    for my $src (qw(P1 P2 BP_Region Uncertain))
                    {
                        my $cnt = $clusters{$chrom}->{$clust_id}->{src}->{$src} ? $clusters{$chrom}->{$clust_id}->{src}->{$src} : 0;
                        push @var_src, $cnt;
                    }
                    my $var_srcs = join "\t", @var_src;
                    
                    ## spectra of variants
                    my %clust_spectra = ();
                    for (my $i=0; $i<@clust_sites; $i++)
                    {
                        my $change = $spectra{$gid}->{$chrom}->{$clust_sites[$i]};
                        $clust_spectra{$change} ++;
                    }
                    
                    my @changes_sort = sort {$clust_spectra{$b} <=> $clust_spectra{$a}} keys %clust_spectra;
                    
                    my $major_type   = $changes_sort[0];
                    my $major_perc   = 100 * $clust_spectra{$major_type} / $clust_size;
                    
                    my @clust_spectra = ();
                    for my $change (sort keys %clust_spectra)
                    {
                        my $count = $clust_spectra{$change};
                        
                        push @clust_spectra, "$change:$count";
                    }
                    
                    my $clust_spectra = join ",", @clust_spectra;
                    
                    ###if ($clust_size == 1) {
                    ###    print STDERR Dumper($clusters{$chrom}->{$clust_id}); exit;
                    ###}
                    
                    print STDOUT "$gid\tCL$out_id\t$chrom\t$clust_start\t$clust_end\t$clust_len"
                               . "\t$clust_size\t$var_srcs\t$clust_spectra\t$major_type\t$major_perc\n";
                }
            }
        }
    }
    
    if ($opts->{fasta}) {
        my $length_stat = Statistics::Descriptive::Full->new();
           $length_stat->add_data(@clust_lengths);
        
        my $length_mean  = $length_stat->mean();
        my $length_stdev = $length_stat->standard_deviation();
        
        my $mut_stat = Statistics::Descriptive::Full->new();
           $mut_stat->add_data(@clust_sizes);
        
        my $mut_mean  = $mut_stat->mean();
        my $mut_stdev = $mut_stat->standard_deviation();
        
        my $clust_count = scalar @clust_lengths;
        my $clust_mut   = sum(@clust_sizes);
        
        print STDOUT "$addition\t$clust_count\t$length_mean\t$length_stdev\t$clust_mut\t$mut_mean\t$mut_stdev\n";
    }
}


