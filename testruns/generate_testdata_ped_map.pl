#!/usr/bin/perl
# generate_testdata.pl Generates test data for imputation from panel
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2017 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

generate_testdata.pl - Generates test data for imputation from panel

=head1 SYNOPSIS

generate_testdata.pl [options]

 Options:
   --hap haplotype file
   --sample sample file
   --legend legend file
   --ped-out output ped file
   --map-out output map file
   --skip number of alleles to skip (Default 50)
   --debug, -d debugging level (Default 0)
   --help, -h display this help
   --man, -m display manual

=head1 OPTIONS

=over

=item B<--hap>

Haplotype file (1000GP_Phase3_chr7.hap.gz)

=item B<--sample>

Sample file (1000GP_Phase3.sample)

=item B<--legend>

Legend file (1000GP_Phase3_chr7.legend.gz)

=item B<--ped-out>

PED output file

=item B<--map-out>

MAP output file

=item B<--skip>

Skip n SNPs before including another in map-out, default 50

=item B<--debug, -d>

Debug verbosity. (Default 0)

=item B<--help, -h>

Display brief usage information.

=item B<--man, -m>

Display this manual.

=back

=head1 EXAMPLES

generate_testdata.pl

=cut

use Scalar::Util qw(looks_like_number);
use vars qw($DEBUG);

my %options = (skip            => 50,
               debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'skip=i',
           'legend=s@',
           'hap=s@',
           'map_out|map-out=s',
           'sample=s',
           'ped_out|ped-out=s',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;

for my $opt (qw(legend hap sample ped_out map_out)) {
    if (not defined $options{$opt}) {
        my $a = $opt;
        $a =~ s/_/-/g;
        push @USAGE_ERRORS,"You must provide the --".$a." option";
    }
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

# output the correct sample file for imputation
my $sample = open_compressed_file($options{sample});
my @samples = read_sample_file($sample);
close($sample);
open(my $map_out,'>:encoding(UTF-8)',$options{map_out});
open(my $ped_out,'>:encoding(UTF-8)',$options{ped_out});


# read in the legend and haps file, and skip $options{skip} SNPs
for my $i (0..$#{$options{legend}}) {
    read_in_legend_hap($options{legend}[$i],
                       $options{hap}[$i],
                       $options{map_out},
                       \@samples,
                      );
}
close($map_out);

output_ped($ped_out,\@samples);
close($ped_out);

sub output_ped {
    my ($ped,$samples) = @_;
    foreach my $sample (@{$samples}) {
        print {$ped} $sample ."\n";
    }
}


sub read_in_legend_hap {
    my ($l,$h,$map,$samples) = @_;
    my $legend = open_compressed_file($l);
    my $hap = open_compressed_file($h);
    my @l_header = split_line($legend);
    my $this_chr;
    my @map_hold;
    while (1) {
        my @l = split_line($legend);
        last unless @l;
        my %l;
        @l{@l_header} = @l;

        my @h = split_line($hap);
        last unless @h;

        my ($rsid) = $l{id} =~ /^(rs\d+|[^:]+:\d+):/;
        if (not defined $this_chr) {
            ($this_chr) = $l{id} =~ /^(\d+|[XY]):/;
        }
        if ($rsid =~ /^rs/) {
            $l{id} = $rsid;
        }
        push @map_hold,"$l{id}\t0\t$l{position}";
        if (defined $this_chr) {
            for my $m_line (@map_hold) {
                print {$map_out} "$this_chr\t$m_line\n";
            }
            @map_hold = ();
        }
        for my $i (0..(@h/2-1)) {
            my $g = -9;
            if (looks_like_number($h[$i*2]) and
                looks_like_number($h[$i*2+1])) {
                $g = $h[$i*2]+$h[$i*2+1];
            }
            if ($g == 0) {
                $samples[$i].=" $l{a0} $l{a0}";
            } elsif ($g == 1) {
                $samples[$i].=" $l{a0} $l{a1}";
            } elsif ($g == 2) {
                $samples[$i].=" $l{a1} $l{a1}";
            } else {
                $samples[$i].=" 0 0";
            }
        }

        last unless skip_n_lines($hap,$options{skip});
        last unless skip_n_lines($legend,$options{skip});
    }
}

sub read_sample_file {
    my ($sample) = @_;

    # throw out the first line
    split_line($sample);
    my @samples;
    my $i = 0;
    while (<$sample>) {
        chomp;
        my ($id,$pop,$group,$sex) = split /\s+/,$_;
        if ($sex eq 'male') {
            $sex = 1;
        } elsif ($sex eq 'female') {
            $sex = 2;
        } else {
            $sex = 0;
        }
        push @samples,"$id $id 0 0 $sex 0";
    }
    return @samples;
}

sub split_line {
    my ($fh) = @_;
    my $a = <$fh>;
    return () unless defined $a;
    chomp $a;
    return split /\s+/,$a;
}

sub skip_n_lines {
    my ($fh,$n) = @_;
    for (1..$n) {
        my $a = <$fh>;
        return 0 if not defined $a;
    }
    return 1;
}

sub open_compressed_file {
    my ($file) = @_;
    my $fh;
    my $mode = '<:encoding(UTF-8)';
    my @opts;
    if ($file =~ /\.gz$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'gzip','-dc';
    }
    if ($file =~ /\.xz$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'xz','-dc';
    }
    if ($file =~ /\.bz2$/) {
        $mode = '-|:encoding(UTF-8)';
        push @opts,'bzip2','-dc';
    }
    open($fh,$mode,@opts,$file) or
        die "Unable to open $file for reading: $!";
    return $fh;
}




__END__
