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
   --sample-out output sample file
   --gens-out genotype output file
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

=item B<--sample-out>

Sample output file

=item B<--gens-out>

Genotype output file

=item B<--skip>

Skip n SNPs before including another in gens-out, default 50

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
           'legend=s',
           'hap=s',
           'gens_out|gens-out=s',
           'sample=s',
           'sample_out|sample-out=s',
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;

for my $opt (qw(legend hap gens_out sample sample_out)) {
    if (not defined $options{$opt}) {
        my $a = $opt;
        $a =~ s/_/-/g;
        push @USAGE_ERRORS,"You must provide the --".$a." option";
    }
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;


# output the correct sample file for imputation
my $sample = open_compressed_file($options{sample});
open(my $sample_out,'>:encoding(UTF-8)',$options{sample_out});
read_and_write_sample_file($sample,$sample_out);
close($sample);
close($sample_out);

# read in the legend and haps file, and skip $options{skip} SNPs, then
# output 1.
my $legend = open_compressed_file($options{legend});
my $hap = open_compressed_file($options{hap});
open(my $gens_out,'>:encoding(UTF-8)',$options{gens_out});
my @l_header = split_line($legend);



while (1) {
    my @l = split_line($legend);
    last unless @l;
    my %l;
    @l{@l_header} = @l;

    my @h = split_line($hap);
    last unless @h;

    my ($rsid) = $l{id} =~ /^(rs\d+|[^:]+:\d+):/;
    # output gens file
    print {$gens_out} "$l{id} $rsid $l{position} $l{a0} $l{a1}";
    for my $i (0..(@h/2-1)) {
        my $g = -9;
        if (looks_like_number($h[$i*2]) and
            looks_like_number($h[$i*2+1])) {
            $g = $h[$i*2]+$h[$i*2+1];
        }
        if ($g == 0) {
            print {$gens_out} " 1 0 0";
        } elsif ($g == 1) {
            print {$gens_out} " 0 1 0";
        } elsif ($g == 2) {
            print {$gens_out} " 0 0 1";
        } else {
            print {$gens_out} " 0 0 0";
        }
    }
    print {$gens_out} "\n";

    last unless skip_n_lines($hap,$options{skip});
    last unless skip_n_lines($legend,$options{skip});
}

sub read_and_write_sample_file {
    my ($sample,$sample_out) = @_;

    print {$sample_out} "ID_1 ID_2 missing sex status\n";
    print {$sample_out} "0 0 0 D B\n";
    # throw out the first line
    split_line($sample);
    my $i = 0;
    while (<$sample>) {
        chomp;
        my ($id,undef) = split /\s+/,$_;
        print {$sample_out} $i++." ".$id. " 0 NA NA\n";
    }
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
