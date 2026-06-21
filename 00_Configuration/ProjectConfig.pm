#===============================================================================
# Methods_Scripts::Config
#
# Loads the project config (Config::General .ini) and substitutes the
# __PROJECT_ROOT__ token in every path with $ENV{PROJECT_ROOT} (or the value
# passed via --ProjectRoot=...).
#
# Usage:
#   use lib '00_Configuration';
#   use Config qw(load_config);
#   my %cfg = load_config('--ConfigFile=00_Configuration/config.GRCh38.ini');
#   print $cfg{paths}{repeats_dir};   # e.g. /your/project/geneconv_complete/repeats
#===============================================================================
package ProjectConfig;

use strict;
use warnings;
use Exporter 'import';
use Config::General;
use File::Spec;
use Cwd qw(abs_path);

our @EXPORT_OK = qw(load_config substitute_root parse_args);

# Replace __PROJECT_ROOT__ inside a string (recursively for hashes/arrays).
sub substitute_root {
    my ($val, $root) = @_;
    if (ref($val) eq 'HASH') {
        $_ = substitute_root($_, $root) for values %$val;
    } elsif (ref($val) eq 'ARRAY') {
        $_ = substitute_root($_, $root) for @$val;
    } elsif (defined $val) {
        $val =~ s/__PROJECT_ROOT__/$root/g;
        return $val;
    }
    return $val;
}

# Consume --ConfigFile= and --ProjectRoot= from @ARGV (so Getopt::Long in the
# caller does not see them as unknown options).
sub parse_args {
    my $cfg_file = '00_Configuration/config.GRCh38.ini';
    my $root     = $ENV{PROJECT_ROOT};
    my @keep;
    for my $a (@ARGV) {
        if    ($a =~ /^--ConfigFile=(.+)$/)  { $cfg_file = $1 }
        elsif ($a =~ /^--ProjectRoot=(.+)$/) { $root     = $1 }
        else  { push @keep, $a }
    }
    @ARGV = @keep;
    die "No project root set. Use --ProjectRoot=... or `export PROJECT_ROOT=...`.\n"
        unless defined $root && length $root;
    die "Config file not found: $cfg_file\n" unless -f $cfg_file;
    # Resolve to an absolute path so child jobs that change cwd can still find it.
    $cfg_file = abs_path($cfg_file);
    $root     = abs_path($root) if -d $root;
    return ($cfg_file, $root);
}

sub load_config {
    my ($cfg_file, $root) = parse_args();
    my $cfg_ref = Config::General->new($cfg_file);
    my %cfg     = $cfg_ref->getall;
    substitute_root(\%cfg, $root);
    $cfg{_meta}{config_file}  = $cfg_file;
    $cfg{_meta}{project_root} = $root;
    return %cfg;
}

1;
