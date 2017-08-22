#!/usr/bin/perl

#The purpos of this script is to select Orthomcl clusters according to certain criteria and produce a multifasta files containing the protein sequences for every cluster meeting the criteria
#Christoph Hahn, 2014
#copied from https://github.com/chrishah/phylog/blob/master/scripts/select_clusters_v2.pl on 22 August 2017

use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Long;
use Cwd qw(abs_path);
use List::MoreUtils 'first_index'; 

my $USAGE = 	"\nUsage: ./select_cluster_v2.pl <parameters>\n
		mandatory parameters:
		--groups	<file>		groups.txt output file from mcl pipeline
		--min		<int>		minimal number of taxa contained in the cluster, default: 8
		--max_median	<int>		maximum median number of sequences per taxon, default: 2
		--max_mean	<int>		maximum mean number of sequences per taxon, default: 5
		
		optional parameters:
		--fasta	<file>		create multifastA files for each cluster, extract from protein fasta file specified
		--out	<string>	name of directory to write results (default: cluster_select_output)
		--critical <file>	list of taxa that cluster must contain (if more than one are given per line, one per line must be present)
		--exclusive		keep only critical taxa in output\n\n";
		

my $command = $0;
for (@ARGV){
        $command .= " $_";
}


my $goodproteins;
my $results_directory = "cluster_select_output";
##############################
my $groupsfile;
my $min_taxon_per_cluster=8;
my $max_median_sequences_per_taxon=2;
my $max_mean_sequence_per_taxon=5;
my $exclude;
my @clusters = ();
my @cluster_ID = ();
my %sequences_per_cluster;
my $percent_meeting;
my @output;
my $critical_taxa;
my @critical;
my @critical_ind;

GetOptions (    "min=i" => \$min_taxon_per_cluster,
                "out=s" => \$results_directory,
                "groups=s" => \$groupsfile,
                "critical=s" => \$critical_taxa,
                "max_median=i" => \$max_median_sequences_per_taxon,
                "max_mean=i" => \$max_mean_sequence_per_taxon,
                "fasta=s" => \$goodproteins,
		"exclusive!" => \$exclude) or die "Incorrect usage!\n$USAGE";

if (!$groupsfile){
	print "$USAGE\n\nplease specify mcl output groups file\n\n";
	exit;
}else{
	$groupsfile=abs_path($groupsfile);
	unless (-e $groupsfile){
		print "$USAGE\nplease specify mcl output groups file\n\n";
		exit;
	}
	my $groupsfile_fh = &read_fh($groupsfile);
	while (<$groupsfile_fh>){
		push (@clusters, "$_");
	}
}
if ((!$min_taxon_per_cluster) || (!$max_median_sequences_per_taxon) || (!$max_mean_sequence_per_taxon)){
	print "$USAGE\n\nNot all necessary parameters specified\n\n";
	exit;
}

if ($critical_taxa){
	$critical_taxa=abs_path($critical_taxa);
	unless (-e $critical_taxa){
		print "$USAGE\n\ncritical taxa file specified, but is it there?\n";
		exit;
	}
	@critical=&read_critical($critical_taxa);
	for (@critical){
		my @temp=split(/\t/);
		if (@temp>0){
			for (@temp){
				push(@critical_ind,$_);
			}
		}else {
			push(@critical_ind,$_);
		}
	}
}
if ($goodproteins){
	$goodproteins=abs_path($goodproteins);
	unless (-e $goodproteins){
		print "$USAGE\n\nIs the proteins file there?\n\n";
		exit;
	}
}

print "\n\n". strftime("%b %e %H:%M:%S", localtime) . "\n";
print "\nFull command run: $command\n";
print "\nall parameters seem to make sense:\n";
print "groups-file: $groupsfile\n";
print "minimum number of taxa: $min_taxon_per_cluster\n";
print "maximum median number of sequences per taxon: $max_median_sequences_per_taxon\n";
print "maximum mean number of sequences per taxon: $max_mean_sequence_per_taxon\n";
if ($critical_taxa){
	print "critical taxa specified in: $critical_taxa\n";
	print "as follows:\n";
	for (@critical){
		print "\t$_\n";
	}
}
if ($goodproteins){
	print "proteins file: $goodproteins\n";
	print "\ncreating results directory: $results_directory\n";
        mkdir $results_directory or die $!;
}

for (@clusters){
	chomp;
	my @critical_temp=@critical;
	my @cluster_line = ();
	my @cluster_elements = ();
	my @cluster_elements_ind = ();
	my @cluster_taxa = ();
	undef @output;
	@cluster_line = split /: /;
	@cluster_elements = (split / /, $cluster_line[1]);
	for (@cluster_elements){
		@cluster_elements_ind = split /\|/;
		push (@cluster_taxa, $cluster_elements_ind[0]);
	}
	%sequences_per_cluster = map { $_ => 1} @cluster_taxa;
	
	my %count;
	map { $count{$_}++} @cluster_taxa;
	my $size = scalar keys %count;
	push (@output, "$cluster_line[0] contains a total of " . scalar @cluster_elements . " sequences");
#	print "\n-------\n$cluster_line[0] contains a total of " . scalar @cluster_elements . " sequences\n";
	my @numbers = ();
	foreach (keys %count){
		push (@output, "$_ : $count{$_}");
#		print "$_ : $count{$_}\n";
		push (@numbers, $count{$_});
		my $current = $_;
#		print "\n\ncurrent taxon is: $current\n";
		my $index = first_index { /$current/ } @critical_temp;
		if ($index>=0){
#			print "$current is found at position $index in @critical_temp\n";
#			print "$critical_temp[$index] will be removed from critical\n";
			splice(@critical_temp,$index,1);
		}
		
	}
	my $mean = scalar @cluster_elements / $size;
	my $median = &median(@numbers);

	push (@output, "number of unique taxa: $size");
#	print "number of unique taxa: $size\n";
	push (@output, "mean number of sequences per taxon: $mean");
#	print "mean number of sequences per taxon: $mean\n";
	push (@output, "median number of sequences per taxon: $median");
        
	if ( (scalar @critical_temp == 0) && ($size >= $min_taxon_per_cluster) && ($median <= $max_median_sequences_per_taxon) && ($mean <= $max_mean_sequence_per_taxon)){
		print "-------\n$cluster_line[0] meets the criteria\n";
#		print "Interesting Taxon ($interesting_taxon) is represented by $expansion copies!!!\n";
               	for (@output){
			print "$_\n";
		}
		print "All critical taxa found!\n";
#		print "-------\n";
		push(@cluster_ID, $cluster_line[0]);
		if ($goodproteins){
			chdir $results_directory or die $!;
			open (OUT,">$cluster_line[0].list");
			my $ex_ID=$cluster_line[0].".fasta:";
			for (@cluster_elements){
				if ($exclude){
					my @taxon = split /\|/;
	                                if ($taxon[0] ~~ @critical_ind){
	                                        print OUT $_ . "\n";
	                                }else{
						$ex_ID.=" $_"; 
					}
				}else{
					print OUT $_ . "\n";
				}
			}
			close OUT;
			if ($exclude){
				print "sequences excluded from the fasta file:\n$ex_ID\n";
			}			
			`select_contigs.pl -n $cluster_line[0].list -l 100 $goodproteins $cluster_line[0].fasta`;
			unlink "$cluster_line[0].list";
			chdir "../" or die $!;
        	}
	}else {
			print "-------\n$cluster_line[0] does not meet the criteria\n";
			for (@output){
	                     	print "$_\n";
               		}
			if (scalar @critical_temp > 0){
				print "not all critical taxa found\n";
			}elsif (scalar @critical_temp ==0){
				print "all critical taxa found\n";
			}
	}
	
	undef %{count};
}
$percent_meeting = sprintf ("%.2f",((scalar @cluster_ID / scalar @clusters) * 100));
print "\ninitial number of clusters: " . scalar @clusters . "\n";
print "number of clusters meeting the criteria: " . scalar @cluster_ID . " ($percent_meeting %)\n";
print strftime("%b %e %H:%M:%S", localtime) . "\n\n";

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub read_critical {
	my $file = shift;
	my @array;
	unless (-e $file){
		print "$USAGE\n\ncritical taxon file is specified but seems not to be there\n\n";
		exit;
	}
        open (IN,"<$file") or die $!;
        for (<IN>){
		chomp;
                push (@array, "$_");
        }
	if (scalar @array == 0){
		print "\n\nfile containing critical taxa seems to be empty..\n\n";
	}
	return @array;
}

sub read_fh {
        my $filename = shift @_;
        my $filehandle;
        if ($filename =~ /gz$/) {
                open $filehandle, "gunzip -dc $filename |" or die $!;
        }
        else {
                open $filehandle, "<$filename" or die $!;
        }
        return $filehandle;
}
