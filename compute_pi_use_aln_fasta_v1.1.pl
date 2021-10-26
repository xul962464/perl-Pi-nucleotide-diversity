#!/usr/bin/env perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my $BEGIN_TIME=time();
my $version="1.0.0";

# ------------------------------------------------------------------
#GetOptions
my $infile;
my $window;
my $step;

GetOptions(
				"help|?" =>\&USAGE,
				"infile:s"=>\$infile,
				"window:s"=>\$window,
				"step:s"=>\$step,
				) or &USAGE;
&USAGE unless ($infile);

#######################################################################################
# ------------------------------------------------------------------
# Main Body
# ------------------------------------------------------------------
$window||=0;
$step||=0;
#get aln seq
my @all = get_aln_seq($infile);

#windows too long
if($window > @{$all[0]}){
	warn "window must less than aln seq length\n";
	exit;
}

#get conserved site 
my ($conserved_site,$for_old_pos) = get_conserved_site(\@all);

#whole pi
my $conserved_len = @{$conserved_site->[0]};
$window ||= $conserved_len;
$step ||= $conserved_len;
my @all_pi_info = computer_pi_with_bin($conserved_site,$conserved_len,$conserved_len);
my ($all_pi,$all_S) = (split/\t/,$all_pi_info[0])[2,3];

#If the entered step size and window are greater than the total length of conservative sites, the total length is calculated
#Same step
if($window > $conserved_len){
	$window = $conserved_len;
}

if($step > $conserved_len){
	$step = $conserved_len;
}

if($window < $step){
	warn" window must large than step\n";
	exit;
}else{
	#Output general information
	print"#infile: ",abs_path($infile),"\n";
	print"#seq number: ",scalar(@all),"\n";
	print"#aln length: ",scalar(@{$all[0]}),"\n";
	print"#conserved length: $conserved_len\n";
	print"#window length: $window\n";
	print"#step size: $step\n";
	print"#Nucleotide diversity, Pi: $all_pi\n";
	print"#Number of polymorphic (segregating) sites, S: $all_S\n";
	print"=" x 35,"\n";
}

#Calculate pi value according to window and step
my @pi_info = computer_pi_with_bin($conserved_site,$window,$step);

#Output information
print"Window\tMidpoint\tPi\tS\n";
print "$_\n" for(@pi_info);
print"=" x 35,"\n";
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";

#sub
sub computer_pi_with_bin{
	my $all_seq = shift;
	my $window = shift;
	my $step = shift;
	
	my $len = @{$all_seq->[0]};
	my $loop = int(($len-$window)/$step);
	my @info;
	for my $d(0..$loop){
		my @bin_site;
		for my $i(0..@$all_seq-1){
			for my $j(0..$window-1){
				$bin_site[$i][$j] = $all_seq->[$i][$d*$step+$j];
			}
		}
		my $start_new = ($d*$step+1);
		my $end_new = $start_new + $window - 1;
		my $mid_new = int(($end_new+$start_new)/2);
		my $mid = $for_old_pos->{$mid_new-1} + 1;
		my $start_old = $for_old_pos->{$start_new-1} + 1;
		my $end_old = $for_old_pos->{$end_new-1} + 1;
		my ($pi,$nu_mutations) = get_conserved_sites_pi(\@bin_site);
		$pi = sprintf"%.5f",$pi;
		$start_old = 1 if($d == 0);
		push @info,"$start_old-$end_old\t$mid\t$pi\t$nu_mutations";
	}

	#Finally there are remaining sequences, less than one window
	if($loop < ($len-$window)/$step){
		my $end_start = ($loop+1)*$step;
		my @bin_site2;
		
		for my $i(0..@$all_seq-1){
			my $m = 0;
			for my $j($end_start..@{$all_seq->[0]}-1){
				$bin_site2[$i][$m++] = $all_seq->[$i][$j];
			}
		}
		my $start_new = ($loop+1)*$step+1;
		my $end_new = @{$all_seq->[0]};
		my $mid_new = int(($end_new+$start_new)/2);
		my $mid = $for_old_pos->{$mid_new-1} + 1;
		my $start_old = $for_old_pos->{$start_new-1} + 1;
		my $end_old = $for_old_pos->{$end_new-1} + 1;
		my ($pi,$nu_mutations) = get_conserved_sites_pi(\@bin_site2);
		$pi = sprintf"%.5f",$pi;

		push @info,"$start_old-$end_old\t$mid\t$pi\t$nu_mutations";
	}
	return @info;
}


sub get_conserved_sites_pi{
	my $all_seq = shift;

	my $sample_nu = @$all_seq;
	my $sum_diff_sites = 0;
	my $len = @{$all_seq->[0]};
	my $nu_mutations = 0;

	die"sample nu less than 2" if($sample_nu <= 1);

	for my $i(0..@{$all_seq->[0]}-1){
		my @sites;
		for my $j(0..@$all_seq-1){
			push @sites , $all_seq->[$j][$i];
		}
		$sum_diff_sites += compute_diff_site(@sites);
		$nu_mutations++ if(compute_diff_site(@sites));
	}
	my $di = $sample_nu * ($sample_nu - 1) / 2;
	my $pi = $sum_diff_sites / $di / $len;
	
	return ($pi,$nu_mutations);
}


sub get_conserved_site{
	my $all_seq = shift;
	my $sample_nu = @$all_seq;
	my @conserved_site;

	my $i2 = 0;
	my %for_old_pos;	#It is used to record the correspondence between the new sequence position and the original sequence position after removing the indel site.

	for my $i(0..@{$all_seq->[0]}-1){
		
		my $indel_flag = 0;
		for my $j(0..@$all_seq-1){
			unless($all_seq->[$j][$i] =~ /[ATGC-]/i){	#Degenerate bases are not supported
				warn "some site is not ATGC in conserved site: $all_seq->[$j][$i](or \L$all_seq->[$j][$i])\n";
#				exit;
			}

			if($all_seq->[$j][$i] eq "-"){
				$indel_flag = 1;
			}
		}
		if($indel_flag == 0){
			$conserved_site[$_][$i2] =  $all_seq->[$_][$i] for (0..@$all_seq-1);
			$for_old_pos{$i2} = $i;
		}else{
			$for_old_pos{$i2} = $i;
			$i2--;
		}
		$i2++;
	}
	return (\@conserved_site,\%for_old_pos);
}


sub compute_diff_site{
	my @sites = @_;

	my $s = 0;
	for (my $i = 0;$i < @sites - 1;$i++) {
		for (my $j = $i+1;$j < @sites;$j++) {
			$s++ if($sites[$i] ne $sites[$j]);
		}
	}
	return $s;
}


sub get_aln_seq{
	my $infile = shift;

	open IN,"$infile" or die"$!:$infile can not open!";
	$/ = ">";<IN>;

	my @all_seq;

	while (<IN>) {
		chomp;
		my ($id,$seq) = (split/\n/,$_,2);
		$seq =~ s/\n|\s//g;
		$seq = "\U$seq";
		$seq =~ s/N/-/g;

		my @seq = split//,$seq;

		push @all_seq,\@seq;
	}
	$/="\n";
	return @all_seq;
}


sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: xul<xul\@genepioneer.cn> <1275875706\@qq.com>
Description:

	This script is used to calculate the  Nucleotide diversities(PI) of aligned sequences.
	The calculation result is consistent with dnasp5.

Usage:
  Options:
	-i --infile		input aln file,fasta format
	-s --step		step	default:length of seq
	-w --window		window	default:length of seq
	-h --help			Help

USAGE
	print $usage;
	exit;
}
