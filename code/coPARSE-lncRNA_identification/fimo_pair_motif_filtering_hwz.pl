#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

my ($infile) = $ARGV[0];
my ($pvalue) = $ARGV[1];
my ($outfile) = $ARGV[2];

#perl motif_conserve.pl hum_mou.txt All75k.HM_align_full.header.out phastCons100way motif_phastCons100way_score.txt

my $usage = "This script is to get model motif.
usage: $0 <motif_file> <path> <out_file>
";
die $usage if $#ARGV<2;


my %minf1 = (); my %minf2 = (); my %zinf1 = (); my %zinf2 = ();
open(FILE1, $infile)||die("open $infile error!\n");
my $sen = "";
while($sen = <FILE1>){
	#HM00055 HUM     118     134     +       16.1137 4.27e-07                CCCCAACTCCTACCCCC
	#HM00056 HUM     121     145     -       13.9466 4.52e-06                CCAGGCTGCCCGGGGGTAGGAGTTG
	chomp($sen);
	#$i = $i + 1;
	my @sent1 = split(/\t/,$sen);
	my $sid = $sent1[0]."|".$sent1[4];
	my $sindex = $sent1[2]*1000000 + $sent1[3];
	if($sent1[6] <= $pvalue){
		if($sent1[1] eq "HUM"){
			$minf1{$sid} = 1;
			if(exists $zinf1{$sindex}){
				push(@{$zinf1{$sindex}}, $sen);
			}else{
				$zinf1{$sindex} = [$sen];
			}
		}else{
			$minf2{$sid} = 1;
			if(exists $zinf2{$sindex}){
				push(@{$zinf2{$sindex}}, $sen);
			}else{
				$zinf2{$sindex} = [$sen];
			}
		}		
	}
}
close FILE1;



open(OUT1, ">", $outfile);

foreach my $ke (sort{$a<=>$b} keys %zinf1){
	my @sent1 = @{$zinf1{$ke}};
	foreach my $ke2 (@sent1){
		my @sent2 = split(/\t/,$ke2);
		my $sid = $sent2[0]."|".$sent2[4];
		if((exists $minf1{$sid})&&(exists $minf2{$sid})){
			print OUT1 $ke2,"\n";
		}
	}
}

foreach my $ke (sort{$a<=>$b} keys %zinf2){
	my @sent1 = @{$zinf2{$ke}};
	foreach my $ke2 (@sent1){
		my @sent2 = split(/\t/,$ke2);
		my $sid = $sent2[0]."|".$sent2[4];
		if((exists $minf1{$sid})&&(exists $minf2{$sid})){
			print OUT1 $ke2,"\n";
		}
	}
}

close OUT1;

exit;


