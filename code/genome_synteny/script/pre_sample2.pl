#!/usr/bin/perl -w
use strict;
use Cwd;
use List::Util qw/max min sum maxstr minstr shuffle/;

#handle the sequence
my $score_file = $ARGV[0];
#my $list_file = $ARGV[1];
my $output_file = $ARGV[1];

#HEK293
my $usage = "This script is to sample the ontology gene list score.
usage: $0 <score_file> <output_file> 
";
die $usage if $#ARGV<1;

open(FILE1, $score_file)||die("open $score_file error!\n");
#open(FILE2, $list_file)||die("open $list_file error!\n");
open(OUT, ">", $output_file);

my $sen = ""; my $seq = ""; my $sna = "";
my %inf = ();
my @sen = (); my @sen1 = (); my @sen2 = ();
my $sid = "";
my $i = 0; my $j = 0; my $r; my $k;

#ENSG00000198888	ENSMUSG00000064341
#while($sen = <FILE2>){
#	chomp($sen);
#	@sen = split(/\t/,$sen);
#	$inf{$sen[0]."_".$sen[1]} = 1;
#}
#close FILE2;
#ENSG00000105352	|	255	244	2	494	32	32	1	63	|	ENSMUSG00000037640	|	298 282	19	588	32	36	1	65	|	0	0	0	113	0	113	2	0	2   15	0	17
#0                  1   2    3  4    5   6   7  8    9  10        11            12   13  14 15   16 17  18  19  20  21  22  23  24   25 26   27 28  29  30  31  32  33
while($sen = <FILE1>){
	chomp($sen);
	@sen = split(/\t/,$sen);
	#print OUT $sen[0],"\t",$sen[11],"\t",$sen[27]-$sen[26],"\t",&TMD($sen[27],$sen[26],$sen[5],$sen[16]),"\t",$sen[33]-$sen[32],"\t",&TMD($sen[33],$sen[32],$sen[9],$sen[20]),"\t",$sen[22],"\t",$sen[23],"\t",&TMD($sen[22],0,$sen[2],$sen[13]),"\t",&TMD($sen[23],0,$sen[3],$sen[14]),"\t",$sen[28],"\t",$sen[29],"\t",&TMD($sen[28],0,$sen[6],$sen[17]),"\t",&TMD($sen[29],0,$sen[7],$sen[18]),"\t",$sen[24],"\t",$sen[25],"\t",&TMD($sen[24],0,$sen[2],$sen[14]),"\t",&TMD($sen[25],0,$sen[3],$sen[13]),"\t",$sen[30],"\t",$sen[31],"\t",&TMD($sen[30],0,$sen[6],$sen[18]),"\t",&TMD($sen[31],0,$sen[7],$sen[17]),"\t1\n";
	print OUT $sen[0],"\t",$sen[11],"\t",$sen[27]-$sen[26],"\t",&TMD($sen[27],$sen[26],$sen[5],$sen[16]),"\t",$sen[33]-$sen[32],"\t",&TMD($sen[33],$sen[32],$sen[9],$sen[20]),"\t",$sen[22],"\t",$sen[23],"\t",&TMD($sen[22],0,$sen[2],$sen[13]),"\t",&TMD($sen[23],0,$sen[3],$sen[14]),"\t",$sen[28],"\t",$sen[29],"\t",&TMD($sen[28],0,$sen[6],$sen[17]),"\t",&TMD($sen[29],0,$sen[7],$sen[18]),"\n";
	#print OUT $sen[0],"\t",$sen[11],"\t",$sen[27]-$sen[26],"\t",($sen[27]-$sen[26])/min($sen[5],$sen[16]),"\t",$sen[33]-$sen[32],"\t",($sen[33]-$sen[32])/min($sen[9],$sen[20]),"\t",$sen[22],"\t",$sen[23],"\t",$sen[22]/min($sen[2],$sen[13]),"\t",$sen[23]/min($sen[3],$sen[14]),"\t",$sen[28],"\t",$sen[29],"\t",$sen[28]/min($sen[6],$sen[17]),"\t",$sen[29]/min($sen[7],$sen[18]),"\t1\n";
	#print OUT $sen[0],"\t",$sen[11],"\t",join("\t",@sen[2..9]),"\t",join("\t",@sen[13..20]),"\t",join("\t",@sen[22..33]),"\t1\n";
	$i = $i + 1;
}
# 
print min(1,2),"\t",min(10,2),"\t",min(32,12),"\n";

sub TMD{
	my $a = shift; my $b = shift; my $c = shift; my $d = shift;
	my $whh = 0.0;
	if (min($c, $d) eq 0){
		$whh = 0.0;
	}else{
		$whh = ($a - $b)/min($c, $d);
	}
	return $whh;
}

close FILE1;
close OUT;
exit;


