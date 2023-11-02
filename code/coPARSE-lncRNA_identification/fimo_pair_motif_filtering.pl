#!usr/bin/perl -w

#BEGIN {
#  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/ScrDB/scripts";
#}

use Myownperlmodule qw(motifpos);
use List::Util qw(min max);

use strict;
#use experimental qw(smartmatch);

my $infimo=shift;
my $pthreshold=shift;

my %redundancy;
open FIMOrep, $infimo;
while(<FIMOrep>){
     next if (/^#/);
     chomp;
     my @F=split(/\t/,$_);
     $redundancy{$F[0]."\t".$F[4]."\t".$F[1]}=1;
     $redundancy{$F[1]}=0;
}
close FIMOrep;

my @deRedundancy;
foreach my $Fk (keys %redundancy){if ($redundancy{$Fk} == 0){push @deRedundancy, $Fk};}
$redundancy{$deRedundancy[0]} = $deRedundancy[1]; $redundancy{$deRedundancy[1]} = $deRedundancy[0];

#print $redundancy{"HUM"},"---", $redundancy{"MOU"}, "\n";

#print join(" :OK: ", @{[keys %redundancy]}),"\n";

my %hash; my %motif;
open FIMO, $infimo;
while(<FIMO>){
     next if (/^#/);
     chomp;
     my @Ln=split(/\t/,$_);
     if (defined $redundancy{$Ln[0]."\t".$Ln[4]."\t".$redundancy{$Ln[1]}}){
         push @{$hash{$Ln[1]}}, $_;
         push @{$motif{$Ln[0]."\t".$Ln[4]}{"Score"}}, $Ln[5];
         push @{$motif{$Ln[0]."\t".$Ln[4]}{"Pvalue"}}, $Ln[6];
#print $Ln[0],"\t",$Ln[4],"\t<==\n";
     }
}
close FIMO;

my %SelectMotif;
foreach my $mkey (sort keys %motif){
     my $tmpSmax=max(@{$motif{$mkey}{"Score"}});
     my $tmpPmin=min(@{$motif{$mkey}{"Pvalue"}});
     if(($tmpSmax >=13)&&($tmpPmin <=$pthreshold)){
#        print $mkey, "=> OK\t$tmpSmax\t$tmpPmin\n";
        $SelectMotif{$mkey}=1;
     }
}

#print join(" :: ", @{[keys %SelectMotif]}),"\n";

foreach my $key (sort keys %hash){
     my @ntmp=@{$hash{$key}}; my @f_ntmp=();
     for(@ntmp){ 
         my @ln=split(/\t/,$_); my $tmpID= $ln[0]."\t".$ln[4];
         if(defined($SelectMotif{$tmpID})){ 
           push @f_ntmp, $_;
         } 
     }
#     print join("\n", @f_ntmp), "\n";
     my $temp=motifpos(\@f_ntmp, 0.8, 1);
     my @temp=@{$temp};
     print join("\n", @temp), "\n";
}
