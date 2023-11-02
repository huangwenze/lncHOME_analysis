#!/usr/bin/perl -w

#BEGIN {
#  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
#}

#This file was used to filter conserved motif and matched to 3th species.


use Lncrnaevo::Myownperlmodule qw(uniq_array);
#use experimental qw(smartmatch);

#use Parallel::ForkManager;
require '/Share/home/zhangqf/usr/perl/lib/site_perl/5.22.0/Parallel/ForkManager.pm';

$scrdir = "/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
$InScoreRank = shift;
$_fa = shift;
$_dir = shift;
$motifdb= shift;


`mkdir -p ${_dir}/compare_3rd/`;
$compare_3rd="${_dir}/compare_3rd";

`rm -rf $compare_3rd/meme_temp; mkdir -p $compare_3rd/meme_temp`;
($inMEME=$_dir)=~s/RUN_(\w+)\/cluster/RUN_$1\/meme_$1\.out\.update/;
$dn_inMEME=`find $compare_3rd/meme_temp/  -maxdepth 1 -type f -name "meme_all.update"`; chomp $dn_inMEME;
if ( ! -e $dn_inMEME ){
     `cp $inMEME $compare_3rd/meme_temp/meme_all.update`;
     $dn_inMEME=`find $compare_3rd/meme_temp/  -maxdepth 1 -type f -name "meme_all.update"`; chomp $dn_inMEME;
}

$faID=$_fa; $faID=~s/\S+\///; $faID=~s/\.fa//;
$_inMEME=`find $motifdb -maxdepth 1 -type f -name "$faID.motifpool.meme"`; chomp $_inMEME;

open IN, $InScoreRank;
while(<IN>){
     next if (/^#/);
     chomp;
     if (/^\//){
        $line1=$_; $line2=<IN>;
        push @TOT, $line1.":".$line2;
     }
}
close IN;

#$mp = Parallel::ForkManager->new(5);

for $ele (@TOT){

#$mp->start and next;

     @TwoParts=split(":", $ele);
     if ($TwoParts[0]=~/^\//){
     @F=split(/\t/,$TwoParts[0]);

     $DP_results=$F[0]; 
     $fileID=$DP_results;   $fileID=~s/.*\///;  $fileID=~s/_\S+\.align//; 
     $fileID=$fileID.".$F[1]";

     #$myarg="";
     #$myarg=getMotifID($TwoParts[1]);                                     ### Using module &getMotifID, note we should use line 57 in file "FIMO_DP_EUC_ThirdPlus_SVG_v2.pl"!!!

     $myarg=""; @myTwoParts=split("\t", $TwoParts[0]);
     $myarg=getMotifID2($myTwoParts[0], $myTwoParts[4], $myTwoParts[5]);         ### Using module &getMotifID2, note we should use line 58 in file "FIMO_DP_EUC_ThirdPlus_SVG_v2.pl"!!!

     `perl $scrdir/grep_meme_recognized_by_mcast_v1.0.pl $dn_inMEME $myarg >$compare_3rd/${faID}.$fileID.pre.motif`;
     `perl $scrdir/motif_update.pl $compare_3rd/${faID}.$fileID.pre.motif $_fa $compare_3rd/${faID}.$fileID.post.motif` if ( not $? );

     `perl $scrdir/grep_meme_recognized_by_mcast_v1.0.pl $_inMEME $myarg >$compare_3rd/${faID}.$fileID.motif` if ( not $? );

     `cat $compare_3rd/${faID}.$fileID.post.motif >>$compare_3rd/${faID}.$fileID.motif` if ( not $? );

     `/Share/home/zhangqf/usr/meme_4.11.2-app/bin/fimo --thresh 1e-4 --verbosity 1 --text $compare_3rd/${faID}.$fileID.motif $_fa|sort -k2,2 -k3,3n -k4,4n >$compare_3rd/${faID}.${fileID}.fimo` if ( not $? );

     `rm $compare_3rd/${faID}.$fileID.pre.motif | mv $compare_3rd/${faID}.$fileID.post.motif $compare_3rd/meme_temp/` if ( not $? );

     }

     if ( not $? ) {
         print STDERR "$compare_3rd/${faID}.${fileID}.fimo finished\n";
     }

#$mp->finish;

}

#$mp->wait_all_children;

if ( not $? ){
    `cat $compare_3rd/meme_temp/${faID}.*.post.motif >$dn_inMEME; rm -f $compare_3rd/meme_temp/${faID}.*.post.motif`;
}



sub mergeMEME{

     $d_inf1=shift; $d_inf2=shift;
     %hash2=(); @orig_meme=();

     open MEME2, $d_inf2;
     while(<MEME2>){
     chomp;
     if(/^MOTIF\s+(\d+)/){
         $motif=$1;
         $hash2{$motif}.=$_."\n";
     }elsif(/^letter-probability matrix:.+ w=\s+(\d+)/){
            $hash2{$motif}.=$_."\n";
            $width=$1;
            for ($i=1;$i<=$width;$i++){
              $newLine=<MEME2>;
              $hash2{$motif}.=$newLine;
            }
     }
     }
     close MEME2;

     open MEME1, $d_inf1;
     while(<MEME1>){
       chomp;
       push @orig_meme, $_;
     }
     close MEME1;

     open OUT,">$d_inf1";

     for (0..$#orig_meme){
       if($orig_meme[$_]=~/^MEME version/){
         print OUT $orig_meme[$_],"\n";
       }elsif($orig_meme[$_]=~/^strands:/){
         print OUT $orig_meme[$_],"\n\n";
       }elsif($orig_meme[$_]=~/^MOTIF\s+(\d+)/){
         $_motif=$1;
         print OUT $orig_meme[$_],"\n" if ( ! defined $hash2{$_motif} );
       }elsif($orig_meme[$_]=~/^letter-probability matrix:.+ w=\s+(\d+)/){
            print OUT $orig_meme[$_],"\n" if ( ! defined $hash2{$_motif} );
            $width=$1;
            for ($i=1;$i<=$width;$i++){
              print OUT $orig_meme[$_+$i],"\n" if ( ! defined $hash2{$_motif} );
            }
       }
     }

     foreach $k (keys %hash2){
         print OUT $hash2{$k};
     }

}

sub getMotifID{

     $line=shift; @G=(); %temph=(); $arg="";
     @G=split(" ",$line); @G=uniq_array(\@G);
     for(@G){ unless(defined($temph{${[split("_",$_)]}[1]})){$arg.=${[split("_",$_)]}[1]." "; $temph{${[split("_",$_)]}[1]}=1; }}
     return $arg;

}

sub getMotifID2{

     $line=shift; $H_se=shift; $M_se=shift;
     %temph=(); $arg="";
     ($sfile=$line)=~s/\.align//; chomp $sfile; $sfile=~s/cluster\//cluster\/Selected_Align\//;
     $outH=$sfile.".HUM.tmp"; open OH, ">$outH";
     $outM=$sfile.".MOU.tmp"; open OM, ">$outM";
     
     open FH, "<$sfile";
     while(<FH>){
       chomp; $nl=$_;
       @myln=split("\t", $nl);
       if($myln[1] eq "HUM"){
          unless(($myln[2]>=${[split("_", $H_se)]}[1])||($myln[3]<=${[split("_", $H_se)]}[0])){
            print OH $nl,"\n";
            unless( defined($temph{$myln[0]}) ){ $arg.=$myln[0]." "; $temph{$myln[0]}=1; }
          }
       }elsif($myln[1] eq "MOU"){
          unless(($myln[2]>=${[split("_", $M_se)]}[1])||($myln[3]<=${[split("_", $M_se)]}[0])){
            print OM $nl,"\n";
            unless( defined($temph{$myln[0]}) ){ $arg.=$myln[0]." "; $temph{$myln[0]}=1; }
          }
       }
     }
     close FH;
     close OH;
     close OM;

     return $arg;

}
