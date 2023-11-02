#!/usr/bin/perl -w

use experimental qw(smartmatch);

($meme,@motif_num)=@ARGV;

open MEME, $meme;
while(<MEME>){
     chomp;
     if(/^MEME version/){
         print $_,"\n";
     }elsif(/^strands:/){
         print $_,"\n\n";
     }elsif(/^MOTIF\s+(\w+)/){
         $motif=$1;
         if($motif~~@motif_num){ print $_,"\n"; } 
     }elsif(/^letter-probability matrix:.+ w=\s+(\d+)/){
         if($motif~~@motif_num){
            print $_,"\n";
            $width=$1;
            for ($i=1;$i<=$width;$i++){ 
              $newLine=<MEME>;
              print $newLine;
            }
         }
     }
}
