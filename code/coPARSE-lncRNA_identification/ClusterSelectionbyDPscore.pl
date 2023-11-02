#!/usr/bin/perl -w

$total_align=shift;
$threshold=shift;
$outfile=shift;



open IN, $total_align;
while(<IN>){
     next if (/^#/);
     chomp;
     if (/^\//){
        $line=$_;
        @ln=split("\t", $line);
        $cluster=$ln[0]; $cluster=~s/\S+\/(C\d+)_\S+/$1/;
        $dp_score=$ln[2];
        $selected{$cluster}=$dp_score;

        #$euc_score=$ln[3];
        #$fir_sted=$ln[4];
        #$sec_sted=$ln[5];
        #$fir_clu_id=$cluster."FIR"; 
        #while (length($fir_clu_id)!=8){
        #    $fir_clu_id="0".$fir_clu_id;
        #}
        #$sec_clu_id=$cluster."SEC";
        #while (length($sec_clu_id)!=8){
        #    $sec_clu_id="0".$sec_clu_id;
        #}
        #$match_score{$fir_clu_id."_".$sec_clu_id}=$ln[2]; 
        #$fir_sted{$fir_clu_id}=${[split("_",$ln[4])]}[0];
        #$sec_sted{$sec_clu_id}=${[split("_",$ln[5])]}[0];

     }
}
close IN;

open OUT, ">$outfile";
foreach $key (keys %selected){
     if ( $selected{$key} >= $threshold ){
        print OUT $key,"\n";
        #print $key,"\n";
     }
}
