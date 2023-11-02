#!/usr/bin/perl -w

BEGIN {
  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
}

use Lncrnaevo::Myownperlmodule qw(getmotif pri_motif motif_define motif_inf motif_n50_inf motif_n50_per);

$input_motif=shift;
$input_fasta=shift;

$ouput_motif=shift; open OUT, ">$ouput_motif";

$in_motif_hash = &getmotif($input_motif);
%in_motif_hash = %{$in_motif_hash};

foreach $key (sort {$a cmp $b} keys %in_motif_hash){
   %key_motif = %{$in_motif_hash{$key}};
   $reinf = &motif_define(\%key_motif, 0.75, $key, $input_motif, $input_fasta);
   %ou_motif_hash = %{$reinf};
   &pri_motif(\%ou_motif_hash, $ouput_motif);
}

close OUT;

