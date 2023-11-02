#!/usr/bin/perl -w

#use Env qw(LD_LIBRARY_PATH);

#BEGIN {
#  $LD_LIBRARY_PATH .= '/bin:/usr/local/gcc-4.9.2/lib64:/opt/openmpi-1.8.7/lib:/opt/intel-xe_2015u3/lib/intel64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOMEZ/usr/openmpi/lib:/path/to/bcl-3.1.0-Linux-x86_64:$HOMEZ/usr/readline-6.3/lib:$HOME/tools/binutils/lib:$HOMEZ/usr/lib64:$HOMEZ/tmp/expat/lib:$HOMEZ/usr/meme/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$HOME/tools/stamp/gsl_dir/gsl/lib:$HOME/tools/binutils/lib:$HOME/usr/lib64:/opt/lsf/9.1/linux2.6-glibc2.3-x86_64/lib:/opt/intel-xe_2015u3/composer_xe_2015.3.187/compiler/lib/intel64:$HOME/tmp/expat/lib:$HOME/tools/meme/lib:$HOMEZ/usr/openmpi/lib:$HOMEZ/shaodi/app/samtools/zlib-1.2.8/build/lib:$HOMEZ/shaodi/app/samtools/ncurses-6.0/build/lib:$HOMEZ/jellyfish/lib:$LD_LIBRARY_PATH';
#  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
#}

use Myownperlmodule qw(ModuleFinder_L2L_preprocess ModuleFinder_L2L ModuleFinder_L2S_HMpreprocess ModuleFinder_L2S_HM DifferentialDistributionAnalysis uniq_array eucdist intersect_two_redun_array);
use List::Util qw(min max);
use experimental qw(smartmatch);

$step=shift; $wind=shift;
$fimo_redundant_motif=shift;     ### full directory

$mydir=basedir($fimo_redundant_motif);
`mkdir -p $mydir/cluster`;


open IN, $fimo_redundant_motif;
# 0224    FIR 1   10  +   10.7241 9.93e-05        GGGCCCAGGG
#HM06407 HUM     50      65      -       15.2667 3.98e-06                CCCCGCTTCCAGCACT
#HM05498 HUM     50      68      +       14.1193 6.38e-06                AGTGCTGGAAGCGGGGATA
while(<IN>){
     next if (/^#/);
     chomp;
     @Ln=split("\t", $_);
     push @{$hash{$Ln[1]}}, $Ln[0]."_".$Ln[4];
     push @{$hash_FL{$Ln[1]}}, join("\t", @Ln);
     push @{$hash_pos{$Ln[1]}}, $Ln[2]."_".$Ln[3];

     push @{$hash_region{$Ln[1]}{"st"}}, $Ln[2];
     push @{$hash_region{$Ln[1]}{"ed"}}, $Ln[3];
}
close IN;

$hum_s=$hum_e=$mou_s=$mou_e=0;
foreach $mykey (keys %hash_region){
     if($mykey eq "HUM"){
     $hum_s=min(@{$hash_region{$mykey}{"st"}});
     $hum_e=max(@{$hash_region{$mykey}{"ed"}});
     $hum_length=$hum_e-$hum_s+1;
     }elsif($mykey eq "MOU"){
     $mou_s=min(@{$hash_region{$mykey}{"st"}});
     $mou_e=max(@{$hash_region{$mykey}{"ed"}});
     $mou_length=$mou_e-$mou_s+1;
     }
}

print $hum_length,"\t",$mou_length,"\n";

$thl=600;
if(($hum_length<$thl)&&($mou_length<$thl)){
     $mode=0;
}elsif(($hum_length>=$thl)&&($mou_length>=$thl)){
     $mode=1;
}elsif(($hum_length>=$thl)&&($mou_length<$thl)){
     $mode=2;
}elsif(($hum_length<$thl)&&($mou_length>=$thl)){
     $mode=3;
}
print $mode,"\n";

if($mode==0){

     $newfimo="C0_HUM_".$hum_s."_".$hum_e."_MOU_".$mou_s."_".$mou_e.".sort";
     `cp $fimo_redundant_motif $mydir/cluster/$newfimo`;

}elsif($mode==1){

     $fimo_final_filter=shift;
     ( $len_hash, $pool_hash, $poolfimo_hash )=ModuleFinder_L2L_preprocess( $step, $wind, \%hash, \%hash_pos, \%hash_FL );
     %mylength=%{$len_hash}; 
     %poolH=%{$pool_hash};
     %poolfimoH=%{$poolfimo_hash};
     ( $seqa, $seqb, $outfinalfile, $_hash_arr_seq1, $_hash_arr_seq2 )=ModuleFinder_L2L( $mydir, $fimo_final_filter, \%poolH, \%poolfimoH );
     %my_hash_arr_seq1=%{$_hash_arr_seq1};
     %my_hash_arr_seq2=%{$_hash_arr_seq2};

     $mydbscanRes=dbscan($outfinalfile, 15);
     $fullcluster=GoupCluster($mydbscanRes) if ( -s $mydbscanRes );
     &MotifRetreiving($fullcluster, $seqa, $seqb, \%my_hash_arr_seq1, \%my_hash_arr_seq2, \%mylength) if (defined $fullcluster);

}elsif($mode>1){

     if ($mode==2){
         $tmpsp="MOU";
     }elsif($mode==3){
         $tmpsp="HUM";
     }
          
     $fimo_final_filter=shift;
     ( $len_hash, $pool_hash, $poolfimo_hash )=ModuleFinder_L2S_HMpreprocess( $tmpsp, $step, $wind, \%hash, \%hash_pos, \%hash_FL );
     %mylength=%{$len_hash};
     %poolH=%{$pool_hash};
     %poolfimoH=%{$poolfimo_hash};
     ( $seqa, $seqb, $outfinalfile, $_hash_arr_seq1, $_hash_arr_seq2 )=ModuleFinder_L2S_HM( $mydir, $fimo_final_filter, $tmpsp, \%poolH, \%poolfimoH );
     %my_hash_arr_seq1=%{$_hash_arr_seq1};
     %my_hash_arr_seq2=%{$_hash_arr_seq2};

     $mydbscanRes=dbscan($outfinalfile, 5);
     $fullcluster=GoupCluster($mydbscanRes) if ( -s $mydbscanRes );
     &MotifRetreiving($fullcluster, $seqa, $seqb, \%my_hash_arr_seq1, \%my_hash_arr_seq2, \%mylength) if (defined $fullcluster);

}




####################
### subfunctions ###
####################

#sub "DBSCAN:Density Based Clustering of Applications with Noise": redundancy_removing_from_"file fimo_final_filter", collecting pair-regions with common motifs, and clustering them
sub dbscan{
   $outfinal=shift; $outfilesize= -s "$outfinal";
   $parameter=shift;
   $fimo_final_filter_dbscan=$outfinal.".dbscan";
   if ( $outfilesize>1500 ){    ### about 10 lines
   `python /home/huangwenze/lncHOME/05MicroH/new_code/exhdbscan.py $outfinal $parameter $fimo_final_filter_dbscan`;
   return  $fimo_final_filter_dbscan if ( not $? );
   }
}

#sub retreive motifs from pair-regions_above without noise
sub GoupCluster{
    $dbscanInput=shift;
    open INF, "$dbscanInput";
    while(<INF>){
      next if (/^#/);
      chomp;
      @NLn=split("\t", $_);
      if ($NLn[2]>0){
         push @{$cluster{$NLn[2]}}, $NLn[0]."_".$NLn[1];
      }
    }
    close INF;   
    return \%cluster;
}

#sub 
sub MotifRetreiving{
    $centermotifs=shift;
    $seq1=shift;   #HUM
    $seq2=shift;
    $t_hash_arr_seq1=shift; %hash_arr_seq1=%{$t_hash_arr_seq1};   #MOU
    $t_hash_arr_seq2=shift; %hash_arr_seq2=%{$t_hash_arr_seq2};
    $t_length=shift; %length=%{$t_length};
    foreach $_key (keys %{$centermotifs}){
       for (@{$$centermotifs{$_key}}){
           push @{$cluster_hash_FIR_common_motifs{$_key}}, @{$hash_arr_seq1{$_}};
           push @{$cluster_hash_SEC_common_motifs{$_key}}, @{$hash_arr_seq2{$_}};
           push @{$cluster_hash_FIR_ord{$_key}}, ${[split("_", $_)]}[0]; ###s=min;e=max+120
           push @{$cluster_hash_SEC_ord{$_key}}, ${[split("_", $_)]}[1];
       }       
    }

    foreach $x (sort {$a<=>$b} keys %cluster_hash_FIR_ord){

           $re= "C".$x. "_HUM_". min(@{$cluster_hash_FIR_ord{$x}}).
                            "_". min( max(@{$cluster_hash_FIR_ord{$x}})+300, $length{$seq1}+30 ).  # +120
                        "_MOU_". min(@{$cluster_hash_SEC_ord{$x}}).
                            "_". min( max(@{$cluster_hash_SEC_ord{$x}})+300, $length{$seq2}+30 );  # +120

           $outf=$mydir."/cluster/".$re; $outfs=$mydir."/cluster/".$re.".sort";
           open $re, ">".$outf;
           print $re join("\n", @{[uniq_array(\@{$cluster_hash_FIR_common_motifs{$x}})]}),"\n",
                     join("\n", @{[uniq_array(\@{$cluster_hash_SEC_common_motifs{$x}})]}),"\n";
           close $re;

           `sort -k2,2 -k3,3n -k4,4n $outf >$outfs; rm -rf $outf`;
           
    }
}

#sub
sub basedir{
    $basedir=shift;
    if($basedir=~/^\S+\//){
    $basedir=~s/(^\S+)\/\S+/$1/g;
    }else{
    $basedir='./';
    }
    return $basedir;
}
