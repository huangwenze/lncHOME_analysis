#!/usr/bin/perl -w

#BEGIN {
#  unshift @INC,"/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
#}

#use Lncrnaevo::Myownperlmodule qw(fimo2Marray eucdist fimolines2blocks motif_include);
use Myownperlmodule qw(fimo2Marray eucdist fimolines2blocks motif_include);
#use experimental qw(smartmatch);

#require '/Share/home/zhangqf/usr/perl/lib/site_perl/5.22.0/Parallel/ForkManager.pm';
#use Parallel::ForkManager;
use List::Util qw(min max);

$InScoreRank = shift;
$_TrXfna = shift;
$_dir = shift;


#$mybindir="/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";

##read and store $InScoreRank (made by get_ML_preFeatures.pl), e.g.:
##/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/MEME_RESULTS_OLD1/DIR_NONHSAG015195_NONMMUG009300/DP_NONHSAG015195_NONMMUG009300.out	0011	2.123	0.01941	2	3
##FIR002_3_-_-9.928 FIR015_1_+_-9.947
##SEC000_3_-_-9.983 SEC003_1_+_-9.983
##FIR002_0 FIR010_307 FIR015_585
##SEC000_0 SEC001_43 SEC003_187
open IN, $InScoreRank;
while(<IN>){
     next if (/^#/);
     chomp;
     if (/^\//){
        $line1=$_; 
        push @TOT, $line1;
     }
}
close IN;

$compare_3rd_fimo=$_dir."/compare_3rd/compare_3rd_fimo.sort.out";

#@species = ("HUM","MOU","COW","OPS","CHC","LIZ","FRG","ZBF");
#@sp_arr  = ("human", "mouse", "cow", "opossum", "chicken", "lizard", "frog", "zebrafish");
#@fsp_arr = ("Homo sapiens", "Mus musculus", "Bos taurus", "Monodelphis domestica", "Gallus gallus", "Anolis carolinensis", "Xenopus tropicalis", "Danio rerio");
#for (0..$#species){ $SP{ $species[$_] }=$species[$_]; }

#$mp = Parallel::ForkManager->new(20);

##Processing $InScoreRank one by one
for $ele (@TOT){

#$mp->start and next;

     @_3rd_cluster=();

     ##specify dir of 3rd fimo file, and read it to a hash
     if ($ele=~/^\//){
      @F=split(/\t/,$ele);
      $DP_results=$F[0]; 
      $fimo_s=$DP_results;   $fimo_s=~s/\.align//;  $fimo_s=~s/cluster\//cluster\/Selected_Align\//; 

      ############$fimo_Hs=$fimo_s.".HUM.fimo.s"; $fimo_Ms=$fimo_s.".MOU.fimo.s";             ### ref.  getMotifID
      $fimo_Hs=$fimo_s.".HUM.tmp"; $fimo_Ms=$fimo_s.".MOU.tmp";                               ### ref.  getMotifID2

      $fileID=$DP_results;   $fileID=~s/.+\///;  $fileID=~s/_\S+\.align//;
      $fileID=$fileID.".".$F[1];
      $_TrXfna=~s/\S+\///; ($sp_abbr=$_TrXfna)=~s/\.fa//;
      $prefix=$sp_abbr.".".$fileID;

      $infimo3rd=$_dir."/compare_3rd/".$prefix.".fimo";
      if ( -s $infimo3rd ){
   
         $infimo3rds=$_dir."/compare_3rd/".$prefix.".fimo.sort";
         `cat $infimo3rd|awk -F '\t' '\$7<0.0001{print \$0}'|sort -k2,2 -k3,3n -k4,4n >$infimo3rds; rm $infimo3rd`;

         $prefixH=$prefix.".HUM"; $prefixM=$prefix.".MOU";

         if ( not $? ) {
         `perl motif_block_merge_edgeR_3rd.pl 10 200 $fimo_Hs $infimo3rds $prefixH`;
         }
         if ( not $? ) {
         `perl motif_block_merge_edgeR_3rd.pl 10 200 $fimo_Ms $infimo3rds $prefixM`;
         }

         if ( not $? ) {
		 $some_dir=$_dir."/compare_3rd/cluster"; $some_dir_file=$some_dir."/".$prefix."*.vs.*_3rd.sort";
		 if (glob("$some_dir_file")){
		 opendir($prefix, $some_dir);
		 @_3rd_cluster = grep { !/^\./ && /^$prefix\..+\.vs\..+\.sort$/ } readdir($prefix);
		 for ( @_3rd_cluster ){ s/^/$some_dir\//; }
		 closedir $prefix; 
                 }
         }

      }else{
        exit;
      }


    for (@_3rd_cluster){ 
#print "KOKO	", $_,"\n";
    if ( -s $_ ) {
#print "OKOK     ", $_,"\n";
     $_3rd_cluster_ele=$_;  #`rm -rf $_3rd_cluster_ele`;
 
     $postfix=$_3rd_cluster_ele; $postfix=~s/\S+vs\.(\S+)\.sort/$1/;
     $mymark=$_3rd_cluster_ele; $mymark=~s/\S+\.(\w+)\.vs\S+\.sort/$1/;
#print "MARK     ", $mymark,"\n";
     `mkdir -p $_dir/DP_EUC_3rd/$prefix`;
     $outputfile="${_dir}/DP_EUC_3rd/$prefix/$postfix"; #`rm -rf ${outputfile}*`;

      $H_motif_str=fimo2Marray("$fimo_Hs", "M"); @H_motif_str=@{$H_motif_str};
      $M_motif_str=fimo2Marray("$fimo_Ms", "M"); @M_motif_str=@{$M_motif_str};
      $TRD_motif_str=fimo2Marray("$_3rd_cluster_ele", "M"); @TRD_motif_str=@{$TRD_motif_str};

      if (scalar(@TRD_motif_str) >= 10){ 

if($mymark eq "HUM"){

       $H_share_motif=intersect_two_redun_array(\@H_motif_str,\@TRD_motif_str);

       if ($H_share_motif >=10){

         $outputfile_H=$outputfile.".HUM";

### ref.          `$mybindir/dynprogress4 $outputfile_H $fimo_Hs $_3rd_cluster_ele 0.01 -t 0.2`;

`cat $fimo_Hs $_3rd_cluster_ele >$outputfile_H` if ( not $? );
`dynprogress3 $outputfile_H 0.01 -t 0.2` if ( not $? );

if( not $? ){
`echo $_3rd_cluster_ele >>$compare_3rd_fimo`;
`cat $_3rd_cluster_ele >>$compare_3rd_fimo`;
`rm $_3rd_cluster_ele`;
}

         print "$outputfile_H\tMicroH_Finished\n" if ( not $? );

         $out_H_align=$outputfile_H."*"; $out_H_fimo=$outputfile_H.".HUM.fimo.s";
	 $filesizeH= -s "$out_H_fimo";
	 if ( $filesizeH<=1500 ){ 
         `rm -rf $out_H_align`; 
         }else{
         $out_H_Falign=$outputfile_H.".align";
         $H_score=getScore($out_H_Falign);
         if ( $H_score<10 ){ `rm -rf $out_H_align`; }
         }

       }     ### testing if shared motif number >3

}elsif($mymark eq "MOU"){

       $M_share_motif=intersect_two_redun_array(\@M_motif_str,\@TRD_motif_str);
       if ($M_share_motif >=10){

         $outputfile_M=$outputfile.".MOU";

### ref.          `$mybindir/dynprogress4 $outputfile_M $fimo_Ms $_3rd_cluster_ele 0.01 -t 0.2`;

`cat $fimo_Ms $_3rd_cluster_ele >$outputfile_M` if ( not $? );
`dynprogress3 $outputfile_M 0.01 -t 0.2` if ( not $? );

if( not $? ){
`echo $_3rd_cluster_ele >>$compare_3rd_fimo`;
`cat $_3rd_cluster_ele >>$compare_3rd_fimo`;
`rm $_3rd_cluster_ele`;
}

         print "$outputfile_M\tMicroH_Finished\n" if ( not $? );

         $out_M_align=$outputfile_M."*"; $out_M_fimo=$outputfile_M.".MOU.fimo.s";
         $filesizeM= -s "$out_M_fimo";
         if ( $filesizeM<=1500 ){ 
         `rm -rf $out_M_align`; 
         }else{
         $out_M_Falign=$outputfile_M.".align";
         $M_score=getScore($out_M_Falign);
         if ( $M_score<10 ){ `rm -rf $out_M_align`; }
         }

       }     ### testing if shared motif number >3

}

      }     ### testing if 3rd fimo file have lines >3

    }
   } 

   }

#$mp->finish;
}

#$mp->wait_all_children;



sub getScore{
    $myFile=shift; $my_score=0;
    open MFH, "<$myFile";
    while(<MFH>){
       chomp; $nl=$_;
       if($nl=~/\//){
       @my_ln=split("\t", $nl);
       $my_score=$my_ln[2];
       }
    }
    close MFH;
    return $my_score;
}

sub intersect_two_redun_array{
    $A=shift; $B=shift; $_common=0;
    %counta=%countb=();
    for $_a_ele (@{$A}){$counta{$_a_ele}++;}
    for $_b_ele (@{$B}){$countb{$_b_ele}++;}
    foreach $_ele_ (keys %counta){
        if(defined($countb{$_ele_})){ $_common+=min($counta{$_ele_},$countb{$_ele_}); }
    }
    return $_common;
}

