#!/usr/bin/perl

package Myownperlmodule;

use strict;
use warnings;
use experimental qw(smartmatch);
use List::Util qw(min max);
use Math::Trig;
#use Math::BigFloat;
#use bignum;
#use POSIX 'log10';

use Exporter qw(import);
our @EXPORT_OK = qw(ModuleFinder_L2L ModuleFinder_L2S_HM ModuleFinder_L2S ModuleFinder_L2L_preprocess ModuleFinder_L2S_HMpreprocess ModuleFinder_L2S_preprocess Marray2fimo fimo2Marray INF grep_meme_recognized_by_mcast getmotif pri_motif motif_define motif_inf motif_n50_inf motif_n50_per DifferentialDistributionAnalysis uniq_array intersect_two_redun_array readfa eucdist dynprog_matrix_pair_fimo dynprog_matrix_selected_and_fimo matrix_norm macro_dynprog fimolines2blocks motif_include motifpos getArrIdx);

###Module List###
# getmotif 
# pri_motif 
# motif_define 
# motif_inf 
# motif_n50_inf 
# motif_n50_per 
# DifferentialDistributionAnalysis
# uniq_array				in: a reference of array;			out: a uniq array
# intersect_two_redun_array		in: two reference of arrays;			out: a reference array of common element, and its number
# readfa				in: fasta file (absolute path);			out: a reference of hash
# eucdist				in: two reference of arrays;			out: a numeric variable
# dynprog_matrix_pair_fimo		in: two reference of arrays (format 1);		out: a numeric variable of "MATCH" score 
# dynprog_matrix_selected_and_fimo	in: two reference of arrays (format 2);		out: a numeric variable of "MATCH" score
# matrix_norm				in: a reference of 2-D hash of matrix;		out: a reference of 2-D hash of normalized matrix
# macro_dynprog				in: two variables of 6bit-string, and a reference hash of scoring matrix, and a reference hash of motif groups;
#                                       out: a variable containing formatted text (dynprog hits), and an array of dynamic programmming matrix
# fimolines2blocks			in: a value(0, 1, 2), and a refernce hash;	out: two or three or four(input an additional "1" or "2") reference of hashs
# motif_include				in: a value(0~1), an array and a refernce hash;	out: 1 or 0
#
#
#
#################
# NOTE: a. defining a variable by "my $v", if $v is not given a value after that, please do not print it (uninitial value)! contrarily, @variable and %variable will be set to "NULL". 
#       b.
#################



sub ModuleFinder_L2L{
    
    ### input: $mydir; $outprefix; %pool; %pool_fimo;
    ### ouput: $dbscan-file; %hash_arr_seq1 or/and 2;

    my $mydir=shift;
    my $outprefix=shift; #my $FINALOUT=$outprefix;
    my $_pool=shift; my %pool=%{$_pool};
    my $_pool_fimo=shift; my %pool_fimo=%{$_pool_fimo};

    my $outfinal=$mydir."/".$outprefix;
    open FINALOUT, ">>$outfinal";
    my (%hash_arr_seq1, %hash_arr_seq2);

    my ($seq1, $seq2)=sort keys %pool;
    my %SEQ1=%{$pool{$seq1}}; my %SEQ2=%{$pool{$seq2}};
    for my $k1 (sort {$a <=> $b} keys %SEQ1){
    my @temp1=@{$SEQ1{$k1}}; my $len_temp1=scalar @temp1; 
    my (@k1, @k2, @temp1NUM, @temp2NUM, @temp_commotif_arr_NUM);
    for my $k2 (sort {$a <=> $b} keys %SEQ2){
      my @temp2=@{$SEQ2{$k2}}; my $len_temp2=scalar @temp2; 
      push @k1, $k1; push @temp1NUM, $len_temp1;
      push @k2, $k2; push @temp2NUM, $len_temp2;
      my ($intersect_arr, $intersect_num)=intersect_two_redun_array(\@temp1, \@temp2);
      my @temp_commotif_arr=@{$intersect_arr};
      push @temp_commotif_arr_NUM, $intersect_num;
      for (@temp_commotif_arr){
           push @{$hash_arr_seq1{$k1."_".$k2}}, @{$pool_fimo{$seq1}{$k1}{$_}};
           push @{$hash_arr_seq2{$k1."_".$k2}}, @{$pool_fimo{$seq2}{$k2}{$_}};
      }
    }
    my @temp1pv=DifferentialDistributionAnalysis(\@temp1NUM,\@temp_commotif_arr_NUM);
    my @temp2pv=DifferentialDistributionAnalysis(\@temp2NUM,\@temp_commotif_arr_NUM);
    for(0..$#temp1pv){
      if((($temp1pv[$_]*$temp2pv[$_])<=0.00001)&&($temp_commotif_arr_NUM[$_]>=5)){
          print FINALOUT "$k1[$_]\t$temp1NUM[$_]\t$k2[$_]\t$temp2NUM[$_]\t$temp_commotif_arr_NUM[$_]\t$temp1pv[$_]\t$temp2pv[$_]\t",$temp1pv[$_]*$temp2pv[$_],"\n";
      }
    } 
    }
    close FINALOUT;

    return ($seq1, $seq2, $outfinal, \%hash_arr_seq1, \%hash_arr_seq2);
}

#
sub ModuleFinder_L2S_HM{

    ### input: $mydir; $outprefix; %pool; %pool_fimo;
    ### ouput: $dbscan-file; %hash_arr_seq1 or/and 2;
    
    my $mydir=shift;
    my $outprefix=shift; #my $FINALOUT=$outprefix;
    my $_sp=shift;
    my $_pool=shift; my %pool=%{$_pool};
    my $_pool_fimo=shift; my %pool_fimo=%{$_pool_fimo};

    my $outfinal=$mydir."/".$outprefix;
    open FINALOUT, ">>$outfinal";
    my (%hash_arr_seq1, %hash_arr_seq2);

    my ($_seq1, $_seq2)=sort keys %pool;
    my $seq1=($_seq1 eq $_sp)?$_sp:$_seq2; 
    my $seq2=($_seq2 eq $_sp)?$_seq1:$_seq2;
    my %SEQ1=%{$pool{$seq1}}; my %SEQ2=%{$pool{$seq2}};
    for my $k1 (sort {$a <=> $b} keys %SEQ1){
    my @temp1=@{$SEQ1{$k1}}; my $len_temp1=scalar @temp1;
    #my (@k1, @k2, @temp1NUM, @temp2NUM, @temp_commotif_arr_NUM);
    my (@k1, @k2, @temp2NUM, @temp_commotif_arr_NUM);
    for my $k2 (sort {$a <=> $b} keys %SEQ2){
      my @temp2=@{$SEQ2{$k2}}; my $len_temp2=scalar @temp2;
      push @k1, $k1; #push @temp1NUM, $len_temp1;
      push @k2, $k2; push @temp2NUM, $len_temp2;
      my ($intersect_arr, $intersect_num)=intersect_two_redun_array(\@temp1, \@temp2);
      my @temp_commotif_arr=@{$intersect_arr};
      push @temp_commotif_arr_NUM, $intersect_num;
      for (@temp_commotif_arr){
           if($seq1 eq "HUM"){
           push @{$hash_arr_seq1{$k1."_".$k2}}, @{$pool_fimo{$seq1}{$k1}{$_}};
           push @{$hash_arr_seq2{$k1."_".$k2}}, @{$pool_fimo{$seq2}{$k2}{$_}};
           }else{
           push @{$hash_arr_seq1{$k2."_".$k1}}, @{$pool_fimo{$seq1}{$k1}{$_}};
           push @{$hash_arr_seq2{$k2."_".$k1}}, @{$pool_fimo{$seq2}{$k2}{$_}};
           }
      }
    }
    #my @temp1pv=DifferentialDistributionAnalysis(\@temp1NUM,\@temp_commotif_arr_NUM);
    my @temp2pv=DifferentialDistributionAnalysis(\@temp2NUM,\@temp_commotif_arr_NUM);
    for(0..$#temp2pv){
      #if((($temp1pv[$_]*$temp2pv[$_])<=0.00001)&&($temp_commotif_arr_NUM[$_]>=5)){
      #    print $FINALOUT "$k1[$_]\t$temp1NUM[$_]\t$k2[$_]\t$temp2NUM[$_]\t$temp_commotif_arr_NUM[$_]\t$temp1pv[$_]\t$temp2pv[$_]\t",$temp1pv[$_]*$temp2pv[$_],"\n";
      #}
      if((($temp2pv[$_])<=0.1)&&($temp_commotif_arr_NUM[$_]>=2)){
          if($seq1 eq "HUM"){
          print FINALOUT "$k1\t$len_temp1\t$k2[$_]\t$temp2NUM[$_]\t$temp_commotif_arr_NUM[$_]\t","n.a.","\t$temp2pv[$_]\t","n.a.","\n";
          }else{
          print FINALOUT "$k2[$_]\t$temp2NUM[$_]\t$k1\t$len_temp1\t$temp_commotif_arr_NUM[$_]\t","n.a.","\t$temp2pv[$_]\t","n.a.","\n";
          }
      }
    }
    }
    close FINALOUT;

    if($seq1 eq "HUM"){
    return ($seq1, $seq2, $outfinal, \%hash_arr_seq1, \%hash_arr_seq2);
    }else{
    return ($seq2, $seq1, $outfinal, \%hash_arr_seq2, \%hash_arr_seq1);
    }
}

#
sub ModuleFinder_L2S{

    ### input: $mydir; $outprefix; $ele_3rd; $_sp; %ppool; %pool; %pool_fimo;
    ### ouput: $dbscan-file; %hash_arr_seq1 or/and 2

    my $mydir=shift;
    my $outprefix=shift;
    my $_sp=shift; my $ele=shift; my $seq2_mark=shift;
    my $_ppool=shift; my %ppool=%{$_ppool};
    my $_pool=shift; my %pool=%{$_pool};
    my $_pool_fimo=shift; my %pool_fimo=%{$_pool_fimo};

    my $outfinal=$mydir."/".$outprefix;
    open FINALOUT, ">>$outfinal";
    my %hash_arr_seq;

    my $seq1=$_sp; my $seq2=$ele; 
#for my $x (keys %ppool){
#    for my $y (keys %{$ppool{$x}}){
#      print $x, "\t:\t", $y, "\t:\t", join(" ", @{$ppool{$x}{$y}}),"\n";
#    }
#}
    my %SEQ1=%{$ppool{$seq1}}; 
    my %SEQ2=%{$pool{$seq2}};
    for my $k1 (sort {$a <=> $b} keys %SEQ1){
    my @temp1=@{$SEQ1{$k1}}; my $len_temp1=scalar @temp1; 
    my (@k2, @temp2NUM, @temp_commotif_arr_NUM); 
    for my $k2 (sort {$a <=> $b} keys %SEQ2){
      my @temp2=@{$SEQ2{$k2}}; my $len_temp2=scalar @temp2; 
      push @k2, $k2; 
      push @temp2NUM, $len_temp2;
      my ($intersect_arr, $intersect_num)=intersect_two_redun_array(\@temp1, \@temp2);
      my @temp_commotif_arr=@{$intersect_arr};
      push @temp_commotif_arr_NUM, $intersect_num;
      for (@temp_commotif_arr){
           my @tmp_pool=@{[map {s/TCONS\S+/$seq2_mark/; $_} @{$pool_fimo{$seq2}{$k2}{$_}}]};
           push @{$hash_arr_seq{$k1."_".$k2}}, @tmp_pool;
      }
    }
    my @temp2pv=DifferentialDistributionAnalysis(\@temp2NUM,\@temp_commotif_arr_NUM);
    for(0..$#temp2pv){
      if((($temp2pv[$_])<=0.1)&&($temp_commotif_arr_NUM[$_]>=2)){
          print FINALOUT "$k1\t$len_temp1\t$k2[$_]\t$temp2NUM[$_]\t$temp_commotif_arr_NUM[$_]\t","n.a.","\t$temp2pv[$_]\t","n.a.","\n";
      }
    }  
    }
    close FINALOUT;

    return ($outfinal, \%hash_arr_seq);
}

#
sub ModuleFinder_L2L_preprocess{

    ### input: $step; $with; %hash; %hash_pos; %hash_FL
    ### ouput: %length; %pool, %pool_fimo

    my $step=shift; my $wind=shift;
    my $hash=shift; my %hash=%{$hash};
    my $hash_pos=shift; my %hash_pos=%{$hash_pos};
    my $hash_FL=shift; my %hash_FL=%{$hash_FL};

    my (%length, %motifnum, %pool, %pool_fimo);

    foreach my $key (keys %hash_pos){
    my @temp=@{$hash_pos{$key}};
    $length{$key}=${[split("_", $temp[-1])]}[0]-${[split("_", $temp[0])]}[0]+30;
    $motifnum{$key}=scalar(@temp)-1;
    }
    my $i;
    my $j=$step;
    my $k=$wind;
    foreach my $nkey (sort keys %hash_pos){
      for ($i=0; $i<=$length{$nkey}; $i+=$j){
        for (0..$motifnum{$nkey}){
          last if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] > ($i+$k));
          if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] >= $i){
            push @{$pool{$nkey}{$i}}, ${$hash{$nkey}}[$_];
            push @{$pool_fimo{$nkey}{$i}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
          }
        }
      }
    }

    return (\%length, \%pool, \%pool_fimo);
}

#
sub ModuleFinder_L2S_HMpreprocess{
    
    ### input: $_sp; $step; $with; %hash; %hash_pos; %hash_FL
    ### ouput: %length; %pool, %pool_fimo

    my $_sp=shift;
    my $step=shift; my $wind=shift;
    my $hash=shift; my %hash=%{$hash};
    my $hash_pos=shift; my %hash_pos=%{$hash_pos};
    my $hash_FL=shift; my %hash_FL=%{$hash_FL};

    my (%length, %motifnum, %pool, %pool_fimo);

    foreach my $key (keys %hash_pos){
    my @temp=@{$hash_pos{$key}};
    $length{$key}=${[split("_", $temp[-1])]}[0]-${[split("_", $temp[0])]}[0]+30;
    $motifnum{$key}=scalar(@temp)-1;
    }
    my ($i, $j, $k);
    foreach my $nkey (sort keys %hash_pos){
      $j=($nkey eq $_sp)?($length{$nkey}+1):$step;
      $k=($nkey eq $_sp)?($length{$nkey}):$wind;
      for ($i=0; $i<=$length{$nkey}; $i+=$j){
        for (0..$motifnum{$nkey}){
          last if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] > ($i+$k));
          if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] >= $i){
            push @{$pool{$nkey}{$i}}, ${$hash{$nkey}}[$_];
            push @{$pool_fimo{$nkey}{$i}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
          }
        }
      }
    }

    return (\%length, \%pool, \%pool_fimo);
}

#
sub ModuleFinder_L2S_preprocess{

    ### input:  $step; $with; %phash; %phash_pos; %hash; %hash_po; %hash_FL
    ### ouput:  %length; %ppool; %pool, %pool_fimo

    my $thl=shift;
    my $step=shift; my $wind=shift;
    my $phash=shift; my %phash=%{$phash};
    #my $phash_pos=shift; my %phash_pos=%{$phash_pos};  
    my $hash=shift; my %hash=%{$hash};
    my $hash_pos=shift; my %hash_pos=%{$hash_pos};
    my $hash_FL=shift; my %hash_FL=%{$hash_FL};

    my (%length, %motifnum, %ppool, %pool, %pool_fimo);

    my $pnkey=${[keys %phash]}[0];
    @{$ppool{$pnkey}{"0"}}=@{$phash{$pnkey}};

    foreach my $key (keys %hash_pos){
      my @temp=@{$hash_pos{$key}};
      $length{$key}=${[split("_", $temp[-1])]}[0]-${[split("_", $temp[0])]}[0]+30;
      $motifnum{$key}=scalar(@temp)-1;
    }

    my $i;
    my $j=$step;
    my $k=$wind;
    foreach my $nkey (sort keys %hash_pos){
      if($length{$nkey} <$thl){
             @{$pool{$nkey}{"0"}}=@{$hash{$nkey}};
             @{$pool_fimo{$nkey}{"0"}}=@{$hash_FL{$nkey}};
      }else{
      for ($i=0; $i<=$length{$nkey}; $i+=$j){
        for (0..$motifnum{$nkey}){
          last if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] > ($i+$k));
          if (${[split("_", ${$hash_pos{$nkey}}[$_])]}[0] >= $i){
             push @{$pool{$nkey}{$i}}, ${$hash{$nkey}}[$_];
             push @{$pool_fimo{$nkey}{$i}{${$hash{$nkey}}[$_]}}, ${$hash_FL{$nkey}}[$_];
          }
        }
      }
      }
    }

   return (\%length, \%ppool, \%pool, \%pool_fimo);
}

#
sub Marray2fimo{
     my $arr=shift; my @arr=@{$arr};
     my $fimoOutFile=shift;
     open OUT, ">$fimoOutFile";
     print OUT join("\n", @arr),"\n";
     close OUT;

     return $fimoOutFile;
}

#
sub fimo2Marray{
     my $fimo=shift; my @arr=();
     my $sig=shift;
     open FIMO, $fimo;
     while(<FIMO>){  
       next if (/^#/);
       chomp;
       my $myln=$_;
       if ($sig=~/M/){
       my @Ln=split(/\t/,$myln);  
       my $temp_motif_str=$Ln[0]."_".$Ln[4];
       push @arr, $temp_motif_str;  
       }elsif($sig=~/F/){
       push @arr, $myln;
       }
     }
     close FIMO;

     return \@arr;
}

#
sub INF{
    my $A=shift; my $C=shift; my $G=shift; my $T=shift;
    my $Ainf=($A==0)?0:$A*log($A)/log(2); my $Cinf=($C==0)?0:$C*log($C)/log(2);
    my $Ginf=($G==0)?0:$G*log($G)/log(2); my $Tinf=($T==0)?0:$T*log($T)/log(2);
    my $H= 2 + ($Ainf+ $Cinf + $Ginf + $Tinf);
    return $H;
}

#  
sub grep_meme_recognized_by_mcast{

my ($i_meme, $o_meme, @motif_num)=@_;
my ($motif, $width, $i, $newLine);

open MEME, $i_meme;
open OUT, ">$o_meme";
while(<MEME>){
     chomp;
     if(/^MEME version/){
         print OUT $_,"\n";
     }elsif(/^strands:/){
         print OUT $_,"\n\n";
     }elsif(/^MOTIF\s+(\w+)/){
         $motif=$1;
         if($motif~~@motif_num){ print OUT $_,"\n"; } 
     }elsif(/^letter-probability matrix:.+ w=\s+(\d+)/){
         if($motif~~@motif_num){
            print OUT $_,"\n";
            $width=$1;
            for ($i=1;$i<=$width;$i++){ 
              $newLine=<MEME>;
              print OUT $newLine;
            }
         }
     }
}

}

#  $motif = &getmotif($motif_db)
sub getmotif{
	my $file = shift;
	my $sen = ""; my %mot = (); my @sen = ();
	my $motif = ""; my $width = 0; my $newLine; my $mver; my $hmot;
	my $i = 0; my $j = 0;
	open(FILE1, $file)||die("open $file error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		if($sen =~m/^MEME version/){
			$mver = $sen."\n";
		}elsif($sen =~m/^strands:/){
			$mver = $mver.$sen."\n\n";
		}
		elsif($sen =~m/^MOTIF\s+(\w+)/){
			$hmot = $sen."\n";;
			$motif = $1;
		}elsif($sen =~m/^letter-probability matrix:.+ w=\s+(\d+)/){
			$width = $1;
			$hmot = $hmot.$sen;
			${$mot{$motif}}{0} = [$mver, $hmot, $width, 0, 0, 0, 0];
			for ($i=1; $i<=$width; $i++){ 
				$newLine=<FILE1>;
				chomp($newLine);
				@sen = split(/\s+/,$newLine);
				${$mot{$motif}}{$i} = [$sen[1],$sen[2],$sen[3],$sen[4]]; 
			}
		}
	}
	close FILE1;
	return (\%mot);
}

#  &pri_motif(\%motif_matrix, $motifID, "hum");
sub pri_motif{
	my $mot = shift; my %mot_inf = %{$mot}; 
        #my $motif_name = shift; my $spe_name = shift;
        #my $mot_file = $motif_name."_".$spe_name;
        my $mot_file = shift;
	open(TUO,">>",$mot_file);
	my $len = ${$mot_inf{0}}[2]; my $i; my $new_var = sprintf("%.3f", ${$mot_inf{0}}[6]);
	print TUO ${$mot_inf{0}}[0], ${$mot_inf{0}}[1]; #" renew= ",${$mot_inf{0}}[3],"/",${$mot_inf{0}}[4],"/",${$mot_inf{0}}[5]," inf= ",$new_var;
	#for($i=1; $i<=$len; $i++){
	#	$new_var = sprintf("#%.3f", &INF(${$mot_inf{$i}}[0], ${$mot_inf{$i}}[1], ${$mot_inf{$i}}[2], ${$mot_inf{$i}}[3]));
	#	print OUT $new_var;
	#}
	print TUO "\n";
	for($i=1; $i<=$len; $i++){
		$new_var = sprintf(" %.6f  %.6f  %.6f  %.6f \n", ${$mot_inf{$i}}[0], ${$mot_inf{$i}}[1], ${$mot_inf{$i}}[2], ${$mot_inf{$i}}[3]);
		print TUO $new_var;
	}
	close TUO;
}

#  $reinf = &motif_define(\%motif_matrix, 0.5, $mot_file, $fafile_to_scan);
sub motif_define{
	my $mot_mat = shift; my $p = shift; my $my_motif_id = shift; my $mo_file = shift; my $fa_file = shift; my $limit = 0.3; my $limit2 = 0.15;
	my %mat = %{$mot_mat}; my %remat = (); my %mot_seq = (); my %ori_mat = %mat; my %max_mat = %mat;
        my $sel_meme_out = $mo_file.".sel.meme";
        &grep_meme_recognized_by_mcast($mo_file, $sel_meme_out, $my_motif_id);
	my $fimo_out = $mo_file.".fimo";
	my $fimo_filter_fa = $fimo_out.".fa";
	my $fimo_filter = $fimo_out.".filter";
	my $sen = ""; my $seq = ""; my $sna = "";
	my @sen = (); my @sen1 = (); my @sen2 = ();
	my $num; my $fi_num; my $keys; my $w1; my $w2; my $sta; my $rem; my $max; my $ke; my $len = ${$mat{0}}[2];
	my $i = 0; my $j = 0; my $k = 0; 
	srand(time());
	foreach $keys (sort {$a<=>$b} keys %mat){
		$remat{$keys} = $mat{$keys};
	}
        `/Share/home/zhangqf2/tools/meme/bin/fimo --verbosity 1 --text $sel_meme_out $fa_file > $fimo_out; rm $sel_meme_out`;
	print "fimo finished!\n" if ( not $? );
	open(FILE1, $fimo_out)||die("open $fimo_out error!\n");
	open(OUT, ">", $fimo_filter_fa);
	$sen = <FILE1>;

	while($sen = <FILE1>){
		$i = $i + 1;
		chomp($sen);
		@sen = split(/\t/,$sen);
		$mot_seq{$i} = [$sen[5], $sen[6], $sen[8], $sen];
		print OUT ">",$i,"\n",$sen[8],"\n";
	}
	$fi_num = $i;
	close FILE1;
	close OUT;
        `/Share/home/zhangqf2/tools/dustmasker -in $fimo_filter_fa -out $fimo_filter -window 30 -level 10 -outfmt fasta`;
	print "dustmasker finished!\n" if ( not $? );
	$num = 0;
	open(FILE1, $fimo_filter)||die("open $fimo_filter error!\n");
	while($sen = <FILE1>){
		chomp($sen);
		if($sen =~m/>([0-9]+)$/){
			$i = $1;
			$sen = <FILE1>;
			chomp($sen);
			my $count = ($sen =~ s/([ATCG])/$1/g);
			my $lsen = length($sen);
			if($count/$lsen >= 0.6){
				$num = $num + 1;
			}else{
				delete $mot_seq{$i};
			}
		}
	}
	close FILE1;
        `rm $fimo_filter_fa`;
        `rm $fimo_filter`;

	%ori_mat = %mat;
	%max_mat = %mat;
	%remat = %mat;
	if( $num eq 0){
		${$remat{0}}[3] = $num;
		${$remat{0}}[4] = $num;
		${$remat{0}}[5] = $fi_num;
	}else{
		$k = 0;
		$w1 = 1 - $p;
		$w2 = $p/$num;
		$sta = &motif_inf(\%mat);
		
		$max = $sta;
		foreach $keys (sort {${$mot_seq{$b}}[0] <=> ${$mot_seq{$a}}[0]} keys %mot_seq){
			$k = $k + 1;
			@sen1 = split(//, uc(${$mot_seq{$keys}}[2]));
			for ($j=1; $j<=$len; $j++){
				$remat{$j} = [0, 0, 0, 0];
				if($sen1[$j-1] eq "A"){
					${$remat{$j}}[0] = ${$remat{$j}}[0] + 1;
				}elsif($sen1[$j-1] eq "C"){
					${$remat{$j}}[1] = ${$remat{$j}}[1] + 1;
				}elsif($sen1[$j-1] eq "G"){
					${$remat{$j}}[2] = ${$remat{$j}}[2] + 1;
				}else{
					${$remat{$j}}[3] = ${$remat{$j}}[3] + 1;
				}
			}
			for ($j=1; $j<=$len; $j++){
				for ($i=0; $i<=3; $i++){
					${$remat{$j}}[$i] = $w2/($w1 + $w2) * ${$remat{$j}}[$i] + $w1/($w1 + $w2) * ${$mat{$j}}[$i];
				}
			}		
			$rem = &motif_inf(\%remat);
			if (&motif_n50_inf(\%mat, \%remat) > 1){
				%max_mat = %remat;
				$max = $rem;
				%mat = %remat;
			}elsif(($rem >= $sta*(1-$limit)) && ($rem >= $max*(1-$limit)) && ($rem >= 0.5) && (&motif_n50_inf(\%ori_mat, \%remat) > (1-$limit2)) && (&motif_n50_inf(\%max_mat, \%remat) > (1-$limit2)) && (&motif_n50_per(\%ori_mat, \%remat) > (1-$limit2)) && (&motif_n50_per(\%max_mat, \%remat) > (1-$limit2))){
				%mat = %remat;
			}else{
				last;
			}
			$w1 = $w1 + $p/$num;
		}
		${$remat{0}}[3] = $k;
		${$remat{0}}[4] = $num;
		${$remat{0}}[5] = $fi_num;
		${$remat{0}}[6] = $rem;
	}
        `rm $fimo_out`;
	
	return (\%remat);
}

#  $rem = &motif_inf(\%remat);
sub motif_inf{
	my $mot_inf = shift;
	my %mot_inf = %{$mot_inf};
	my $i = 0; my $j; my $key; my $sum = 0.0; my $len = ${$mot_inf{0}}[2];
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		$sum = $sum + &INF(${$mot_inf{$key}}[0], ${$mot_inf{$key}}[1], ${$mot_inf{$key}}[2], ${$mot_inf{$key}}[3]);
	}
	$sum = $sum/$i;
	return $sum;
}

#  &motif_n50_inf(\%ori_mat, \%remat)
sub motif_n50_inf{
	my $inf1 = shift; my $inf2 = shift;
	my %mot_inf1 = %{$inf1}; my %mot_inf2 = %{$inf2}; my %pinf = ();
	my $i = 0; my $j; my $key; my $sum = 0.0; my $len = ${$mot_inf1{0}}[2]; my $s1 = 0.0; my $s2 = 0.0;
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		$pinf{$key}= [&INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]), &INF(${$mot_inf2{$key}}[0], ${$mot_inf2{$key}}[1], ${$mot_inf2{$key}}[2], ${$mot_inf2{$key}}[3])];
		$sum = $sum + &INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]);
	}
	$sum = $sum/2;
	foreach $key ( sort {${$pinf{$b}}[0] <=> ${$pinf{$a}}[0]} keys %pinf){
		$s1 = $s1 + ${$pinf{$key}}[0];
		$s2 = $s2 + ${$pinf{$key}}[1];
		if($s1 >= $sum){
			last;
		}
	}
	return $s2/$s1;
}

#  &motif_n50_per(\%ori_mat, \%remat) 
sub motif_n50_per{
	my $inf1 = shift; my $inf2 = shift;
	my %mot_inf1 = %{$inf1}; my %mot_inf2 = %{$inf2}; my %pinf = ();
	my $i = 0; my $j; my $key; my $sum1 = 0.0; my $sum2 = 0.0; my $len = ${$mot_inf1{0}}[2]; my $s1 = 0.0; my $s2 = 0.0;
	for ($key=1; $key<=$len; $key++){
		$i = $i + 1;
		$pinf{$key}= [&INF(${$mot_inf1{$key}}[0], ${$mot_inf1{$key}}[1], ${$mot_inf1{$key}}[2], ${$mot_inf1{$key}}[3]), &INF(${$mot_inf2{$key}}[0], ${$mot_inf2{$key}}[1], ${$mot_inf2{$key}}[2], ${$mot_inf2{$key}}[3])];
		$sum1 = $sum1 + ${$pinf{$key}}[0];
		$sum2 = $sum2 + ${$pinf{$key}}[1];
	}
	foreach $key ( sort {${$pinf{$b}}[0] <=> ${$pinf{$a}}[0]} keys %pinf){
		$s1 = $s1 + ${$pinf{$key}}[0];
		$s2 = $s2 + ${$pinf{$key}}[1];
		if($s1 >= $sum1/2){
			last;
		}
	}
	return ($s2/$sum2)/($s1/$sum1);
}


#DifferentialDistributionAnalysis(negative bonuly distribution). INPUT: two references of two array. OUT: an array 
sub DifferentialDistributionAnalysis{
    my $arr1=shift; my @ARR1=@{$arr1};
    my $arr2=shift; my @ARR2=@{$arr2};
    my @outARR=();
    my $sum1=0; my $sum2=0;
    for(@ARR1){ $sum1+=$_; }
    for(@ARR2){ $sum2+=$_; }
    my $p25=$sum2/$sum1;

    for(0..$#ARR1){
       my $i=$_;
       my $p=0;
       my $r=$ARR2[$i]; my $x=$ARR1[$i];
       if($r==0){
          $outARR[$i]=1;
       }elsif($r==1){
          for(1..$x){
              my $j=$_;
              $p+= $p25*(1-$p25)**($j-1);
          }
          $outARR[$i]=$p;
       }else{
          $p=$p+$p25**($r);      ### if $r==$x, the p just equals to the "$p" [ $p25**($r) ]
          if($r<$x){
            for(($r+1)..$x){
               my $k=$_;
               my $tmpsumA=0; my $tmpsumB=0;
               for(($k-$r+1)..($k-1)){ $tmpsumA+=log($_)/log(10); }
               for(1..($r-1)){ $tmpsumB+=log($_)/log(10); }
               #my $a = Math::BigFloat->new("$p25"); my $b = Math::BigFloat->new("$r"); my $ab=sprintf("%e", $a**$b);
               #my $ten=Math::BigFloat->new("10");
               #my $p75 = 1-$p25; my $kr= $k-$r;
               #my $x = Math::BigFloat->new("$p75"); my $y = Math::BigFloat->new("$kr"); my $xy=sprintf("%e", $x**$y);
               ###$a=sprintf("%e", $p25**($r)); $b=sprintf("%e", (1-$p25)**($k-$r)); print $p25, " <-> ", $r,"\n",$a," <=> ",$b,"\n";
               $p+= 10**( $tmpsumA - $tmpsumB + $r*log($p25)/log(10) + ($k-$r)*log(1-$p25)/log(10) );
               #$p+= 10**( $tmpsumA - $tmpsumB + log($a**$b)/log(10) + log($x**$y)/log(10) );
            }
          }
          $outARR[$i]=$p;
       }
    }

    return @outARR;
}

#uniq a array. INPUT: a reference of array. OUTPUT: an array.
sub uniq_array{
    my $arr=shift; my @color=@{$arr};
    my (%seen, @unique);
    foreach my $value (@color) {
        if (!$seen{$value}) {
           push @unique, $value;
           $seen{$value} = 1;
        }
    }
    return @unique;
}

#retreiving shared elements(maybe duplicated) from two redundant arrays, and calculating intersection number. INPUT: two reference of arrays to be compared.
sub intersect_two_redun_array{
    my $A=shift; my $B=shift; my $common=0;
    my (%counta, %countb, @common);
    for my $a_ele (@{$A}){$counta{$a_ele}++;}
    for my $b_ele (@{$B}){$countb{$b_ele}++;}
    foreach my $ele (keys %counta){
        if(defined($countb{$ele})){
          push @common, $ele;
          $common+=min($counta{$ele},$countb{$ele});
          }
    }
    return (\@common, $common);
}

#read standard fasta file to hash. INPUT: fasta file (absolute path). OUTPUT: a reference of hash.
sub readfa{
    my $fa=shift; my (%fa, $id);
    open(FA,$fa);
    while(<FA>){
        chomp;
        if(/^>(\S+)/){
             $id=$1;
        }else{
             $fa{$id}.=$_;
        }
    }
    close FA;
    return \%fa;
}

#calculating euc-distance from two reference arrays of the same length. INPUT: two reference of arrays.
sub eucdist {
    my $dist;
    my $l_dis=shift; my $r_dis=shift; 
    my @lef_dis=@{$l_dis}; my @rit_dis=@{$r_dis};
        for (my $index = 0; $index < scalar(@lef_dis); $index++)
            { $dist += ($lef_dis[$index] - $rit_dis[$index]) * ($lef_dis[$index] - $rit_dis[$index]); }
        $dist = sprintf "%.3f", sqrt($dist)/scalar(@lef_dis);
        return $dist;
}


#Dynamic Programming for pair of refernce arrays containing multiple fimo-like lines (belong to one "BLOCK" like FIR000), such as:
####################   5 FIR     41      47      +       13.1143 4.39e-05                GCTTGGC     ############### 
#INPUT: two reference of arrays to calculate their "MATCH" score
sub dynprog_matrix_pair_fimo{

    my $A=shift; my $B=shift;
    my ($match,$mismatch)=(1,-1);
    my (@tempA1,@tempA2,@singleMotifA1,@singleMotifA2);
    @tempA1=@{$A}; @tempA2=@{$B};

    my (%A1Motif,%A2Motif,%A1Prob,%A2Prob,%count,%A1count,%A2count,);
    for (0..$#tempA1){ @singleMotifA1=split("\t", $tempA1[$_]); $A1Motif{$_}=$singleMotifA1[0]."_".$singleMotifA1[4]; $A1Prob{$_}=log($singleMotifA1[6])/log(10); }
    for (0..$#tempA2){ @singleMotifA2=split("\t", $tempA2[$_]); $A2Motif{$_}=$singleMotifA2[0]."_".$singleMotifA2[4]; $A2Prob{$_}=log($singleMotifA2[6])/log(10); }

    my ($SizeScale,$GapScale)=(1,0.5);   ### how to combine the Prob size (e.g. e-5) and Prob gap (or distance of two Prob values)  ###

    foreach my $A1ranknum (sort keys %A1Motif){
     foreach my $A2ranknum (sort keys %A2Motif){
        if ($A1Motif{$A1ranknum} eq $A2Motif{$A2ranknum}){
            my $ProbGap=abs($A1Prob{$A1ranknum}-$A2Prob{$A2ranknum});
            my $ProbSize=($A1Prob{$A1ranknum}>$A2Prob{$A2ranknum})?abs($A2Prob{$A2ranknum}):abs($A1Prob{$A1ranknum});
            my $SizeMinusGapScore=$SizeScale*$ProbSize-$GapScale*$ProbGap;
            if (defined($count{$A1Motif{$A1ranknum}})){ if ($SizeMinusGapScore < $count{$A1Motif{$A1ranknum}}){
                                                           $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
            }else{ $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
        } else {
            my $A1ProbSize=abs($A1Prob{$A1ranknum});
            if (defined($A1count{"A1MaxMotif"})){ if ($A1ProbSize > $A1count{"A1MaxMotif"}){ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            }else{ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            my $A2ProbSize=abs($A2Prob{$A2ranknum});
            if (defined($A2count{"A2MaxMotif"})){ if ($A2ProbSize > $A2count{"A2MaxMotif"}){ $A2count{"A2MaxMotif"}=$A2ProbSize; }
            }else{ $A2count{"A2MaxMotif"}=$A2ProbSize; }
        }
     }
    }

    my $score;
    if (!%count){
       if ((%A1count) && (%A2count)){
         $score=$SizeScale*max($A1count{"A1MaxMotif"},$A2count{"A2MaxMotif"});
         $score=$mismatch*$score;
       }
    }else {
       foreach my $key (keys %count){
         $score+=$count{$key};
       }
       $score=$match*$score;
    }

    return sprintf "%.3f",$score;
    undef %count; undef %A1count; undef %A2count;

}

#Dynamic Programming for pair of refernce arrays, respectively containing (one) selected motif types [such as: FIR000_5_+_log(4.39e-05)...], and (another) multiple fimo-like lines (belong to one "BLOCK" like FIR000), such as:
#####################   5 FIR     41      47      +       13.1143 4.39e-05                GCTTGGC     ############### 
##INPUT: two reference of arrays to calculate their "MATCH" score
sub dynprog_matrix_selected_and_fimo{

    my $A=shift; my $B=shift;
    my ($match,$mismatch)=(1,-1);
    my (@tempA1,@tempA2,@singleMotifA1,@singleMotifA2);
    my (%A1Motif,%A2Motif,%A1Prob,%A2Prob,%count,%A1count,%A2count);
    @tempA1=@{$A}; @tempA2=@{$B};

    for (0..$#tempA1){ @singleMotifA1=split("_", $tempA1[$_]); $A1Motif{$_}=$singleMotifA1[1]."_".$singleMotifA1[2]; $A1Prob{$_}=$singleMotifA1[3]; }
    for (0..$#tempA2){ @singleMotifA2=split("\t", $tempA2[$_]); $A2Motif{$_}=$singleMotifA2[0]."_".$singleMotifA2[4]; $A2Prob{$_}=log($singleMotifA2[6])/log(10); }

    my ($SizeScale,$GapScale)=(1,0.5);   ### how to combine the Prob size (e.g. e-5) and Prob gap (or distance of two Prob values)  ###

    foreach my $A1ranknum (sort keys %A1Motif){
     foreach my $A2ranknum (sort keys %A2Motif){
        if ($A1Motif{$A1ranknum} eq $A2Motif{$A2ranknum}){
            my $ProbGap=abs($A1Prob{$A1ranknum}-$A2Prob{$A2ranknum});
            my $ProbSize=($A1Prob{$A1ranknum}>$A2Prob{$A2ranknum})?abs($A2Prob{$A2ranknum}):abs($A1Prob{$A1ranknum});
            my $SizeMinusGapScore=$SizeScale*$ProbSize-$GapScale*$ProbGap;
            if (defined($count{$A1Motif{$A1ranknum}})){ if ($SizeMinusGapScore < $count{$A1Motif{$A1ranknum}}){
                                                           $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
            }else{ $count{$A1Motif{$A1ranknum}}=$SizeMinusGapScore; }
        } else {
            my $A1ProbSize=abs($A1Prob{$A1ranknum});
            if (defined($A1count{"A1MaxMotif"})){ if ($A1ProbSize > $A1count{"A1MaxMotif"}){ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            }else{ $A1count{"A1MaxMotif"}=$A1ProbSize; }
            my $A2ProbSize=abs($A2Prob{$A2ranknum});
            if (defined($A2count{"A2MaxMotif"})){ if ($A2ProbSize > $A2count{"A2MaxMotif"}){ $A2count{"A2MaxMotif"}=$A2ProbSize; }
            }else{ $A2count{"A2MaxMotif"}=$A2ProbSize; }
        }
     }
    }

    my $score;
    if (!%count){
       if ((%A1count) && (%A2count)){
         $score=$SizeScale*max($A1count{"A1MaxMotif"},$A2count{"A2MaxMotif"});
         $score=$mismatch*$score;
       }
    }else {
       foreach my $key (keys %count){
         $score+=$count{$key};
       }
       $score=$match*$score;
    }

    return sprintf "%.3f",$score;

}

#Normalizing scoring matrix. INPUT: a reference of 2-D hash of matrix of "MATCH" score
sub matrix_norm{

    my $macro_matrix=shift; my %_macro_matrix=%{$macro_matrix};
    my ($_M, $_m)=(6, -6);
    my $slop=1; ######### bigger $slop, more step-like; when $slop=1, sigmoid min/max key -> value: -6 -> 0 and 6->1  ###########
    my (@t_array, %Macro_matrix);
    my $Macro_matrix=\%Macro_matrix;
    for (values %_macro_matrix){ push @t_array, (values %{$_});}
    my ($_max, $_min)=(max(@t_array), min(@t_array));
    ($_max, $_min)=(10, -10);    #renew to constant values
    foreach my $k1 (sort {$a cmp $b} keys %_macro_matrix){
      my %temphash=%{$_macro_matrix{$k1}};
      foreach my $k2 (sort {$a cmp $b} keys %temphash){
        #$Macro_matrix->{$k1}{$k2}= sprintf "%.3f", tanh(($temphash{$k2}-$_min)/($_max-$_min)*($_M-$_m)+$_m);        ##### tanh transformation
        #$Macro_matrix->{$k1}{$k2}= sprintf "%.3f", 1/(1+exp(-$slop*(($temphash{$k2}-$_min)/($_max-$_min)*($_M-$_m)+$_m)))+$_m/($_M-$_m);   ##### sigmoid transformation to range of max and min values
        $Macro_matrix->{$k1}{$k2}= sprintf "%.3f", 2*1/(1+exp(-$slop*((($temphash{$k2}-4)-$_min)/($_max-$_min)*($_M-$_m)+$_m)))-1;       ##### sigmoid transformation to range of constant (-10,10)
      }
    }
    
    return \%Macro_matrix;

}

#Generating syntenic motifs with high scores. 
#INPUT: two variables of string like "FIR000FIR001...", and 
#       a reference hash of scoring matrix (suggesting a normalization before use) 
#       a reference hash of motif groups (merged "blocks") [ such as: FIR000 => qw(1_+ 2_- 4+) ]
sub macro_dynprog{

   my $seq1=shift; my $seq2=shift;
   my $macro_matrix=shift; my $Macro_matrix={%{$macro_matrix}};
   my $hash_motif_blocks=shift; my %_hash=%{$hash_motif_blocks};
   my ($GAP, $MATCH, @matrix);

   $matrix[0][0]{score}   = 0;
   $matrix[0][0]{pointer} = "none";
   for(my $j = 1; $j <= length($seq1)/6; $j++) {
     $matrix[0][$j]{score}   = 0;
     $matrix[0][$j]{pointer} = "none";
   }
   for(my $i = 1; $i <= length($seq2)/6; $i++) {
     $matrix[$i][0]{score}   = 0;
     $matrix[$i][0]{pointer} = "none";
   }

   my %cmp_hash;
   for(my $i = 1; $i <= length($seq2)/6; $i++) {
     for(my $j = 1; $j <= length($seq1)/6; $j++) {
        my ($diagonal_score, $left_score, $up_score);

        my $letter1 = substr($seq1, ($j-1)*6, 6);
        my $letter2 = substr($seq2, ($i-1)*6, 6);
           $MATCH=$Macro_matrix->{$letter1}{$letter2};
           $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;

        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;

        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next;
        }

        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                if ($diagonal_score >0){
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        } else {
            if ($up_score >= $left_score) {
                if ($up_score >0){
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
                }
            }
            else {
                if ($left_score >0){
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
                }
            }
        }

        $cmp_hash{$i."_".$j} = $matrix[$i][$j]{score};
     }
   }

   # trace-back
   my $max_i     = 0;
   my $max_j     = 0;

   my (@pool, @align1, @align2, $result);

   foreach my $_key_ (sort { $cmp_hash{$b} <=> $cmp_hash{$a} } keys %cmp_hash){

   unless ($_key_~~ @pool){
     ($max_i, $max_j)=split("_", $_key_);
     if($matrix[$max_i][$max_j]{pointer} eq "diagonal"){
        my $j = $max_j;
        my $i = $max_i;

        while (1) {
          last if $matrix[$i][$j]{pointer} eq "none";
          push @pool, $i."_".$j;
          if ($matrix[$i][$j]{pointer} eq "diagonal") {
             push @align1, substr($seq1, ($j-1)*6, 6);
             push @align2, substr($seq2, ($i-1)*6, 6);
             $i--; $j--;
          }elsif ($matrix[$i][$j]{pointer} eq "left") {
             push @align1, substr($seq1, ($j-1)*6, 6);
             push @align2, "------";
             $j--;
          }elsif ($matrix[$i][$j]{pointer} eq "up") {
             push @align1, "------";
             push @align2, substr($seq2, ($i-1)*6, 6);
             $i--;
          }
        }

        @align1 = reverse @align1;
        @align2 = reverse @align2;

        $result .= ">>>Score".$cmp_hash{$_key_}."\n".join("",@align1)."\n".join("",@align2)."\n\n";

        my (@motifID, $sys);

        foreach (0..$#align1){

          $sys=$_;
          $result.=$align1[$sys].":";

          if (defined ($_hash{$align1[$sys]})){
             foreach (@{$_hash{$align1[$sys]}}){
               @motifID=split("\t",$_);
               $result.=$motifID[0]."_".$motifID[4].",";
             }
          }else{
             $result.="------";
          }

          $result.="\t\t\t";
          $result.=$align2[$sys].":";

          if (defined ($_hash{$align2[$sys]})){
             foreach (@{$_hash{$align2[$sys]}}){
               @motifID=split("\t",$_);
               $result.=$motifID[0]."_".$motifID[4].",";
             }
          }else{
             $result.="------";
          }

          $result.="\n";

        }

        $result.="\n";

        undef @align1;
        undef @align2;
     }
   }

   }

   return ($result, @matrix);

}

#input pairwise filtered fimolines (in hash), and a threshold to cluster neighbour fimolines (neighboring overlap ratio) into blocks, and a value 0 or 1 as a flag to change run mode
#INPUT: a reference of hash, and a value (0,1,2), 0.5 or 0.8 recommended, as well as a vaule 0 (output two ref hash) or 1 (output three ref hash) or 2 (output four ref hash).
sub fimolines2blocks{
#open OUT, ">>/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/MicroH/CaseG/gDIR001/NONHSAT142366/DIR_NONHSAT142366_ENSMUST00000182174/DIR_NONHSAT142366-linker-ENSMUST00000182174/cluster/ttt";
  my (%t_hash,%t__hash); 
  my $inhash=shift; my $signal=shift;
  #my $thresholdratio=shift; my $pRate=shift;
  my $thresholdratio=0.8; my $pRate=1;
  my %inhash=%{$inhash};
  my %motif_hash if ($signal>=1);
  my %motif_hash_prob if ($signal>=2);

  foreach my $key (sort keys %inhash){
     my @ntmp=@{$inhash{$key}};
#print OUT "PREmotifpos[-->start\n", join("\n", @ntmp), "\n-->end]\n";
     my $temp=motifpos(\@ntmp, 1, 1); 
     my @temp=@{$temp}; 
     #my @temp=@ntmp;
#print OUT "\nmotifpos[-->start\n", join("\n", @temp), "\n-->end]\n";
     my @first=split(/\t/,$temp[0]);
     my ($sp,$st,$ed)=($first[1],$first[2],$first[3]);
     my $i="000";
      push @{$t_hash{$key.$i}}, $temp[0];                                  
      push @{$motif_hash{$key.$i}}, $first[0]."_".$first[4] if ($signal);
      push @{$motif_hash_prob{$key.$i."_".$first[0]."_".$first[4]}}, sprintf "%.3f", log($first[6])/log(10) if ($signal==2); 
      $t__hash{$key.$i}=$sp."_".$st."_".$ed;
     foreach (1..$#temp){
     my @tempArray=split(/\t/,$temp[$_]);
       if (($tempArray[1] eq $sp) && ($tempArray[2] <= $ed)){
          if (motif_include($thresholdratio, $pRate, \@tempArray, \@{$t_hash{$tempArray[1].$i}})){
              push @{$t_hash{$key.$i}}, $temp[$_];                
              push @{$motif_hash{$key.$i}}, $tempArray[0]."_".$tempArray[4] if ($signal);
              push @{$motif_hash_prob{$key.$i."_".$tempArray[0]."_".$tempArray[4]}}, sprintf "%.3f", log($tempArray[6])/log(10) if ($signal==2);
             ($sp,$st,$ed)=($tempArray[1],min($st,$tempArray[2]),max($ed,$tempArray[3]));
              $t__hash{$key.$i}=$sp."_".$st."_".$ed;
          }else{
             $i++;
              push @{$t_hash{$key.$i}}, $temp[$_];             
              push @{$motif_hash{$key.$i}}, $tempArray[0]."_".$tempArray[4] if ($signal);
              push @{$motif_hash_prob{$key.$i."_".$tempArray[0]."_".$tempArray[4]}}, sprintf "%.3f", log($tempArray[6])/log(10) if ($signal==2);
              $t__hash{$key.$i}=$tempArray[1]."_".$tempArray[2]."_".$tempArray[3];
             ($sp,$st,$ed)=($tempArray[1],$tempArray[2],$tempArray[3]);
          }
       }else{
          $i++;
           push @{$t_hash{$key.$i}}, $temp[$_];            
           push @{$motif_hash{$key.$i}}, $tempArray[0]."_".$tempArray[4] if ($signal);
           push @{$motif_hash_prob{$key.$i."_".$tempArray[0]."_".$tempArray[4]}}, sprintf "%.3f", log($tempArray[6])/log(10) if ($signal==2);
           $t__hash{$key.$i}=$tempArray[1]."_".$tempArray[2]."_".$tempArray[3];
          ($sp,$st,$ed)=($tempArray[1],$tempArray[2],$tempArray[3]);
       }
     }
   }

  if ($signal==1){
    return (\%t_hash, \%t__hash, \%motif_hash);
  }elsif ($signal==2){
    return (\%t_hash, \%t__hash, \%motif_hash, \%motif_hash_prob);
  }else {
    return (\%t_hash, \%t__hash, "OK");
  }

}

#assess if a fimo line (in array) can be included into a "already defined cluster/block"(in hash, value of which is an array of fimolines)
#INPUT: threshold, candidate fimoline(split into an array), and a refernce of hash
sub motif_include{
  my $_ratio=shift; my $pRate=shift;
  my $cand_arr=shift; my @cand_arr=@{$cand_arr};
  my ($cand_s, $cand_e)=($cand_arr[2], $cand_arr[3]);
  my $cand_len=$cand_e-$cand_s+1;
  my $target_arrlib=shift; my @target_arrlib=@{$target_arrlib};
  my ($overlaplen, $resultSignal1, $resultSignal0)=(0, 0, 0);

  for (@target_arrlib){
    my @temp=split(/\t/, $_);
    my ($target_s, $target_e)=($temp[2], $temp[3]);
    my $target_len=$target_e-$target_s+1;

    if ($target_s<=$cand_s){
        $overlaplen=($target_e<=$cand_e)?$target_e-$cand_s+1:$cand_e-$cand_s+1;
    }else{
        $overlaplen=($target_e<=$cand_e)?$target_e-$target_s+1:$cand_e-$target_s+1;
    }

    if ( ($overlaplen/$cand_len>=$_ratio)||($overlaplen/$target_len>=$_ratio) ){
       $resultSignal1++;
    }else{
       $resultSignal0++;
    }
  }

  return ( $resultSignal1/($resultSignal1+$resultSignal0)>=$pRate )?1:0;
}

#
######################   5 FIR     41      47      +       13.1143 4.39e-05                GCTTGGC     ###############
sub motifpos{

  my $fimoArr=shift; my @fimoArr=@{$fimoArr}; my $t_ratio=shift; my $pRate=shift;
  my (%hash, %sorthash, %keepfimo, @outArr);
  for my $ln (@fimoArr){
      my @ln=split(/\t/,$ln);
      $hash{$ln[1]}{$ln[0]}{$ln}{"score"}=$ln[5];
      $hash{$ln[1]}{$ln[0]}{$ln}{"pos"}=$ln[2];
      $hash{$ln[1]}{$ln[0]}{$ln}{"strand"}=$ln[4];
      push @{$sorthash{$ln[1]}{$ln[0]}}, $ln;
  }
#open OUT, ">>/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/tttmp";
  foreach my $sp (keys %hash){
  foreach my $motifid (keys %{$hash{$sp}}){
      my %tmphash = %{$hash{$sp}{$motifid}};
	  my @tmparr  = @{$sorthash{$sp}{$motifid}};
	  #foreach my $fimoline (sort{$tmphash{$b}<=>$tmphash{$a}} keys %tmphash){
	  foreach my $fimoline (sort { $tmphash{$b}{"score"} <=> $tmphash{$a}{"score"} or
                                       $tmphash{$a}{"pos"} <=> $tmphash{$b}{"pos"} or
                                       $tmphash{$a}{"strand"} cmp $tmphash{$b}{"strand"} or
                                       $a cmp $b
                                     } keys %tmphash){
	  if ( $fimoline ~~ @tmparr ){
		   my @arr_range; push @arr_range, $fimoline;
		   my @arr_range_idx;
		   my @fimoline=split("\t", $fimoline);
               $keepfimo{$sp}{$fimoline}=1;
	       #$keepfimo{$fimoline}{"k1"}=$fimoline[2];
               #$keepfimo{$fimoline}{"k2"}=$fimoline[3];
           my $idx = getArrIdx($fimoline, \@tmparr); push @arr_range_idx, $idx;
           my $dist=max($idx-0, $#tmparr-$idx);
#print OUT "0=>", $idx,"=>",$fimoline,"\n";
           for(my $i=1; $i<=$dist; $i++){
             if($idx-$i>=0){
	         my $lef_idx=$idx-$i;
                 if (motif_include($t_ratio, $pRate, \@{[split("\t", $tmparr[$lef_idx])]}, \@arr_range)){
		      push @arr_range, $tmparr[$lef_idx];
		      push @arr_range_idx, $lef_idx;
#print OUT "1=>", $lef_idx,"\n";
		 }
	     }
	     if($idx+$i<=$#tmparr){
	         my $rit_idx=$idx+$i;
                 if (motif_include($t_ratio, $pRate, \@{[split("\t", $tmparr[$rit_idx])]}, \@arr_range)){
		      push @arr_range, $tmparr[$rit_idx];
		      push @arr_range_idx, $rit_idx;
#print OUT "2=>", $rit_idx,"\n";
		 }
	     }
           }
#print OUT "3=>", scalar(@tmparr),"=>",join(" ", @arr_range_idx),"\n";
	   for(sort { $b <=> $a } @arr_range_idx){
	      splice(@tmparr, $_, 1);
           }
	  }
	  }
  }
  }
  #@outArr = sort {
  #      $keepfimo{$a}{"k1"} <=> $keepfimo{$b}{"k1"} or
  #      $keepfimo{$a}{"k2"} <=> $keepfimo{$b}{"k2"} or
  #      $a cmp $b
  #} keys %keepfimo;

  foreach my $mykey (@fimoArr){
     my @F=split(/\t/,$mykey);
     if (defined $keepfimo{$F[1]}{$mykey}){
        push @outArr, $mykey;
     }
  }

  return \@outArr;

}


###
sub getArrIdx{
  my $element=shift; my $arr=shift; my @arr=@{$arr};
  my $idx;
  for(0..$#arr){
     if ("$arr[$_]" eq "$element"){
	    $idx=$_;
	 }
  }
  return $idx;
}


#
#sub {
#
#}



1;
