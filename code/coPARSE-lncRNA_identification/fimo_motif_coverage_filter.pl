#!/usr/bin/perl -w 

#BEGIN {
#  unshift @INC, "/Share/home/zhangqf2/xiongtl/xiongtl/data/annotation/H_M_lncFIMO/Gencode_MicroH/FinalRun/ScrDB/scripts";
#}
use List::Util qw(min max);
use Myownperlmodule qw(Marray2fimo fimo2Marray motifpos);

$FimoFilterSort = shift;

$readfimo=fimo2Marray($FimoFilterSort, "F");
$tmpArr=motifpos(\@{$readfimo}, 0.01, 1);
$FimoFilterSort=Marray2fimo(\@{$tmpArr}, $FimoFilterSort);

$motifCovRatio = 10; $motifNum = 10;

open INF, "<$FimoFilterSort";
while (<INF>){
    chomp;
    if (!/^#/){
    @F=split("\t", $_);
    push @{$hashsortLn{$F[1]}}, $F[3]-$F[2]+1;
    push @{$hashsortSt{$F[1]}}, $F[2];
    push @{$hashsortEd{$F[1]}}, $F[3];
    }
}

@mykeys = keys %hashsortLn;
if ( scalar(@mykeys) >1 ){
   $mykey = ( scalar(@{$hashsortLn{$mykeys[0]}}) < scalar(@{$hashsortLn{$mykeys[1]}}) ) ? $mykeys[0] : $mykeys[1];
}else{
   $mykey = $mykeys[0];
}

@Ln = @{$hashsortLn{$mykey}};
@St = @{$hashsortSt{$mykey}};
@Ed = @{$hashsortEd{$mykey}};

#$CovRange = max(@Ed)-min(@St)+1;

$Sum=0; $Num=0;
foreach ( @Ln ){ $Sum+=$_; $Num+=1; }
$MotifAveLen=$Sum/$Num;

=head
$s=$St[0]; $e=$Ed[0];
$UniqCov=$e-$s+1;
foreach ( 1..$#Ln ){ 
    $mys=$St[$_]; 
    $mye=$Ed[$_];
    if ($mys <= $e){
        $UniqCov+=($mye<=$e)?0:($mye-$e+1);
        $e=max($e, $mye);
    }else {
        $UniqCov+=$mye-$mys+1;
        $s=$mys; $e=$mye;
    }
}
=cut

%covArr=();
foreach ( 0..$#Ln ){ 
    for ($i=$St[$_]; $i<=$Ed[$_]; $i++){
         $covArr{$i}="YES";
    }
}
for ($j=min(@St); $j<=max(@Ed); $j++){
    if ((defined $covArr{$j}) && ($covArr{$j} eq "YES")){
        $UniqCov+=1;
    }
} 

if ( ($UniqCov >= $motifCovRatio * $MotifAveLen) && ($Num >= $motifNum) ){              #CovRange > $UniqCov
   #print "$CovRange, $Sum, $Num, $MotifAveLen, $UniqCov\n";
   #print "Yes\n";
}else{
   #print "$CovRange, $Sum, $Num, $MotifAveLen, $UniqCov\n";
   #print "No\n";
   #`rename ".sort" ".SORT.DEL" $FimoFilterSort`;
   `rm -rf $FimoFilterSort`;
}
