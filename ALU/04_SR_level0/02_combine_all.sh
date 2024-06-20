outpath=$1
ver=$2
masterpath=$3
sub=$4
while read line; do
    $masterpath/ALU/04_SR_level0/03_combine_SR.pl $line $ver $outpath
done < $outpath/$sub/list 

