pdir=./
odir=../../data/dic
msample=100
nsample=100


dataset=$1
ncell=$2
epsilon=$3
type=$4
s3_weight=$5

suffix="logs3_weight${s3_weight}"

fout=$odir/dic_${type}_1-100_${suffix}_$dataset
cat /dev/null > $fout

for i in {1..100..1}
do
    echo $i
    dic=`julia $pdir/compute_DIC_${type}_logdist.jl $dataset $ncell $epsilon $msample $nsample $s3_weight`
    echo "$i $dataset $ncell $epsilon $dic" >> $fout
done
