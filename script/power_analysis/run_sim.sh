########## use ABC rejection to do power analysis on simulated data #########
# generate simulations
pdir=./
only_mut=0
for ncell in {20..60..20}
do
    echo $ncell
    julia $pdir/generate_sim_sstat.jl $ncell 1 100000 $only_mut
    julia $pdir/generate_sim_sstat.jl $ncell 0 100000 $only_mut
done


# model selection with DIC
ntest=100
ndic=1
num_sample=100
s3_weight=2.0
only_mut=0
for ncell in {20..60..20}
do
    echo $ncell
    # neutral model
    julia compute_DIC_sim_abcrej_logdist_vscale_withbratio.jl $ntest $ndic $num_sample $ncell $s3_weight $only_mut

    # selection model
    for f in {2..8..2}
    do
        # fitness=`echo "scale=1; -$f/10" | bc`
        fitness=`awk -v var=$f 'BEGIN {printf "%.1f", -var/10}'`
        echo $fitness
        julia compute_DIC_sim_abcrej_logdist_vscale_withbratio.jl $ntest $ndic $num_sample $ncell $s3_weight $only_mut $fitness
    done
done
