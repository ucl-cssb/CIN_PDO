using Distributions
using DelimitedFiles
using Distances
using DataFrames
using Statistics
# using StatsBase


# Compute DIC for ABC rejection on simulated data directly based on simulated summary statistics under different models


# For a given set of parameters, simulate 10 target datasets
ntest = parse(UInt32, ARGS[1])
# For a simulated target data, compute DIC for 100 times
ndic = parse(UInt32, ARGS[2])
num_sample = parse(UInt32, ARGS[3])
ncell = parse(UInt32, ARGS[4])
s3_weight = parse(Float64, ARGS[5])
target_fitness = 0
if(size(ARGS,1) > 5)
    target_fitness = parse(Float64, ARGS[6])
end


# Simulate target data
# ntest = 1
# ndic = 1
# num_sample = 10
# ncell = 20
# target_fitness = -0.2

println(ncell)
println(target_fitness)
println(s3_weight)

# agg_func = Statistics.median
# agg_func = median
agg_func = mean


# true values of the parameters (varing fitness values)
mu_sim1 = 0.1
mu_sim2 = 0.1
brate_sim = 0.4

# 200 used in the DIC paper
nsample = num_sample
msample = num_sample

nsims=100000
acc=0.001


max_mu = 0.5

min_brate = 0.3
max_brate = 0.5

# very unlikely values
LNL_MIN = -1e20

SMALL_VAL = 1e-9


# change directory accordingly
dir = "../../"
progSim = dir * "/code/simpdo"
simdir = dir * "/data/power_analysis/"


rsuffix = "_logs3_weight" * string(s3_weight) * ".tsv"

skip = 3


############################## Simulate target data ##############################
function getTargetdata_neutral(seed, mu_sim1, mu_sim2, brate_sim)
    mu1 = mu_sim1
    mu2 = mu_sim2
    birth_rate = brate_sim
    death_rate = 0

    suffix=string(seed)

    cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --suffix $suffix`

    res = chomp(read(cmdSim, String))
    arr = split(res,"\t")
    targetdata = map(x->parse(Float64,x), arr)

    targetdata = targetdata[1:3]

    return targetdata
end


function getTargetdata_withsel(seed, mu_sim1, mu_sim2, brate_sim, fitness)
    mu1 = mu_sim1
    mu2 = mu_sim2
    birth_rate = brate_sim
    death_rate = 0
    suffix=string(seed)

    fitval = "1 " * string(fitness)
    model = 1

    cmdSim=`$progSim -o $dir -e $ncell --model $model --fitness "$fitval" --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --suffix $suffix`

    res = chomp(read(cmdSim, String))
    arr = split(res,"\t")
    targetdata = map(x->parse(Float64,x), arr)

    targetdata = targetdata[1:3]

    return targetdata
end


############################## standard routine to compute DIC ##############################
# Compute surrogate likelihood
function compute_surrogate_likelihood(simdata, targetdata, scales)
    # likelihood = exp.(-0.5 * sqeuclidean(log.(targetdata), log.(simdata))) / sqrt(2 * pi)
    likelihood = exp.(-0.5 * sqeuclidean((targetdata), (simdata))) / sqrt(2 * pi)
    return likelihood
end


function compute_lnl_nsample(params, nsample, targetdata, scales)
  #println(params)
  mu1 = params[1]
  mu2 = params[2]
  birth_rate = params[3]
  death_rate = 0

  if size(params, 1) > 3
      fitval = "1 " * string(params[4])
      model = 1
  else
      fitval = "1 0"
      model = 0
  end
  # println(model)
  # println(fitval)

  sum_likelihood = 0.0
  for i in 1:nsample
      seed=rand(UInt64, 1)[1]
      suffix=string(seed)

      # call simulation program to generate data
      cmdSim=`$progSim -o $dir -e $ncell --model $model --fitness "$fitval" --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --suffix $suffix`

      res = chomp(read(cmdSim, String))
      arr = split(res,"\t")
      simdata = map(x->parse(Float64,x), arr)
      simdata = simdata[1:3]
      # println(simdata)
      # simdata = vcat(simdata[1:2], log(simdata[3]/ncell))

      # println("original simdata ", simdata)
      # simdata[simdata.==0] .= SMALL_VAL
      # simdata = vcat(log(simdata[1] + 1)/scales[1], log(simdata[2] + 1)/scales[2], log(simdata[3])/scales[3])
      # simdata = vcat(log(simdata[1])/scales[1], log(simdata[2])/scales[2], log(simdata[3])/scales[3])
      simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]))
      # println("log simdata ", simdata)

      ll = compute_surrogate_likelihood(simdata, targetdata, scales)
      #println(string(simdata) * " with likelihood " * string(ll))

      sum_likelihood += ll
  end

  #println("sum_likelihood of nsample: " * string(sum_likelihood))
  if sum_likelihood > 0
      lnl = log(sum_likelihood / nsample)
  else
      lnl = LNL_MIN
  end

  return lnl
end


function compute_dev_param(params, nsample, targetdata, scales)
    dev = -2 * compute_lnl_nsample(params, nsample, targetdata, scales)
    return dev
end


function compute_devbar(abcres, type, msample, nsample, targetdata, scales)
    sum_likelihood = 0.0
    resampled = sample(1:size(abcres,1), msample)

    for i in 1:msample
        sel = abcres[resampled[i],:]

        if type == "neutral"
            params = [sel[1], sel[2], sel[3]]
        else
            params = [sel[1], sel[2], sel[3], sel[4]]
        end
        lnl = compute_lnl_nsample(params, nsample, targetdata, scales)
        sum_likelihood += lnl
        # println("sum_likelihood for msample " * string(i) * " is " * string(sum_likelihood))
    end

    devbar = -2 * sum_likelihood / msample

    return devbar
end



function compute_dic(devbar, dev_param)
    dic = 2 * devbar - dev_param
    return dic
end



# type: the type under which ABC is run
# agg_func: can be mean or meadian
function get_dic_by_type(type, abcres, msample, nsample, targetdata, agg_func, scales)
    # Compute posterior mean
    pst_mu1 = agg_func(abcres[:,1])
    pst_mu2 = agg_func(abcres[:,2])
    pst_brate = agg_func(abcres[:,3])

    if type == "neutral"
        pparams = [pst_mu1, pst_mu2, pst_brate]
    else
        pst_fitness = agg_func(abcres[:,4])
        pparams = [pst_mu1, pst_mu2, pst_brate, pst_fitness]
    end
    # println(pparams)
    # println("original targetdata ", targetdata)
    # targetdata[targetdata.==0] .= SMALL_VAL
    # targetdata = vcat(log(targetdata[1] + 1)/scales[1], log(targetdata[2] + 1)/scales[2], log(targetdata[3])/scales[3])
    # targetdata = vcat(log(targetdata[1])/scales[1], log(targetdata[2])/scales[2], log(targetdata[3])/scales[3])
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]))
    # println("log targetdata ", targetdata)

    dev_param = compute_dev_param(pparams, nsample, targetdata, scales)
    devbar = compute_devbar(abcres, type, msample, nsample, targetdata, scales)
    dic = compute_dic(devbar, dev_param)
    # println(devbar)
    # println(dev_param)
    # println(dic)

    return dic
end

############################## Get DIC with ABC REJ ##############################

# Get variance of summary statistics from prior distributions (simulations)
# use log scale due to skewness and huge range of s3
function get_sstat_var(simdata_sstat,)
    simdata_all[simdata_sstat.==0] .= SMALL_VAL
    # sstart = 6
    # var1 = std(log.(simdata_all[:, sstart] .+ 1))
    # var2 = std(log.(simdata_all[:, sstart+1] .+ 1))
    var1 = std(log.(simdata_sstat[:, 1]))
    var2 = std(log.(simdata_sstat[:, 2]))
    # var1 = 1
    # var2 = 1
    var3 = std(log.(simdata_sstat[:, 3]))

    return  [var1, var2, var3]
end


# acceptance rate: acc
# should generate simulations for different number of cells
# separate abc from dic computation to check posterior distribution
function get_abc_res(type, simdata_all, sstart, scales, ncell, nsims, msample, nsample, targetdata, agg_func, acc)
    # println("running abc rej")
    # targetdata = vcat(log(targetdata[1] + 1)/scales[1], log(targetdata[2] + 1)/scales[2], log(targetdata[3])/scales[3])
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]))

    dim_stat = size(simdata_all, 2)  # dimension of simulation results
    # Compute the distance between simulated data and targetdata
    dists = zeros(Float64, nsims)
    for i in 1:nsims
        simdata1 = simdata_all[i, :]

        simdata = simdata1[sstart:sstart+2]
        # println("original simdata ", simdata)
        # simdata = vcat(log(simdata[1] + 1)/scales[1], log(simdata[2] + 1)/scales[2], log(simdata[3])/scales[3])
        simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]))

        # just used for sorting, no need to use log distance. Use original distance to be consistent with previous results. Biased with longer branches
        # println("log targetdata ", targetdata)
        # println("log simdata ", simdata)
        dist = evaluate(Euclidean(1e-16), (targetdata), (simdata))
        # dist = evaluate(Euclidean(1e-16), log.(targetdata), log.(simdata))
        dists[i] = dist
    end

    simdata_all = hcat(simdata_all, dists)

    # Sort data by distance and take top acceptance rate as rejection results
    ndim = dim_stat + 1
    naccept = Int(size(simdata_all, 1) * acc)
    abcres = sortslices(simdata_all, dims = 1, by = x -> x[ndim])[1:naccept,:]

    return abcres
end



# write results line by line to reduce memory usage and get results in time
fout = simdir * "/dic_ncell" * string(ncell) * "_fitness" * string(target_fitness) * "_ntest" * string(ntest) * "_ndic" * string(ndic) * "_nsample" * string(num_sample) * "_f" * string(agg_func) * "_acc" * string(acc) *  rsuffix

fsim_suffix = "_ncell" * string(ncell)  * "_n" * string(nsims) * "_maxmu" * string(max_mu) * "_minb" * string(min_brate) * "_maxb" * string(max_brate) * ".tsv"

type = "neutral"
fsim = simdir * "/sim_data_" * type * fsim_suffix
simdata_all_neutral = readdlm(fsim)
sstart_neutral = 5  # start index of summary statistics, from 5 when estimating birth rate

type = "selection"
fsim = simdir * "/sim_data_" * type * fsim_suffix
simdata_all_selection = readdlm(fsim)
sstart_selection = 6

# no need to do scaling when assuming neutral model
scales_neutral = [1, 1, 1]
scales_selection = [1, 1, 1]


if target_fitness == 0
    println("target data under neutral evolution")
    # res_all = zeros(Float64, ntest * ndic, 10)
    open(fout, "a") do f
        for i in 1:ntest
            seed = rand(UInt64, 1)[1]
            targetdata = getTargetdata_neutral(seed, mu_sim1, mu_sim2, brate_sim)
            # println(targetdata)

            for j in 1:ndic
                type = "neutral"
                abcres_neutral = get_abc_res(type, simdata_all_neutral, sstart_neutral, scales_neutral, ncell, nsims, msample, nsample, targetdata, agg_func, acc)
                dic_neutral = get_dic_by_type(type, abcres_neutral, msample, nsample, targetdata, agg_func, scales_neutral)

                type = "selection"
                abcres_withsel = get_abc_res(type, simdata_all_selection, sstart_selection, scales_selection, ncell, nsims, msample, nsample, targetdata, agg_func, acc)
                dic_withsel = get_dic_by_type(type, abcres_withsel, msample, nsample, targetdata, agg_func, scales_selection)

                # println(string(dic_neutral) * " " * string(dic_withsel))
                res = vcat(ncell, mu_sim1, mu_sim2, brate_sim, target_fitness, seed, targetdata, j, dic_neutral, dic_withsel)
                # println(res)

                k = (i - 1) * ndic + j
                println(k)
                # res_all[k,:] = res
                writedlm(f, transpose(res))
            end
        end
    end
else
    println("target data with selection")
    # res_all = zeros(Float64, ntest * ndic, 11)
    open(fout, "a") do f
        for i in 1:ntest
            seed = rand(UInt64, 1)[1]
            targetdata = getTargetdata_withsel(seed, mu_sim1, mu_sim2, brate_sim, target_fitness)
            # println(targetdata)
            for j in 1:ndic
                type = "neutral"
                abcres_neutral = get_abc_res(type, simdata_all_neutral, sstart_neutral, scales_neutral, ncell, nsims, msample, nsample, targetdata, agg_func, acc)
                dic_neutral = get_dic_by_type(type, abcres_neutral, msample, nsample, targetdata, agg_func, scales_neutral)

                type = "selection"
                abcres_withsel = get_abc_res(type, simdata_all_selection, sstart_selection, scales_selection, ncell, nsims, msample, nsample, targetdata, agg_func, acc)
                dic_withsel = get_dic_by_type(type, abcres_withsel, msample, nsample, targetdata, agg_func, scales_selection)

                # println(string(dic_neutral) * " " * string(dic_withsel))
                res = vcat(ncell, mu_sim1, mu_sim2, brate_sim, target_fitness, seed, targetdata, j, dic_neutral, dic_withsel)
                # println(res)

                k = (i - 1) * ndic + j
                println(k)
                # res_all[k,:] = res
                writedlm(f, transpose(res))
            end
        end
    end
end
