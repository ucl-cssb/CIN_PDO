using Distributions
using DelimitedFiles
using Distances
using CSV
using DataFrames
using Statistics
using StatsBase

# Compute DIC for SMC ABC inference under neutral model based on euclidean distance (four summary statistics)


# read from file
dataset=ARGS[1]
ncell=ARGS[2]
epsilon=ARGS[3]
msample = parse(UInt32, ARGS[4])
nsample = parse(UInt32, ARGS[5])
s3_weight = parse(Float64, ARGS[6])
only_mut = parse(UInt32, ARGS[7])


# dataset="dataset_20190409_2"
# ncell=13
# epsilon=0.1
# msample=2
# nsample=2


type="neutral"

# println(dataset)
# println(ncell)
# println(epsilon)


# change directory accordingly
dir = "../../"
progSim=dir * "/code/simpdo"
frealCnvec = dir * "/data/4trees/" * dataset * "_rvec.txt"
frealTree = dir * "/data/4trees/" * dataset * "_divtree.txt"
frealTree_ratio = dir * "/data/4trees/" * dataset * "_divtree_blenratio.txt"

dir_smc = dir * "/data/abcsmc" * type


suffix = dataset * "_n" * string(ncell) * "_epsilon" * string(epsilon) * "_logs3_weight" * string(s3_weight) * "_logs4_mut" * string(only_mut)
if only_mut == 1
    suffix = suffix * "_s4plus1"
end
fres =  dir_smc *  "/smc-joint-" * suffix


targetdataCnvec = readdlm(frealCnvec)[:, 1]


# skip header and the first three cell divisions
targetTree = readdlm(frealTree, skipstart=7)[:, 3]
tblen = sum(targetTree)*3/1440/2
push!(targetdataCnvec, tblen)


# skip header: "type"     "time"     "after"     "before"     "after21"     "after22"
targetTree_ratio = readdlm(frealTree_ratio)[2:end, :]
# remove rows with before == 0
ratios_all = targetTree_ratio[targetTree_ratio[:,4] .> 0,:]
if only_mut == 1
    ratios = ratios_all[ratios_all[:,1] .!= "N",:]
else
    ratios = ratios_all
end
if size(ratios, 1) != 0
    aratio = mean(ratios[:,4] ./ ratios[:,3])
    push!(targetdataCnvec, aratio)

    targetdata = targetdataCnvec
else
    exit()
end
# end

if only_mut == 1
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]), log(1+targetdata[4]))
else
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]), log(targetdata[4]))
end

# println(targetdata)

# ftargetdata = dir * "/realdata/" * dataset * "_stat.txt"
# writedlm(ftargetdata, transpose(targetdata))


smcres = DataFrame(CSV.File(fres, comment="#"));


# Compute posterior mean
pst_mu1 = mean(smcres.parameter1, weights(smcres.weights))
pst_mu2 = mean(smcres.parameter2, weights(smcres.weights))
pst_brate = mean(smcres.parameter3, weights(smcres.weights))

params = [pst_mu1, pst_mu2, pst_brate]


# Compute surrogate likelihood
function compute_surrogate_likelihood(simdata, targetdata)
    likelihood = exp.(-0.5 * sqeuclidean(targetdata, simdata)) / sqrt(2 * pi)
    return likelihood
end


# Get the distance between simulated and real data (read from files) from copy numbers
function compute_lnl_nsample(params, nsample)
  #println(params)
  mu1 = params[1]
  mu2 = params[2]
  birth_rate = params[3]
  death_rate = 0
  skip = 3

  sum_likelihood = 0.0
  for i in 1:nsample
      seed=rand(UInt32, 1)[1]

      # call simulation program to generate data
      cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --only_mut $only_mut`

      res = chomp(read(cmdSim, String))
      arr = split(res,"\t")
      simdata = map(x->parse(Float64,x), arr)
      simdata = simdata[1:4]

      # avoid 0 when taking log
      if only_mut == 1
          simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]), log(1+simdata[4]))
      else
          simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]), log(simdata[4]))
      end

      #println(simdata)
      sum_likelihood += compute_surrogate_likelihood(simdata, targetdata)
  end

  lnl = log(sum_likelihood / nsample)

  return lnl
end


function compute_dev_param(params, nsample)
    dev = -2 * compute_lnl_nsample(params, nsample)
    return dev
end


function compute_devbar(smcres, msample, nsample)
    sum_likelihood = 0.0
    resampled = wsample(1:nrow(smcres), smcres.weights, msample)

    for i in 1:msample
        sel = smcres[resampled[i],:]
        params = [sel.parameter1, sel.parameter2, sel.parameter3]
        #println(params)
        sum_likelihood += compute_lnl_nsample(params, nsample)
    end

    devbar = -2 * sum_likelihood / msample
    return devbar
end


function compute_dic(devbar, dev_param)
    dic = 2 * devbar - dev_param
    return dic
end

dev_param = compute_dev_param(params, nsample)
devbar = compute_devbar(smcres, msample, nsample)
dic = compute_dic(devbar, dev_param)

println(dic)
