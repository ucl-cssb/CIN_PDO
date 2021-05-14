using ApproxBayes
using Distributions
using DelimitedFiles
using Distances

# Run ABC assuming neutral model, with four summary statistics (including branch length ratios when available)

# export JULIA_NUM_THREADS=32

# read from file
dataset = ARGS[1]
ncell = ARGS[2]
epsilon = ARGS[3]
s3_weight = parse(Float64, ARGS[4])
only_mut = parse(UInt32, ARGS[5])

# dataset="dataset_20190409_1"
# ncell=25

# dataset="dataset_20190409_2"
# ncell=13

println(dataset)
println(ncell)

# epsilon = 0.1
min_brate = 0.2
max_brate = 2


# change directory accordingly
dir = "../../"
progSim=dir * "/code/simpdo"
frealCnvec = dir * "/data/4trees/" * dataset * "_rvec.txt"
frealTree = dir * "/data/4trees/" * dataset * "_divtree.txt"
frealTree_ratio = dir * "/data/4trees/" * dataset * "_divtree_blenratio.txt"

simdir=dir * "/data/abcsmc"


if only_mut == 1
    freal_sstat = dir * "/data/4trees/" * dataset * "_sstat4mut.txt"
else
    freal_sstat = dir * "/data/4trees/" * dataset * "_sstat4.txt"
end


fout="smc-joint-" * dataset * "_n" * string(ncell) * "_epsilon" * string(epsilon) * "_logs3_weight" * string(s3_weight) * "_logs4_mut" * string(only_mut)
if only_mut == 1
    fout = fout * "_s4plus1"
end

# if isfile(freal_sstat)
#     targetdata = readdlm(freal_sstat)
# else

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
    writedlm(freal_sstat, transpose(targetdata))
else
    exit()
end
# end

println(targetdata)



# Get the distance between simulated and real data (read from files) from copy numbers
function tumourABCneutralcn(params, constants, targetdata)
  # call simulation program to generate data
  # print(params)
  mu1=params[1]
  mu2=params[2]
  birth_rate = params[3]
  death_rate = 0

  ncell=constants[1]

  seed=rand(UInt32, 1)[1]

  # skip firsth 3 branches due to missing time measurements
  skip = 3

  cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --only_mut $only_mut`
  res = chomp(read(cmdSim, String))
  arr = split(res,"\t")
  simdata = map(x->parse(Float64,x), arr)
  simdata = simdata[1:4]

  if only_mut == 1
      simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]), log(1+simdata[4]))
  else
      simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]), log(simdata[4]))
  end

  r = euclidean(targetdata, simdata)
  # r = evaluate(Euclidean(1e-16), log.(targetdata), log.(simdata))

  return r, 1
end


setup = ABCSMC(tumourABCneutralcn, #simulation function
    3, # number of parameters
    parse(Float64, epsilon), # target ϵ
    Prior([Uniform(0, 1), Uniform(0, 1), Uniform(min_brate, max_brate)]), #Prior for each of the parameters
    maxiterations = 10*10^6, #Maximum number of iterations before the algorithm terminates
    nparticles = 500,
    α = 0.3,
    ϵ1 = 10000.0,
    convergence = 0.05,
    constants = [ncell],
  )


if only_mut == 1
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]), log(1+targetdata[4]))
else
    targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]), log(targetdata[4]))
end

smc = runabc(setup, targetdata, verbose = true, progress = true, parallel = true)

writeoutput(smc, dir=simdir, file=fout)
