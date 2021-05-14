using ApproxBayes
using Distributions
using DelimitedFiles
using Distances

# Run ABC assuming neutral model, with three summary statistics

# export JULIA_NUM_THREADS=32

# read from file
dataset = ARGS[1]
ncell = ARGS[2]
epsilon = ARGS[3]
s3_weight = parse(Float64, ARGS[4])

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
simdir=dir * "/data/abcsmc"
frealCnvec = dir * "/data/4trees/" * dataset * "_rvec.txt"
frealTree = dir * "/data/4trees/" * dataset * "_divtree.txt"
freal_sstat = dir * "/data/4trees/" * dataset * "_sstat3.txt"

fout="smc-joint-" * dataset * "_n" * string(ncell) * "_epsilon" * string(epsilon) * "_logs3_weight" * string(s3_weight)


targetdataCnvec = readdlm(frealCnvec)[:, 1]

# skip header and the first three cell divisions
targetTree = readdlm(frealTree, skipstart=7)[:, 3]
tblen = sum(targetTree)*3/1440/2
push!(targetdataCnvec, tblen)
targetdata = targetdataCnvec
println(targetdata)


writedlm(freal_sstat, transpose(targetdata))


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
  suffix=string(seed) * "_" * dataset * "_n" * string(ncell)

  skip = 3

  cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --suffix $suffix`
  res = chomp(read(cmdSim, String))
  arr = split(res,"\t")
  simdata = map(x->parse(Float64,x), arr)
  simdata = simdata[1:3]

  # simdata = vcat(log(simdata[1] + 1), log(simdata[2] + 1), log(simdata[3]))
  # simdata[simdata.=0] .= SMALL_VAL
  # simdata = vcat(log(simdata[1]), log(simdata[2]), log(simdata[3]))
  simdata = vcat((simdata[1]), (simdata[2]), s3_weight * log(simdata[3]))

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

targetdata = vcat((targetdata[1]), (targetdata[2]), s3_weight * log(targetdata[3]))

smc = runabc(setup, targetdata, verbose = true, progress = true, parallel = true)

writeoutput(smc, dir=simdir, file=fout)
