using ApproxBayes
using Distributions
using DelimitedFiles
using Distances

# Run ABC assuming neutral model on copy number data only, with two summary statistics

# export JULIA_NUM_THREADS=32

# read from file
dataset=ARGS[1]
ncell=ARGS[2]
epsilon=ARGS[3]
birth_rate = ARGS[4]

# dataset="dataset_04022020"
# ncell=37
# epsilon=0.01
# birth_rate = 0.5

println(dataset)
println(ncell)


# change directory accordingly
dir = "../../"
progSim=dir * "/code/simpdo"
simdir=dir * "/data/abcsmc"
freal = dir * "/data/4trees/" * dataset * "_rvec.txt"


targetdata = readdlm(freal)[:, 1]
println(targetdata)


# Get the distance between simulated and real data (read from files) from copy numbers
function tumourABCneutralcn(params, constants, targetdata)
  # call simulation program to generate data
  # print(params)
  mu1=params[1]
  mu2=params[2]
  death_rate = 0

  ncell=constants[1]

  seed=rand(UInt32, 1)[1]
  suffix=string(seed) * "_" * dataset * "_n" * string(ncell)

  skip = 0

  cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --suffix $suffix`
  res = chomp(read(cmdSim, String))
  arr = split(res,"\t")
  simdata = map(x->parse(Float64,x), arr)
  simdata = simdata[1:2]

  # r = euclidean(targetdata, simdata)
  r = evaluate(Euclidean(1e-16), targetdata, simdata)

  return r, 1
end


setup = ABCSMC(tumourABCneutralcn, #simulation function
    2, # number of parameters
    parse(Float64, epsilon), # target ϵ
    Prior([Uniform(0, 1), Uniform(0, 1)]), #Prior for each of the parameters
    maxiterations = 9*10^6, #Maximum number of iterations before the algorithm terminates
    nparticles = 500,
    α = 0.3,
    ϵ1 = 10000.0,
    convergence = 0.05,
    constants = [ncell],
  )

smc = runabc(setup, targetdata, verbose = true, progress = true, parallel = false)

fout="smc-cna-" * dataset * "_n" * string(ncell) * "_epsilon" * string(epsilon)* "_brate" * string(birth_rate)
writeoutput(smc, dir=simdir, file=fout)
