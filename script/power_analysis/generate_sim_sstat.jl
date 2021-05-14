using Distributions
using DelimitedFiles
using DataFrames


# Generate n simulations under different models and record their summary statistics

ncell = parse(UInt32, ARGS[1])
model = parse(UInt32, ARGS[2])
nsims = parse(UInt32, ARGS[3])
only_mut = 0
if size(ARGS, 1) > 3
    only_mut = parse(UInt32, ARGS[4])
end
# ncell = 20
# model = 1
# nsims = 2


# change directory accordingly
dir = "../../"
progSim = dir * "/code/simpdo"
simdir = dir * "/data/power_analysis/"


skip = 3

min_mu = 0
max_mu = 0.5
# mean_mu = 0.1
# sigma_mu = 0.01

min_brate = 0.3
max_brate = 0.5
# mean_brate = 0.4
# sigma_brate = 0.02

min_fitness = -1
max_fitness = 0
# mean_fitness = -0.5
# if size(ARGS, 1) > 3
#     mean_fitness = parse(UInt32, ARGS[4])
# end
# sigma_fitness = 0.5



function gen_simulations_neutral(fout, ncell, nsims, dim_stat)
    res_sim = zeros(Float64, nsims, dim_stat)
    for i in 1:nsims
        mu1 = rand(Uniform(0,max_mu), 1)[1]
        mu2 = rand(Uniform(0,max_mu), 1)[1]
        birth_rate = rand(Uniform(min_brate, max_brate), 1)[1]
        # mu1 = rand(truncated(Normal(mean_mu, sigma_mu), min_mu, max_mu), 1)[1]
        # mu2 = rand(truncated(Normal(mean_mu, sigma_mu), min_mu, max_mu), 1)[1]
        # birth_rate = rand(truncated(Normal(mean_brate, sigma_brate), min_brate, max_brate), 1)[1]
        # birth_rate = brate_sim
        death_rate = 0

        seed=rand(UInt32, 1)[1]

        cmdSim=`$progSim -o $dir -e $ncell --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --only_mut $only_mut`

        res = chomp(read(cmdSim, String))
        arr = split(res,"\t")
        simdata = map(x->parse(Float64,x), arr)

        params = [mu1, mu2, birth_rate, seed]
        res_sim[i,:] = vcat(params, simdata)
    end
    writedlm(fout, res_sim)
end


function gen_simulations_withsel(fout, ncell, nsims, dim_stat)
    res_sim = zeros(Float64, nsims, dim_stat)
    for i in 1:nsims
        mu1 = rand(Uniform(0,max_mu), 1)[1]
        mu2 = rand(Uniform(0,max_mu), 1)[1]
        birth_rate = rand(Uniform(min_brate, max_brate), 1)[1]
        # mu1 = rand(truncated(Normal(mean_mu, sigma_mu), min_mu, max_mu), 1)[1]
        # mu2 = rand(truncated(Normal(mean_mu, sigma_mu), min_mu, max_mu), 1)[1]
        # birth_rate = rand(truncated(Normal(mean_brate, sigma_brate), min_brate, max_brate), 1)[1]
        # birth_rate = brate_sim
        fitness = rand(Uniform(-1, 0), 1)[1]
        # fitness = rand(truncated(Normal(mean_fitness, sigma_fitness), min_fitness, max_fitness), 1)[1]
        death_rate = 0

        seed=rand(UInt32, 1)[1]

        fitval = "1 " * string(fitness)
        model = 1

        cmdSim=`$progSim -o $dir -e $ncell --model $model --fitness "$fitval" --chr_prob $mu1 --arm_prob $mu2 --birth_rate $birth_rate --death_rate $death_rate --skip $skip --seed $seed --only_mut $only_mut`

        res = chomp(read(cmdSim, String))
        arr = split(res,"\t")
        simdata = map(x->parse(Float64,x), arr)

        params = [mu1, mu2, birth_rate, fitness, seed]
        res_sim[i,:] = vcat(params, simdata)
    end
    writedlm(fout, res_sim)
end


# suffix = string(ncell) * "_n" * string(nsims) * "_tnorm.tsv"
# used for current figure
# suffix = string(ncell) * "_n" * string(nsims) * "_maxmu" * string(max_mu) * "_minb" * string(min_brate) * "_maxb" * string(max_brate) * ".tsv"
# adding average of branch ratios
if only_mut == 1
    suffix = string(ncell) * "_n" * string(nsims) * "_maxmu" * string(max_mu) * "_minb" * string(min_brate) * "_maxb" * string(max_brate) * "_withbratio_mut.tsv"
else
    suffix = string(ncell) * "_n" * string(nsims) * "_maxmu" * string(max_mu) * "_minb" * string(min_brate) * "_maxb" * string(max_brate) * "_witheblen.tsv"
end

# Simulate
if model == 0
    println("simulate data under neutral evolution")
    fout = simdir * "sim_data_neutral_ncell" * suffix
    # fout = simdir * "sim_data_neutral_ncell" *
    # dim_stat=11
    dim_stat=13
    gen_simulations_neutral(fout, ncell, nsims, dim_stat)
else
    println("simulate data with selection")
    fout = simdir * "sim_data_selection_ncell" * suffix
    # dim_stat=12
    dim_stat=14
    gen_simulations_withsel(fout, ncell, nsims, dim_stat)
end
