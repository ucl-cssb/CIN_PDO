# sample commands to run ABC SMC on real data


########## for datasets without trees #############
epsilon=0.015
brate=0.5   # fix birth rate, negligible affect on inferrence of mutation rate

julia abc/runABC_SMC.jl dataset_20190310 72 $epsilon $brate &
julia abc/runABC_SMC.jl dataset_20190331 67 $epsilon $brate &
julia abc/runABC_SMC.jl dataset_04022020 37 $epsilon $brate &



########## for datasets with trees #############
s3_weight=2
only_mut=0
epsilon=0.2

julia runABC_SMC-joint-logdist-withbratio.jl dataset_20181127 24 $epsilon $s3_weight $only_mut
julia runABC_SMC-joint-withsel-logdist-withbratio.jl dataset_20181127 24 $epsilon $s3_weight $only_mut

julia runABC_SMC-joint-logdist-withbratio.jl dataset_20190409_1 25 $epsilon $s3_weight $only_mut
julia runABC_SMC-joint-withsel-logdist-withbratio.jl dataset_20190409_1 25 $epsilon $s3_weight $only_mut

julia runABC_SMC-joint-logdist-withbratio.jl dataset_20200304 19 $epsilon $s3_weight $only_mut
julia runABC_SMC-joint-withsel-logdist-withbratio.jl dataset_20200304 19 $epsilon $s3_weight $only_mut

# for a short tree in which branch length ratio is impossible to compute
epsilon=0.1
julia runABC_SMC-joint-logdist.jl dataset_20190409_2 13 $epsilon $s3_weight
julia runABC_SMC-joint-withsel-logdist.jl dataset_20190409_2 13 $epsilon $s3_weight
