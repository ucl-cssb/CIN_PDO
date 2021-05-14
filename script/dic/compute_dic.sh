# sample commands to compute DIC on real data

# type="neutral"
type="withsel"

s3_weight=2.0


only_mut=0
epsilon=0.2

dataset=dataset_20200304
ncell=19
bash sub_dic_real_logdist_withbratio.sh $dataset $ncell $epsilon $type $s3_weight $only_mut

dataset=dataset_20190409_1
ncell=25
bash sub_dic_real_logdist_withbratio.sh $dataset $ncell $epsilon $type $s3_weight $only_mut

dataset=dataset_20181127
ncell=24
bash sub_dic_real_logdist_withbratio.sh $dataset $ncell $epsilon $type $s3_weight $only_mut


# for a short tree in which branch length ratio is impossible to compute
epsilon=0.1

dataset=dataset_20190409_2
ncell=13
bash sub_dic_real_logdist.sh $dataset $ncell $epsilon $type $s3_weight
