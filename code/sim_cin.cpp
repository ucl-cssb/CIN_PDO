#include <cstdlib>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "cell.hpp"


using namespace std;


typedef map<pair<int, int>, int> pcn;   // copy number at a position


// Read posterior distributions of parameters from file
// parameter1	parameter2	distance	weights
void read_params_posterior(string fname, vector<double>& param1, vector<double>& param2, vector<double>& param3, vector<double>& weights, int verbose = 0){
   ifstream fin(fname);
   string line;
   // int i = 0;
   while (getline(fin, line)) {
       // i++;
       // cout << "line " << i << endl;
       // skip comment lines
       if(!isdigit(line.c_str()[0])){
           continue;
       }
       // Split line into tab-separated parts
       vector<string> parts;
       // parameter1	parameter2	distance	weights
       boost::split(parts, line, boost::is_any_of("\t"));
       // cout << parts.size() << "\t" << parts[0] << "\t" << parts[1] << "\t" << parts[2] << endl;

       param1.push_back(stod(parts[0]));
       param2.push_back(stod(parts[1]));
       if(parts.size()>4){
           param3.push_back(stod(parts[2]));    // birth rates
           weights.push_back(stod(parts[4]));
       }else{
           weights.push_back(stod(parts[3]));
       }
   }
   fin.close();

   // Test reading input
   if(verbose > 0){
       cout << "Reading file " << fname << endl;
       cout << param1.size() << "\t" << param2.size() << "\t" << param3.size() << "\t" << weights.size() << "\n";
   }
   // for(int i = 0; i < param1.size(); i++){
   //     cout << param1[i] << "\t" << param2[i] << "\t" << weights[i] << "\n";
   // }
}


void get_params_file(string dir_param, vector<string>& filenames, int verbose = 0){
    for (boost::filesystem::directory_iterator itr(dir_param); itr!=boost::filesystem::directory_iterator(); ++itr)
    {
        string fname = itr->path().filename().string();
        // cout << itr->path().filename() << ' '; // display filename only
        // if (is_regular_file(itr->status())) cout << " [" << file_size(itr->path()) << ']';
        // cout << '\n';
        std::string prefix = "smc-";
        if(fname.compare(0, prefix.size(), prefix) == 0){
            filenames.push_back(itr->path().string());
        }
    }
    if(verbose > 0){
        cout << "Parameter files: " << endl;
        for(int i = 0; i < filenames.size(); i++){
            cout << filenames[i] << endl;
        }
        cout << endl;
    }
}


// Simulate the growth of one cell into a single organoid (clone)
Clone run_sim(double birth_rate, double death_rate, int Nend, double arm_prob, double chr_prob, double multi_prob, int skip, unsigned long rseed, string outdir, int genotype_diff = 0, double chr_weight = 1.0, int chr_sel=0, int use_ablen=0, string suffix="", string prefix="", int model = 0, string fitness_vals="", int num_subclone = 0, double tmin = 0, double tmax = 0, double min_clone_freq = 0, double max_clone_freq = 0, int num_clonal_mutation = 0, string file_cmut = "", int verbose = 0){
    Clone clone;

    double lambda = birth_rate - death_rate;
    double tend = log(Nend)/(lambda);

    double mutation_rate = arm_prob + chr_prob + multi_prob;

    vector<double> fitness;
    if(fitness_vals!=""){
        stringstream ss(fitness_vals);
        int num;   ss >> num;
        for(int i = 0; i < num; i++){
            double d1;
            ss >> d1;
            fitness.push_back(d1);
        }
    }

    vector<double> time_occur;

    if(verbose > 0){
        cout << "\nSimulating tumor growth" << endl;
        cout << "Random seed: " << rseed << endl;
        cout << "Model of evolution: " << model << endl;
        cout << "Initial Net growth rate: " << lambda << endl;
        if(fitness_vals != ""){
          cout << "Fitness values: " << fitness_vals << endl;
        }
        cout << "Mutation rate: " << mutation_rate << endl;
        cout << "\tChr-level CNA rate: " << chr_prob << endl;
        cout << "\tArm-level CNA rate: " << arm_prob << endl;
        cout << "\tMultipolar spindle (hopeful monster) rate: " << multi_prob << endl;
        cout << "\tEffective mutation rate (μ/β): " << mutation_rate/((birth_rate - death_rate)/birth_rate) << endl;
        cout << "Estimated simulation finish time (tumor doublings): " << tend << endl;
    }

    if(tmax > tend){
        tmax = tend;
    }

    if(mutation_rate > 0){
        if(verbose == 1) cout << "\nSimulating CIN" << endl;
        clone.grow_with_cnv(num_subclone, num_clonal_mutation, model, fitness, time_occur, Nend, birth_rate, death_rate, mutation_rate, arm_prob, chr_prob, multi_prob, genotype_diff, chr_weight, chr_sel, file_cmut, verbose);
    }
    else{
        clone.grow(num_subclone, num_clonal_mutation, fitness, time_occur, Nend, birth_rate, death_rate, mutation_rate, verbose);
    }

    // not output anything when simulating bulk samples
    if(verbose < 0) return clone;

    if(verbose > 0){
        const char* path = outdir.c_str();
        boost::filesystem::path dir(path);
        if(boost::filesystem::create_directory(dir))
        {
            if(verbose == 1) cerr << "Directory Created: " << outdir <<endl;
        }
        if(suffix==""){
            string sep = "-";
            suffix = sep + to_string(num_subclone) + sep + to_string(Nend) + sep + to_string(int(birth_rate*10)) + sep + to_string(int(death_rate*10)) + sep + to_string(int(mutation_rate)) + sep + to_string(num_clonal_mutation) + sep + prefix + sep + to_string(int(tmin)) + sep + to_string(int(rseed));
        }
        string outfile = "";

        // outfile = outdir + "curr_vaf" + suffix;
        // cout << "Computing VAF" << endl;
        // map<int, double> mut_freq = clone.get_allele_freq();
        // clone.print_map(mut_freq, outfile);

        outfile = outdir + "summary" + suffix;
        cout << "Printing summary" << endl;
        clone.print_summary(outfile);

        outfile = outdir + "end_cells_cn" + suffix + ".txt";
        clone.print_obs_cn(clone.curr_cells, outfile, verbose);

        // outfile = outdir + "end_cells_cn" + suffix + "_rvec.txt";
        // clone.write_avg_reciprocal_cn(clone.curr_cells, outfile, false);

        // // The proportion of chr-level and arm-level (p,q) events for each chromosome
        // outfile = outdir + "end_cells_cn" + suffix + "_prop.txt";
        // clone.write_prop(clone.curr_cells, outfile);
        //
        // outfile = outdir + "end_cells_cn" + suffix + "_cnvec.txt";
        // clone.write_avg_cn(clone.curr_cells, outfile);
        //
        // outfile = outdir + "end_cells_cn" + suffix + "_uvec.txt";
        // clone.write_avg_uniq_cn(clone.curr_cells, outfile);
        //
        // outfile = outdir + "end_cells_cn" + suffix + "_avec.txt";
        // clone.write_aggregated_uniq_cn(clone.curr_cells, outfile);
        //
        // outfile = outdir + "end_cells_cn" + suffix + "_lvec.txt";
        // clone.write_relative_reciprocal_cn(clone.curr_cells, outfile, false);

        // Output total number of de novo events in the population
        // outfile = outdir + "summary_stats" + suffix + ".txt";
        // clone.write_summary_stats(sum_stats, outfile);

        if(verbose > 1){
            cout << "Printing out informaton of all cells in a tumor clone" << endl;
            // cout << "file name" << outfile << endl;
            clone.print_lineage(clone.cells, outdir, suffix);

            outfile = outdir + "all_cells_tree" + suffix + ".nwk";
            clone.write_newick(clone.cells, outfile);

            outfile = outdir + "all_cells_tree" + suffix + ".txt";
            clone.write_tree(clone.cells, outfile);
        }
    }

    // vector<double> sum_stats;
    cout.precision(9);
    vector<int> avg_cn(2, 0);
    clone.get_avg_reciprocal_cn(clone.curr_cells, avg_cn, false);
    for(int i = 0; i < avg_cn.size(); i++){
        double ac = (double) avg_cn[i] / clone.curr_cells.size();
        // sum_stats.push_back(ac);
        cout << ac << "\t";
    }

    if(use_ablen == 1){     // Output all the tree lenghts
        vector<double> tratio;
        if(fitness_vals != ""){
            // clone.get_treelen_ratios(tratio, clone.cells, skip);
            // cout << "\t" << tratio[0] << "\t" << tratio[1];
            vector<double> blens;
            clone.get_treelen_vec(blens, clone.cells, skip);
            cout << blens[0];
            for(int i = 1; i < blens.size(); i++){
                cout << "\t" << blens[i];
            }
            // cout << endl;
            cout << "\t";
        }
    }else{
        double tlen = clone.get_treelen(clone.cells, skip);
        // sum_stats.push_back(tlen);
        // cout << tlen << endl;
        cout << tlen << "\t";
    }

    int ndenovo = 0;
    pcn clone_cnp;      // Group CNPs by arms
    pcn clone_cnc;
    for(auto cell : clone.curr_cells) {
        cell.set_obs_cn();
        for(auto cp : cell.obs_cn_profile){
            pair<int, int> pos = cp.first;
            if(cp.second != 0){
                ndenovo++;

                if(pos.second==0){  // divide whole chromosome to two arms
                    pair<int, int> pos1(pos.first, 1);
                    clone_cnp[pos1] += abs(cp.second);
                    clone_cnc[pos1] += 1;
                    pair<int, int> pos2(pos.first, 2);
                    clone_cnp[pos2] += abs(cp.second);
                    clone_cnc[pos2] += 1;
                }else{
                    clone_cnp[pos] += abs(cp.second);
                    clone_cnc[pos]++;
                }
            }
        }
    }
    int naler_pos_l1sample = 0;
    double sum_avg_cn = 0;
    map<int, int> nalter_chr;   // count the number of altered position for each chromosome (at most 2)
    for(auto cp : clone_cnp){
        pair<int, int> pos = cp.first;
        assert(pos.second!=0);
        nalter_chr[pos.first]++;
        if(clone_cnp[pos] > 1) naler_pos_l1sample += 1;
        // double avg_cn = (double) clone_cnp[pos] / clone_cnc[pos];     // taking average over altered cells
        double avg_cn = (double) clone_cnp[pos] / clone.curr_cells.size();     // taking average over all cells in an organoid
        sum_avg_cn += abs(avg_cn);
    }

    int nalter = 0;
    for(auto nc : nalter_chr){
        assert(nc.second <= 2);
        nalter += nc.second;
    }

    // double prop_alter = (double) nalter / (NUM_CHR * 2);
    // Output summary statistics
    // cout << prop_alter << endl;
    cout << nalter << "\t" << sum_avg_cn  << "\t" << ndenovo << "\t" << naler_pos_l1sample << endl;

    return clone;

}



// set bulk parameters for one run according to random choice
void set_bulk_params_random(vector<double>& params,  double min_birth_rate, double max_birth_rate, double min_arm_prob, double max_arm_prob, double min_chr_prob, double max_chr_prob, double min_mut_rate, double max_mut_rate, double avg_mut_rate, double max_multi_prob) {
    // double mut_rate = runiform(r, min_mut_rate, max_mut_rate);
    double mut_rate_init = gsl_ran_exponential(r, avg_mut_rate);
    double mut_rate = mut_rate_init > max_mut_rate? max_mut_rate : mut_rate_init;

    // // select arm_prob since chr_prob is biased towards small values in real data sets
    // double max_arm_rate = mut_rate < max_arm_prob? mut_rate : max_arm_prob;
    // double min_arm_rate = max_arm_rate < min_arm_prob? max_arm_rate : min_arm_prob;
    // arm_prob = runiform(r, min_arm_rate, max_arm_rate);
    // // chr_prob = runiform(r, min_chr_prob, max_arm_prob);
    // double chr_prob_init = mut_rate - arm_prob;
    // chr_prob = chr_prob_init > max_chr_prob? max_chr_prob:chr_prob_init;

    // select chr_prob since chr_prob is often less than arm_prob in real data sets
    double max_chr_rate = mut_rate < max_chr_prob? mut_rate : max_chr_prob;
    double min_chr_rate = max_chr_rate < min_chr_prob? max_chr_rate : min_chr_prob;
    double chr_prob = runiform(r, min_chr_rate, max_chr_rate);
    // chr_prob = runiform(r, min_chr_prob, max_arm_prob);
    double arm_prob_init = mut_rate - chr_prob;
    double arm_prob = arm_prob_init > max_arm_prob? max_arm_prob : arm_prob_init;

    double multi_prob = 0;
    if(max_multi_prob > 0){
        double multi_prob = runiform(r, 0, max_multi_prob);
    }

    mut_rate = chr_prob + arm_prob + multi_prob;
    double birth_rate = runiform(r, min_birth_rate, max_birth_rate);

    params.push_back(chr_prob);
    params.push_back(arm_prob);
    params.push_back(multi_prob);
    params.push_back(mut_rate);
    params.push_back(birth_rate);
}

// set bulk parameters for one run according to weighted choice from posterior distributions
void set_bulk_params(vector<double>& params, const vector<string>& filenames,  double min_birth_rate, double max_birth_rate, double max_multi_prob, int verbose = 0){
    // randomly pick up a dataset to choose parameters from
    int i = runiform(r, 0, filenames.size());

    string fname = filenames[i];
    vector<double> param1;
    vector<double> param2;
    vector<double> param3;
    vector<double> weights;
    read_params_posterior(fname, param1, param2, param3, weights, verbose);

    double* arr_weights = &weights[0];
    gsl_ran_discrete_t*  dis = gsl_ran_discrete_preproc(weights.size(), arr_weights);
    int sel = gsl_ran_discrete(r, dis);

    double chr_prob, arm_prob, mut_rate, birth_rate, multi_prob = 0;
    chr_prob = param1[sel];
    arm_prob = param2[sel];

    if(param3.size() > 0){
        birth_rate = param3[sel];
    }else{
        birth_rate = runiform(r, min_birth_rate, max_birth_rate);
    }
    assert(birth_rate > 0);

    // Multipolar rate:
    // 1/24=0.0417  divisions in 20190409_1
    // 1/23=0.0435  divisions in 20181127
    // 1/18=0.0556  divisions in 20200403
    if(max_multi_prob > 0){
        // multi_prob = runiform(r, 0, max_multi_prob);
        if(fname.find("20190409_1")!=std::string::npos){
            multi_prob = gsl_ran_exponential(r, 0.0417);
        }else if(fname.find("20181127")!=std::string::npos){
            multi_prob = gsl_ran_exponential(r, 0.0435);
        }else if(fname.find("20200403")!=std::string::npos){
            multi_prob = gsl_ran_exponential(r, 0.0556);
        }else{
            multi_prob = 0;
        }
    }

    // Used for printing
    mut_rate = chr_prob + arm_prob + multi_prob;

    params.push_back(chr_prob);
    params.push_back(arm_prob);
    params.push_back(multi_prob);
    params.push_back(mut_rate);
    params.push_back(birth_rate);

    if(verbose > 0)
        cout << "Selected parameters: " << chr_prob << "\t" << arm_prob << "\t" << multi_prob << "\t" << mut_rate << "\t" << birth_rate << endl;
}


// set single organoid parameters for one run according to weighted choice from posterior distributions
void set_organoid_params(vector<double>& params, string filename, double& birth_rate, int verbose = 0){
    vector<double> param1;
    vector<double> param2;
    vector<double> param3;
    vector<double> weights;
    read_params_posterior(filename, param1, param2, param3, weights, verbose);

    double* arr_weights = &weights[0];
    gsl_ran_discrete_t*  dis = gsl_ran_discrete_preproc(weights.size(), arr_weights);
    int sel = gsl_ran_discrete(r, dis);

    double chr_prob, arm_prob;
    chr_prob = param1[sel];
    arm_prob = param2[sel];

    if(param3.size() > 0){
        birth_rate = param3[sel];
    }
    assert(birth_rate > 0);

    // mut_rate = chr_prob + arm_prob + multi_prob;

    params.push_back(chr_prob);
    params.push_back(arm_prob);
    // params.push_back(mut_rate);
    params.push_back(birth_rate);

    if(verbose > 0)
        cout << "Selected parameters: " << chr_prob << "\t" << arm_prob << "\t" << birth_rate << endl;
}


void get_bulk_stat(pcn& bulk_cnp, pcn& bulk_cnc, int nsample, int ndenovo, int verbose = 0){
    // Take the average of each bulk sample to get the CN profiles of all bulk samples
    double sum_avg_cn = 0;
    // double total_cn = 0;
    map<int, int> nalter_chr;   // count the number of altered position for each chromosome (at most 2)
    int naler_pos_l1sample = 0;    // count the number of  altered position for each chromosome which appear in more than one samples
    // cout << "Bulk CNPs: " << endl;
    assert(bulk_cnp.size() <= 2 * NUM_CHR);
    for(auto cp : bulk_cnp){
        pair<int, int> pos = cp.first;
        assert(pos.second!=0);
        nalter_chr[pos.first]++;
        if(bulk_cnc[pos]>1) naler_pos_l1sample++;
        // double avg_cn = (double) bulk_cnp[pos] / bulk_cnc[pos];     // taking average over all samples
        double avg_cn = (double) bulk_cnp[pos] / nsample;     // taking average over all samples
        // if(verbose>0) cout << pos.first + 1 << "\t" << pos.second << "\t" << cp.second << "\t" << avg_cn << "\t" << bulk_cnc[pos] << endl;
        sum_avg_cn += abs(avg_cn);
        // total_cn += cp.second;
    }

    int nalter_bulk = 0;
    for(auto nc : nalter_chr){
        assert(nc.second <= 2);
        nalter_bulk += nc.second;
    }

    // Compute the proportion of altered positions
    if(verbose>0){
        cout << "Number of CNP in bulk sample is " << bulk_cnp.size() << endl;
        cout << "Number of altered CNP by arm in bulk sample is " << nalter_bulk << endl;
    }
    // double prop_alter = (double) nalter_bulk / (NUM_CHR * 2);
    cout << nalter_bulk << "\t" << sum_avg_cn << "\t" << ndenovo  << "\t" << naler_pos_l1sample << endl;
}


void sim_bulk(int nsample, double min_birth_rate, double max_birth_rate, double death_rate, int min_ncell, int max_ncell, double min_arm_prob, double max_arm_prob, double min_chr_prob, double max_chr_prob, double min_mut_rate, double max_mut_rate, double avg_mut_rate, double max_multi_prob, int skip, unsigned long rseed, string dir_param, string outdir, int genotype_diff = 0, double chr_weight = 1.0, int chr_sel=0, int use_ablen=0, string suffix="", string prefix="", int model = 0, string fitness_vals="", int num_subclone = 0, double tmin = 0, double tmax = 0, double min_clone_freq = 0, double max_clone_freq = 0, int num_clonal_mutation = 0, string file_cmut = "", int verbose = 0){
    pcn bulk_cnp;   // used to compute the average of bulk CNPs
    pcn bulk_cnc;
    int ndenovo = 0;
    // upport for moving iostreams was added to GCC 5.1 (not work on cs cluster with default gcc)
    ofstream fcn(" "), fcn_rate(" ");
    vector<string> filenames;

    get_params_file(dir_param, filenames, verbose);

    // Suppress all the output when verbose = -1
    if(verbose >= 0){
        string outfile = outdir + "bulk_cn" + suffix + ".txt";
        fcn = ofstream(outfile);
        // Output the parameters used for each simulation for double checking
        string outfile_rate = outdir + "bulk_rate" + suffix + ".txt";
        fcn_rate = ofstream(outfile_rate);
    }

    string orig_suffix = suffix;
    double avg_bulk_birth_rate = 0, avg_bulk_mut_rate = 0, avg_bulk_chr_rate = 0, avg_bulk_arm_rate = 0, avg_bulk_mp_rate = 0;
    double avg_bulk_ncell = 0, avg_bulk_num_mut = 0, avg_bulk_first_mut_time = 0, bulk_first_mut_count = 0;

    for(int i = 0; i < nsample; i++){
        int ncell = runiform(r, min_ncell, max_ncell);

        vector<double> params;
        // set_bulk_params_random(params, min_birth_rate, max_birth_rate, min_arm_prob, max_arm_prob, min_chr_prob, max_chr_prob, min_mut_rate, max_mut_rate, avg_mut_rate, max_multi_prob);
        set_bulk_params(params, filenames, min_birth_rate, max_birth_rate, max_multi_prob, verbose);

        double chr_prob = params[0];
        double arm_prob = params[1];
        double multi_prob = params[2];
        double mut_rate = params[3];
        double birth_rate = params[4];

        avg_bulk_birth_rate += birth_rate;
        avg_bulk_mut_rate += mut_rate;
        avg_bulk_chr_rate += chr_prob;
        avg_bulk_arm_rate += arm_prob;
        avg_bulk_mp_rate += multi_prob;

        if(verbose>0){
            cout << "Simulating bulk sample " << i+1 << " with " << ncell << " cells" << endl;
            cout << "\tBirth rate: " << birth_rate << endl;
            cout << "\tMutation rate: " << mut_rate << endl;
            cout << "\tChr-level CNA rate: " << chr_prob << endl;
            cout << "\tArm-level CNA rate: " << arm_prob << endl;
            cout << "\tMultipolar spindle (hopeful monster) rate: " << multi_prob << endl;
            // cout << "\tDeath rate: " << death_rate << endl;
            // cout << "\tEffective mutation rate (μ/β): " << cell.mutation_rate / ((cell.birth_rate-cell.death_rate)/cell.birth_rate) << endl;
        }
        suffix = orig_suffix + "_bulk" + to_string(i+1);
        Clone s1 = run_sim(birth_rate, death_rate, ncell, arm_prob, chr_prob, multi_prob, skip, rseed + i, outdir, genotype_diff, chr_weight, chr_sel, use_ablen, suffix, prefix, model, fitness_vals, num_subclone, tmin, tmax, min_clone_freq, max_clone_freq, num_clonal_mutation, file_cmut, verbose);

        avg_bulk_ncell += ncell;
        avg_bulk_num_mut += s1.num_novel_mutation;
        avg_bulk_first_mut_time += s1.time_first_mut;
        if(s1.time_first_mut==1) bulk_first_mut_count += 1;
        if(verbose >= 0)
            fcn_rate << ncell << "\t" << s1.num_novel_mutation << "\t" << s1.time_first_mut << "\t" << birth_rate << "\t" << mut_rate << "\t" << chr_prob << "\t" << arm_prob << "\t" << multi_prob << endl;


        // by default, C++ creates an empty map with value for unknown key being 0
        pcn s1_cnp_split;
        pcn s1_cnc_split;
        pcn s1_cnp_merged;
        pcn s1_cnc_merged;

        // Use s1_cnp_split to merge events on one chromosome across cells in the organoid
        for(auto cell : s1.curr_cells) {
            cell.set_obs_cn();
            for(auto cp : cell.obs_cn_profile){
                pair<int, int> pos = cp.first;
                if(cp.second != 0){
                    if(pos.second==0){  // divide whole chromosome to two arms
                        pair<int, int> pos1(pos.first, 1);
                        s1_cnp_split[pos1] += cp.second;
                        s1_cnc_split[pos1] += 1;

                        pair<int, int> pos2(pos.first, 2);
                        s1_cnp_split[pos2] += cp.second;
                        s1_cnc_split[pos2] += 1;
                    }else{
                        s1_cnp_split[pos] += cp.second;
                        s1_cnc_split[pos] += 1;
                    }
                }
            }
        }
        assert(s1_cnp_split.size() <= 2 * NUM_CHR);

        // Take the average of each organoid sample to get the CN profiles of a bulk sample
        for(auto cp : s1_cnp_split){
            pair<int, int> pos = cp.first;
            assert(pos.second != 0);
            // double acn = (double) s1_cnp_split[pos] / s1_cnc_split[pos];
            double acn = (double) s1_cnp_split[pos] / ncell;        // The detected CN should be based on total number of cells
            // int racn = int(acn);
            // if(round_cn==1)
            int racn = round(acn);
            // if(i==0){
            //     cout << i + 1 << "\t" << pos.first + 1 << "\t" << pos.second << "\t" << s1_cnp_split[pos] << "\t" << s1_cnc_split[pos] << "\t" << acn << "\t" << racn << endl;
            // }
            if(racn == 0) continue;

            // chr-level event on one chr has been recorded
            pair<int, int> pos_chr(pos.first, 0);
            if(s1_cnp_merged.find(pos_chr) != s1_cnp_merged.end()){
                continue;
            }

            pair<int, int> pos2;
            if(pos.second==1){
                pos2 = make_pair(pos.first, 2);
            }else{
                assert(pos.second==2);
                pos2= make_pair(pos.first, 1);
            }
            if(s1_cnp_split.find(pos2) != s1_cnp_split.end()){
                double acn2 = (double) s1_cnp_split[pos2] / s1_cnc_split[pos2];
                // int racn2 = int(acn2);
                // if(round_cn==1){
                int racn2 = round(acn2);
                // }
                if(racn == racn2){
                    s1_cnp_merged[pos_chr] = racn;
                    bulk_cnp[pos] += abs(racn);
                    bulk_cnc[pos] += 1;
                    bulk_cnp[pos2] += abs(racn);
                    bulk_cnc[pos2] += 1;
                    continue;
                }
            }
            s1_cnp_merged[pos] = racn;
            bulk_cnp[pos] += abs(racn);
            bulk_cnc[pos] += 1;
        }

        int nalt = 0;
        for(auto cp : s1_cnp_merged){
            pair<int, int> pos = cp.first;
            int rmcn = s1_cnp_merged[pos];
            if(rmcn == 0) continue;
            nalt++;
            // Output the CN of each bulk sample
            if(verbose >= 0)
                fcn << i + 1 << "\t" << pos.first + 1 << "\t" << pos.second << "\t" << rmcn << endl;
        }

        ndenovo += nalt;

        if(verbose>0) cout << "Number of altered CNP in sample " << i+1 << " is " << nalt << endl;
    }  // Simulate a bulk sample with nsample bulk organoid data

    cout << avg_bulk_ncell / nsample << "\t" << avg_bulk_num_mut / nsample << "\t" << avg_bulk_first_mut_time / nsample << "\t" << bulk_first_mut_count << "\t" << avg_bulk_birth_rate / nsample << "\t" << avg_bulk_mut_rate / nsample << "\t" << avg_bulk_chr_rate / nsample << "\t" << avg_bulk_arm_rate / nsample << "\t" << avg_bulk_mp_rate / nsample << endl;

    get_bulk_stat(bulk_cnp, bulk_cnc, nsample, ndenovo, verbose);

    if(verbose >= 0){
        fcn.close();
        fcn_rate.close();
    }
}


// TODO: pass parameters
int main(int argc, char const *argv[]) {
    Clone clone;

    int mode;
    int nsample;
    int min_ncell, max_ncell;
    double min_mut_rate, max_mut_rate, avg_mut_rate;
    double min_arm_prob, max_arm_prob;
    double min_chr_prob, max_chr_prob;
    double max_multi_prob;
    double min_birth_rate, max_birth_rate;
    string dir_param;

    int Nend;
    double birth_rate, death_rate;
    double mutation_rate, arm_prob, chr_prob, multi_prob;
    int multi_nchr;
    string file_param;

    int model;
    string fitness_vals;
    int num_subclone;
    double tmin, tmax;
    // Mutation rate per division per genome
    double min_clone_freq, max_clone_freq;
    int num_clonal_mutation;

    int chr_sel; // level of selection
    int genotype_diff;
    double chr_weight;
    int use_ablen;
    int skip; // number of cells to skip when computing division time
    string outdir, prefix, suffix; // output

    unsigned long seed;

    string tree_file, file_cmut;

    int verbose;

    namespace po = boost::program_options;

    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help,h", "produce help message")
      ;

    po::options_description required("Required parameters");
    required.add_options()
      ("odir,o", po::value<string>(&outdir)->required()->default_value("./"), "output directory")
       ;

    po::options_description optional("Optional parameters");
    optional.add_options()
      ("mode", po::value<int>(&mode)->default_value(0), "mode of simulation. 0: single cell data; 1: bulk data")

      ("birth_rate", po::value<double>(&birth_rate)->default_value(1), "birth rate")
      ("death_rate", po::value<double>(&death_rate)->default_value(0), "death rate")
      ("Nend,e", po::value<int>(&Nend)->default_value(100), "size of final cell populations")

      ("nsample", po::value<int>(&nsample)->default_value(100), "number of bulk samples")
      ("min_ncell", po::value<int>(&min_ncell)->default_value(100), "minimum number of cells in a bulk sample")
      ("max_ncell", po::value<int>(&max_ncell)->default_value(100), "maximum number of cells in a bulk sample")
      ("min_mut_rate", po::value<double>(&min_mut_rate)->default_value(0.02), "minimum CNA rate")
      ("max_mut_rate", po::value<double>(&max_mut_rate)->default_value(0.2), "maximum CNA rate")
      ("avg_mut_rate", po::value<double>(&avg_mut_rate)->default_value(0.1), "average CNA rate")
      ("min_chr_prob", po::value<double>(&min_chr_prob)->default_value(0.01), "minimum chr-level CNA rate")
      ("max_chr_prob", po::value<double>(&max_chr_prob)->default_value(0.2), "maximum chr-level CNA rate")
      ("min_arm_prob", po::value<double>(&min_arm_prob)->default_value(0.02), "minimum arm-level CNA rate")
      ("max_arm_prob", po::value<double>(&max_arm_prob)->default_value(0.2), "maximum arm-level CNA rate")
      ("max_multi_prob", po::value<double>(&max_multi_prob)->default_value(0.06), "maximum multiple simultaneous chr-level CNA rate (hopeful monster)")
      ("min_birth_rate", po::value<double>(&min_birth_rate)->default_value(0.2), "minimum cell division rate")
      ("max_birth_rate", po::value<double>(&max_birth_rate)->default_value(0.5), "maximum cell division rate")
      ("dir_param", po::value<string>(&dir_param)->default_value(""), "directory containing posterior distributions of inferred parameters")

      ("model", po::value<int>(&model)->default_value(0), "model of evolution. 0: neutral; 1: gradual; 2: punctuated; 3: positive")
      ("fitness_vals", po::value<string>(&fitness_vals)->default_value(""), "fitness values of mutatants")
      ("chr_sel", po::value<int>(&chr_sel)->default_value(0), "level of selection. 0: selection on all CNAs; 1: selection on chr-level CNAs")
      ("genotype_diff", po::value<int>(&genotype_diff)->default_value(0), "whether or not to use genotype difference (L1 distance) in simulating selection. 0: no; 1: yes")
      ("chr_weight", po::value<double>(&chr_weight)->default_value(1.0), "relative fitness of chr-level CNAs")

      ("use_ablen", po::value<int>(&use_ablen)->default_value(0), "tree summary statistics. 0: half tree length; 1: vector of all branch lengths")

      ("num_subclone,n", po::value<int>(&num_subclone)->default_value(0), "number of subclones to simulate")
      ("min_clone_freq", po::value<double>(&min_clone_freq)->default_value(0), "the minimal frequency of a subclone")
      ("max_clone_freq", po::value<double>(&max_clone_freq)->default_value(0), "the maximal frequency of a subclone")
      ("tmin", po::value<double>(&tmin)->default_value(0), "the earliest time that a subclone occurs")
      ("tmax", po::value<double>(&tmax)->default_value(0), "the latest time that a subclone occurs")

      ("num_clonal_mutation", po::value<int>(&num_clonal_mutation)->default_value(0), "number of clonal mutations")
      ("file_cmut", po::value<string>(&file_cmut)->default_value(""), "the TSV file which contains clonal CNVs, with three columns (chr, arm, cn)")

      ("arm_prob, a", po::value<double>(&arm_prob)->default_value(0.5), "arm-level CNA rate per cell division")
      ("chr_prob, c", po::value<double>(&chr_prob)->default_value(0.5), "chr-level CNA rate per cell division")
      ("multi_prob, m", po::value<double>(&multi_prob)->default_value(0), "multiple simultaneous chr-level CNA rate per cell division")
      ("multi_nchr", po::value<int>(&multi_nchr)->default_value(16), "number of chromosomes affected by multipolar division")
      ("file_param", po::value<string>(&file_param)->default_value(""), "the file containing posterior distributions of inferred parameters")

      ("prefix,p", po::value<string>(&prefix)->default_value("run1"), "prefix of output file (it will be sim-data-N if not specified")
      ("suffix,s", po::value<string>(&suffix)->default_value(""), "suffix of output file")
      ("skip", po::value<int>(&skip)->default_value(3), "number of cells to skip when computing division time")

      ("seed", po::value<unsigned long>(&seed)->default_value(0), "seed used for generating random numbers")
      ("verbose", po::value<int>(&verbose)->default_value(0), "verbose level (0: default, 1: print information of final cells; 2: print information of all cells)")
      ;

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(required).add(optional);
    po::variables_map vm;

    try {
        po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
        if(vm.count("help")){
            cout << cmdline_options << endl;
            return 1;
        }
        if(vm.count("version")){
            cout << "sim_cin [version 0.1], a program to simulate copy number variations along a cell division tree" << endl;
            return 1;
        }
        po::notify(vm);
    }
    catch (const exception& e) {
          cerr << e.what() << endl;
          return 1;
    }


    unsigned long rseed = setup_rng(seed);

    if(mode == 0){
        if(file_param!=""){
            vector<double> params;
            set_organoid_params(params, file_param, birth_rate, verbose);

            chr_prob = params[0];
            arm_prob = params[1];
            birth_rate = params[2];
        }
        run_sim(birth_rate, death_rate, Nend, arm_prob, chr_prob, multi_prob, skip, rseed, outdir, genotype_diff, chr_weight, chr_sel, use_ablen, suffix, prefix, model, fitness_vals, num_subclone, tmin, tmax, min_clone_freq, max_clone_freq, num_clonal_mutation, file_cmut, verbose);
    }else{
        sim_bulk(nsample, min_birth_rate, max_birth_rate, death_rate, min_ncell, max_ncell, min_arm_prob, max_arm_prob, min_chr_prob, max_chr_prob, min_mut_rate, max_mut_rate, avg_mut_rate, max_multi_prob, skip, rseed, dir_param, outdir, genotype_diff, chr_weight, chr_sel, use_ablen, suffix, prefix, model, fitness_vals, num_subclone, tmin, tmax, min_clone_freq, max_clone_freq, num_clonal_mutation, file_cmut, verbose);
    }

    return 0;
}
