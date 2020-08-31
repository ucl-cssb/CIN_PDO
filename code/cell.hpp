#ifndef CELL_HPP
#define CELL_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <random>
#include <string>

#include <assert.h>

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "stats.hpp"


using namespace std;



const int NUM_CHR = 23;
const int NORM_PLOIDY = 2;

// Number of chromsomes affected in a multipolar event
const int MULTI_NCHR = 16;

class Mutation
{
public:
    int mut_ID;

    int type;   // Mutation type
    int size;   // 1 for SNV
    double time_occur;

    int cell_ID;
    int chr;  // chromosome on which the Mutation occurs
    int arm;
    int reciprocal;

    int start;  // start position
    int end;    // end position

    double vaf;     // for SNV
    int number;     // snv number

    Mutation(){
            mut_ID = 0;
            time_occur = 0;
            vaf = 0;
            number = 1;
    }


    Mutation(int mut_ID, double time_occur){
            this->mut_ID = mut_ID;
            this->time_occur = time_occur;
            this->vaf = 0;
            this->number = 1;
    }

    Mutation(int mut_ID, double time_occur, int chr, int arm, int type, int reciprocal){
            this->mut_ID = mut_ID;
            this->time_occur = time_occur;

            this->vaf = 0;
            this->number = 1;

            this->chr = chr;
            this->arm = arm;
            this->type = type;
            this->reciprocal = reciprocal;
    }

    Mutation(int mut_ID, double time_occur, int chr, int start, int end, int arm, int type, int reciprocal){
            this->mut_ID = mut_ID;
            this->time_occur = time_occur;
            this->vaf = 0;
            this->number = 1;

            this->chr = chr;
            this->start = start;
            this->end = end;
            this->arm = arm;
            this->type = type;
            this->reciprocal = reciprocal;
    }

    // ~Mutation() = default;
    // Mutation(const Mutation& other) = default;
    // Mutation(Mutation&& other) = default;
    // Mutation& operator=(const Mutation& other) = default;
    // Mutation& operator=(Mutation&& other) = default;
};


class Cell
{
public:
    int cell_ID;
    int parent_ID;
    int clone_ID;
    vector<int> daughters;

    double birth_rate;
    double death_rate;

    double mutation_rate;
    double arm_prob;
    double chr_prob;
    double multi_prob;
    int multi_nchr;

    double ploidy;  // change ploidy to reflect WGD
    int num_division;
    double time_occur;

    int flag;   // whether the cell is alive or not. 0:new, 1:divided, -1: death
    double fitness;

    double pos_x;
    double pos_y;
    double pos_z;

    int wgd;
    // copy number for chr, arm; relative so that it is safe to ignore clonal events
    map<pair<int, int>, int> cn_profile;        // only applied to each arm
    map<pair<int, int>, int> obs_cn_profile;    // grouped by chr-/arm- level
    vector<float> cn_all;   // copy number for all regions

    Cell *sibling = NULL;

    double chr_probs[NUM_CHR];
    double arm_probs[2] = {0.5, 0.5};
    //chr gain, chr loss, arm gain, arm loss, (arm) doubling
    vector<double> cna_type_probs{0.3, 0.3, 0.2, 0.2, 0};

    vector<Mutation> mutations;

    ~Cell() = default;
    Cell(const Cell& other) = default;
    Cell(Cell&& other) = default;
    Cell& operator=(const Cell& other) = default;
    Cell& operator=(Cell&& other) = default;

    Cell() {
            cell_ID = 0;
            parent_ID = 0;
            clone_ID = 0;

            birth_rate = log(2);
            death_rate = 0;

            mutation_rate = 0;
            arm_prob = 0;
            chr_prob = 0;
            multi_prob = 0;
            multi_nchr = MULTI_NCHR;

            ploidy = 2;
            num_division = 0;
            time_occur = 0;
            flag = 0;

            // parent = NULL;
            for(int i = 0; i < NUM_CHR; i++){
                chr_probs[i] = double(1/NUM_CHR);
                for(int j=0; j < 3; j++){
                    pair<int, int> pos(i,j);
                    cn_profile[pos] = 0;
                }
            }
    }

    Cell(int cell_ID, int parent_ID) {
            this->cell_ID = cell_ID;
            this->parent_ID = parent_ID;
            this->clone_ID = 0;

            this->birth_rate = log(2);
            this->death_rate = 0;

            this->mutation_rate = 0;
            this->arm_prob = 0;
            this->chr_prob = 0;
            this->multi_prob = 0;
            this->multi_nchr = MULTI_NCHR;

            this->ploidy = 2;
            this->num_division = 0;
            this->time_occur = 0;
            this->flag = 0;

            for(int i = 0; i < NUM_CHR; i++){
                chr_probs[i] = double(1/NUM_CHR);
                for(int j=0; j < 3; j++){
                    pair<int, int> pos(i,j);
                    cn_profile[pos] = 0;
                }
            }
    }

    Cell(int cell_ID, int parent_ID, double birth_rate, double death_rate, double mutation_rate, double ploidy, double time_occur){
            this->cell_ID = cell_ID;
            this->parent_ID = parent_ID;
            this->clone_ID = 0;

            this->birth_rate = birth_rate;
            this->death_rate = death_rate;

            this->mutation_rate = mutation_rate;
            this->arm_prob = 0;
            this->chr_prob = 0;
            this->multi_prob = 0;
            this->multi_nchr = MULTI_NCHR;

            this->ploidy = ploidy;
            this->time_occur = time_occur;
            this->flag = 0;

            for(int i = 0; i < NUM_CHR; i++){
                chr_probs[i] = double(1/NUM_CHR);
                for(int j=0; j < 3; j++){
                    pair<int, int> pos(i,j);
                    cn_profile[pos] = 0;
                }
            }
    }


    Cell* get_parent(vector<Cell>& cells){
        for(int i = 0; i < cells.size(); i++){
            Cell* cell = &cells[i];
            if(cell->cell_ID == parent_ID) return cell;
        }
        return NULL;
    }

    /*
    mut_ID -- the ID of last mutation
    generate a random number of mutations
    */
    int generate_mutations(double mutation_rate, int& mut_ID, double time_occur){
            // poisson_distribution<int> pois(mutation_rate);
            // int nu = pois(eng);
            int nu = gsl_ran_poisson(r, mutation_rate);
            // cout << "Generating " << nu << " mutations" << endl;
            for (int j=0; j < nu; j++) {
                    mut_ID += 1;
                    // cout << mut_ID << "\t" << time_occur << endl;
                    Mutation mut(mut_ID, time_occur);
                    this->mutations.push_back(mut);
            }
            // cout << this->mutations.size() << endl;
            return nu;
    }


    /*
    generate a random number of CNAs based on a single mutation rate and type probabilities
    */
    int generate_CNV_bytype(int& mut_ID, double time_occur, int verbose = 0){
        assert(mutation_rate > 0);
        // poisson_distribution<int> pois(mutation_rate);
        // int nu = pois(eng);
        int nu = gsl_ran_poisson(r, mutation_rate);
        if(verbose > 1) cout << "Generating " << nu << " mutations under rate " << mutation_rate << " in cell " << cell_ID << endl;

        cna_type_probs[0] = cna_type_probs[1] = chr_prob/2;
        cna_type_probs[2] = cna_type_probs[3] = arm_prob/2;

        // cout << "Probabilities of each event type: ";
        // for(int i = 0; i < 5; i++){
        //     cout << "\t" << cna_type_probs[i];
        // }
        // cout << endl;

        gsl_ran_discrete_t* dis_chr = gsl_ran_discrete_preproc(NUM_CHR, chr_probs);
        gsl_ran_discrete_t* dis_arm = gsl_ran_discrete_preproc(2, arm_probs);
        int chr, arm, reciprocal = 1;

        for (int j=0; j < nu; j++) {
                // chr = gsl_ran_discrete(r, dis_chr);
                // double u = runiform(r, 0, 1);
                // if(u<0.5)   reciprocal = 1;

                // randomly select a type
                int e = rchoose(r, cna_type_probs);
                mut_ID += 1;

                switch (e)
                {
                case 0:
                    {
                        if(verbose > 1) cout << "\t\tchromosomal gain" << endl;
                        generate_chr_gain(mut_ID, dis_chr, reciprocal);
                        break;
                    }
                case 1:
                    {
                        if(verbose > 1) cout << "\t\tchromosomal loss" << endl;
                        generate_chr_loss(mut_ID, dis_chr, reciprocal);
                        break;
                    }
                case 2:
                    {
                        if(verbose > 1) cout << "\t\tchromosomal arm gain" << endl;
                        generate_arm_gain(mut_ID, dis_chr, dis_arm, reciprocal);
                        break;
                    }
                case 3:
                    {
                        if(verbose > 1) cout << "\t\tchromosomal arm loss" << endl;
                        generate_arm_loss(mut_ID, dis_chr, dis_arm, reciprocal);
                        break;
                    }
                case 4:
                    {
                        if(verbose > 1) cout << "\t\twhole arm region doubling" << endl;
                        generate_arm_gain(mut_ID, dis_chr, dis_arm, reciprocal);
                        generate_arm_gain(mut_ID, dis_chr, dis_arm, reciprocal);
                        break;
                    }
                default:
                    break;
                }
        }
        // cout << this->mutations.size() << endl;
        return nu;
    }


    /*
    generate a random number of CNAs
    */
    vector<int> generate_CNV(int& mut_ID, double time_occur, int verbose = 0){
        assert(arm_prob >= 0 && chr_prob >= 0);
        // A chr is not applicable to mutation if it has 0 copy. Since only relative changes are considered, assume there are sufficient copies
        gsl_ran_discrete_t* dis_chr = gsl_ran_discrete_preproc(NUM_CHR, chr_probs);
        gsl_ran_discrete_t* dis_arm = gsl_ran_discrete_preproc(2, arm_probs);
        int chr, arm, reciprocal = 1;
        double u = 0;

        int nu1 = gsl_ran_poisson(r, chr_prob);
        if(verbose > 1 && nu1 > 0) cout << "Generating " << nu1 << " chr-level mutations under rate " << chr_prob << " in cell " << cell_ID << endl;
        if(nu1 > 0){
            for (int j=0; j < nu1; j++) {
                u = runiform(r, 0, 1);
                if(u < 0.5){  // gain
                    generate_chr_gain(mut_ID, dis_chr, reciprocal, verbose);
                }else{
                    generate_chr_loss(mut_ID, dis_chr, reciprocal, verbose);
                }
            }
        }

        int nu2 = gsl_ran_poisson(r, arm_prob);
        if(verbose > 1 && nu2 > 0) cout << "Generating " << nu2 << " arm-level mutations under rate " << arm_prob << " in cell " << cell_ID << endl;
        if(nu2 > 0){
            for (int j=0; j < nu2; j++) {
                u = runiform(r, 0, 1);
                if(u < 0.5){  // gain
                    generate_arm_gain(mut_ID, dis_chr, dis_arm, reciprocal, verbose);
                }else{
                    generate_arm_loss(mut_ID, dis_chr, dis_arm, reciprocal, verbose);
                }
            }
        }

        int nu3 = 0;
        if(multi_prob>0){
            nu3 = gsl_ran_poisson(r, multi_prob);
            if(verbose > 1 && nu3 > 0) cout << "Generating " << nu3 << " multipolar divisions under rate " << multi_prob << " in cell " << cell_ID << endl;
            if(nu3 > 0){
                for (int j=0; j < nu3; j++) {
                    generate_multipolar(mut_ID, dis_chr, multi_nchr, reciprocal, verbose);
                }
            }
        }

        int nu = nu1 + nu2 + nu3;
        vector<int> nnu{nu, nu1, nu2, nu3};
        // cout << this->mutations.size() << endl;
        return nnu;
    }


    // generate fixed number of mutations
    void generate_mutations_fixed(int& mut_ID, double time_occur, int num_mut){
            for (int j=0; j < num_mut; j++) {
                    mut_ID += 1;
                    // cout << mut_ID << "\t" << time_occur << endl;
                    Mutation mut(mut_ID, time_occur);
                    this->mutations.push_back(mut);
            }
            // cout << this->mutations.size() << endl;
    }


    double get_pos_fitness(){
        double fitness = 0;
        double fmax = 0;
        if(this->ploidy>2){
            // fmax = this->ploidy / 2;
            fmax = 1;
        }
        // if(this->wgd == 1){
        //     fitness += 1;
        // }
        // uniform_real_distribution<double> runifu(0,fmax);
        // fitness = runifu(eng);

        fitness = runiform(r, 0, fmax);

        return fitness;
    }

    /*
    Increase the number of mutations undergoing WGD
    */
    void update_mut_count(int multiple){
        for (auto mut : this->mutations){
            mut.number *= multiple;
        }
    }

    int get_num_mut(){
        int sum=0;
        for (auto mut : this->mutations){
            sum += mut.number;
        }
        return sum;
    }


    /*
       This method generates CNVs, whose size follow exponential distribution
     */
    // int generate_mutations_SV(double mutation_rate, int mut_ID, double time_occur){
    //     poisson_distribution<int> pois(mutation_rate);
    //     int nu = pois(eng);
    //
    //
    //     // cout << "Generating " << nu << " mutations" << endl;
    //     for (int j=0; j < nu; j++){
    //         mut_ID += 1;
    //         // cout << mut_ID << "\t" << time_occur << endl;
    //         Mutation mut(mut_ID, time_occur);
    //         this->mutations.push_back(mut);
    //     }
    //     // cout << this->mutations.size() << endl;
    //     return nu;
    // }

    bool is_uniq_chr(int chr, vector<int>& mut_chrs){
        bool uniq = true;
        for(int j = 0; j < mut_chrs.size(); j++){
            if(mut_chrs[j] == chr){
                uniq = false;
                break;
            }
        }
        return uniq;
    }

    void generate_multipolar(int& mut_ID, gsl_ran_discrete_t* dis_chr, int multi_nchr, int reciprocal, int verbose = 0){
    // int generate_chr_gain(int& mut_ID, int chr, int reciprocal, int verbose = 0){
        // cout << "gain event on cell " << cell_ID << endl;
        vector<int> mut_chrs;
        for (int i = 0; i < multi_nchr; i++){
            int chr = gsl_ran_discrete(r, dis_chr);
            while(!(is_uniq_chr(chr, mut_chrs))){
                chr = gsl_ran_discrete(r, dis_chr);
            }
            mut_chrs.push_back(chr);
            // dis_chr[chr] = 0;

            pair<int, int> pos1(chr, 1);
            cn_profile[pos1]++;
            pair<int, int> pos2(chr, 2);
            cn_profile[pos2]++;

            Mutation mut(mut_ID, time_occur, chr, 0, 1, reciprocal);
            this->mutations.push_back(mut);

            if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1;

            // && sibling->cn_profile[pos]>0
            if(sibling!=NULL){
                sibling->cn_profile[pos1]--;
                sibling->cn_profile[pos2]--;
                Mutation mutr(mut_ID, time_occur, chr, 0, -1, reciprocal);
                sibling->mutations.push_back(mutr);
                if(verbose > 1) cout << " reciprocal event on cell " << sibling->cell_ID << endl;
            }

            mut_ID += 1;
        }

        if(verbose > 1) cout << endl;
    }

    // Assume that one chr can only has only type of event at one division
    void generate_chr_gain(int& mut_ID, gsl_ran_discrete_t* dis_chr, int reciprocal, int verbose = 0){
    // int generate_chr_gain(int& mut_ID, int chr, int reciprocal, int verbose = 0){
        // cout << "gain event on cell " << cell_ID << endl;
        int chr = gsl_ran_discrete(r, dis_chr);
        // dis_chr[chr] = 0;

        pair<int, int> pos1(chr, 1);
        cn_profile[pos1]++;
        pair<int, int> pos2(chr, 2);
        cn_profile[pos2]++;

        Mutation mut(mut_ID, time_occur, chr, 0, 1, reciprocal);
        this->mutations.push_back(mut);

        if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1;

        // && sibling->cn_profile[pos]>0
        if(reciprocal && sibling!=NULL){
            sibling->cn_profile[pos1]--;
            sibling->cn_profile[pos2]--;
            Mutation mutr(mut_ID, time_occur, chr, 0, -1, reciprocal);
            sibling->mutations.push_back(mutr);
            if(verbose > 1) cout << " reciprocal event on cell " << sibling->cell_ID << endl;
        }

        mut_ID += 1;

        if(verbose > 1) cout << endl;
    }

    void generate_chr_loss(int& mut_ID, gsl_ran_discrete_t* dis_chr, int reciprocal, int verbose = 0){
        // cout << "loss event on cell " << cell_ID << endl;
        int chr = gsl_ran_discrete(r, dis_chr);
        // dis_chr[chr] = 0;

        pair<int, int> pos1(chr, 1);
        cn_profile[pos1]--;
        pair<int, int> pos2(chr, 2);
        cn_profile[pos2]--;

        Mutation mut(mut_ID, time_occur, chr, 0, -1, reciprocal);
        this->mutations.push_back(mut);

        if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1;

        //  && sibling->cn_profile[pos]>0
        if(reciprocal && sibling!=NULL){
            sibling->cn_profile[pos1]++;
            sibling->cn_profile[pos2]++;
            Mutation mutr(mut_ID, time_occur, chr, 0, 1, reciprocal);
            sibling->mutations.push_back(mutr);
            if(verbose > 1) cout << " reciprocal event on cell " << sibling->cell_ID << endl;
        }

        mut_ID += 1;

        if(verbose > 1) cout << endl;
    }


    void generate_arm_gain(int& mut_ID, gsl_ran_discrete_t* dis_chr, gsl_ran_discrete_t* dis_arm, int reciprocal, int verbose = 0){
        // cout << "gain event on cell " << cell_ID << endl;
        int chr = gsl_ran_discrete(r, dis_chr);
        // // dis_chr[chr] = 0;
        //
        int arm = gsl_ran_discrete(r, dis_arm) + 1;

        pair<int, int> pos(chr, arm);
        // cout << "\tcopy number on cell " << cell_ID << " chr " << chr+1 << " arm " << arm << " is " << cn_profile[pos] << endl;
        // while(cn_profile[pos]==0){
        //     // cout << "\tcopy number on cell " << cell_ID << " chr " << chr+1 << " arm " << arm << " is " << cn_profile[pos] << endl;
        //     chr = gsl_ran_discrete(r, dis_chr);
        //     pos = make_pair(chr, arm);
        // }
        // // gain is only feasible with at least one copy
        // assert(cn_profile[pos] > 0);
        cn_profile[pos]++;

        Mutation mut(mut_ID, time_occur, chr, arm, 1, reciprocal);
        this->mutations.push_back(mut);

        if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1 << " arm " << arm+1;

        // && sibling->cn_profile[pos]>0
        if(reciprocal && sibling!=NULL){
            sibling->cn_profile[pos]--;
            Mutation mutr(mut_ID, time_occur, chr, arm, -1, reciprocal);
            sibling->mutations.push_back(mutr);
            if(verbose > 1) cout << " reciprocal event on cell " << sibling->cell_ID << endl;
        }

        mut_ID += 1;

        if(verbose > 1) cout << endl;
    }

    void generate_arm_loss(int& mut_ID, gsl_ran_discrete_t* dis_chr, gsl_ran_discrete_t* dis_arm, int reciprocal, int verbose = 0){
        // cout << "loss event on cell " << cell_ID << endl;
        int chr = gsl_ran_discrete(r, dis_chr);
        int arm = gsl_ran_discrete(r, dis_arm) + 1;

        pair<int, int> pos(chr, arm);
        // cout << "\tcopy number on cell " << cell_ID << " chr " << chr+1 << " arm " << arm << " is " << cn_profile[pos] << endl;
        // while(cn_profile[pos]==0){
        //     // cout << "\tcopy number on cell " << cell_ID << " chr " << chr+1 << " arm " << arm << " is " << cn_profile[pos] << endl;
        //     chr = gsl_ran_discrete(r, dis_chr);
        //     pos = make_pair(chr, arm);
        // }
        // // loss is only feasible with at least one copy
        // assert(cn_profile[pos] > 0);
        cn_profile[pos]--;

        Mutation mut(mut_ID, time_occur, chr, arm, -1, reciprocal);
        this->mutations.push_back(mut);

        if(verbose > 1) cout << "\tmutation on cell " << cell_ID << " chr " << chr+1 << " arm " << arm+1;

        //  && sibling->cn_profile[pos]>0
        if(reciprocal && sibling!=NULL){
            sibling->cn_profile[pos]++;
            Mutation mutr(mut_ID, time_occur, chr, arm, 1, reciprocal);
            sibling->mutations.push_back(mutr);
            if(verbose > 1) cout << " reciprocal event on cell " << sibling->cell_ID << endl;
        }

        mut_ID += 1;

        if(verbose > 1) cout << endl;
    }

    // only output CNAs
    void set_obs_cn(){
        obs_cn_profile.clear();
        for(int c = 0; c < NUM_CHR; c++){
            bool has_cna = false;
            if(cn_profile[pair<int, int>(c,1)] == cn_profile[pair<int, int>(c,2)] && cn_profile[pair<int, int>(c,1)] != 0){
                obs_cn_profile[pair<int, int>(c,0)] = cn_profile[pair<int, int>(c,1)];
                has_cna = true;
                continue;
            }
            if(cn_profile[pair<int, int>(c,1)]!=0){
                obs_cn_profile[pair<int, int>(c,1)] = cn_profile[pair<int, int>(c,1)];
                has_cna = true;
            }
            if(cn_profile[pair<int, int>(c,2)]!=0){
                obs_cn_profile[pair<int, int>(c,2)] = cn_profile[pair<int, int>(c,2)];
                has_cna = true;
            }
            if(!has_cna){  // no cn changes
                obs_cn_profile[pair<int, int>(c,0)] = 0;
            }
        }
    }


    // Compute the genotype differences relative to the starting cell (the sum of absolute copy number changes over all positions)
    double get_cn_diff(double chr_weight = 1.0){
        double diff = 0;
        set_obs_cn();
        for(auto cn : obs_cn_profile){
            // cout << cn.first.first << "\t" << cn.first.second << "\t" << cn.second << endl;
            int type = cn.first.second;
            if(type == 0){
                // cout << "adding weight to chr-level CNAs " << "\t" << chr_weight << endl;
                diff += chr_weight * abs(cn.second);
            }else{
                diff += abs(cn.second);
            }
        }
        // cout << "sum of CN changes: " << diff << endl;
        // diff = (double)diff/NUM_CHR;
        // cout << "average of CN changes: " << diff << endl;

        return diff;
    }


    void set_cn_all(){
        set_obs_cn();

        int k = 0;
        for(int i = 0; i < NUM_CHR; i++){
            for(int j=0; j < 3; j++){
                cn_all.push_back(0);
            }
        }
        for(auto cp: obs_cn_profile){
            int m = cp.first.first * 3 + cp.first.second;
            cn_all[m] = cp.second;
        }

        // cout << "relative copy number of cell " << cell_ID << " is ";
        // for(int i = 0; i < cn_all.size(); i++){
        //     cout << "\t" << cn_all[i];
        // }
        // cout << endl;
    }

    // output the copy numbers of cells. When all = 1, printing all the cells
    // When mutation rate is low, each pos can have at most one event
    void write_obs_cn(ofstream& fout, int all = 1){
        // fout << "\tCNAs in cell " << cell_ID << " with flag " << flag << endl;
        if(all == 0 && flag != 0) return;
        set_obs_cn();
        // convert cn_profile to final profile: assign chr-level event if cn[c][p]==cn[c][q] != 0
        for(auto cn : obs_cn_profile){
            fout << cell_ID << "\t" << cn.first.first + 1 << "\t" << cn.first.second << "\t" << cn.second << endl;
        }
    }

    // output the copy number profiles of all cells.
    void print_cn_profile(){
        // fout << "\tCNAs in cell " << cell_ID << " with flag " << flag << endl;
        for(auto cn : cn_profile){
            if(cn.second!=0){
                cout << cell_ID << "\t" << cn.first.first + 1 << "\t" << cn.first.second << "\t" << cn.second << endl;
            }
        }
    }
};



/*
to represent a population of cells
*/
class Clone
{
public:
    int clone_ID;

    vector<Cell>  cells;    // all the cells in the history, used for checking lineage history
    vector<Cell>  curr_cells;   // only available cells at present

    // variables to store subclone informaton
    vector<int> subclone_ID;  // ID of subclones
    map<int, double> subclone_time;   // Time subclone emerges
    // map<int, double> subclone_freq; // Subclone frequency
    // map<int, int> subclone_mut; // Number of mutations in subclone
    map<int, double> subclone_fitness;
    // map<int, int> subclone_division;
    // map<int, int> subclone_size;
    map<int, int> subclone_parent;
    map<int, int> subclone_psize;
    // map<int, double> mut_freq;

    int num_clonal_mutation;
    map<pair<int, int>, int> clonal_cn_profile;
    set<tuple<int, int, int>> uniq_cns;
    int num_novel_mutation;

    int tot_division;    // total number of divisions in the population
    int time_first_mut;     // which division the first mutation occurs

    // int num_cell;
    // double ploidy;
    // double time_occur;
    double time_end;    // ending time of the simulation
    double frequency;
    double fitness;
    int model;  // the model of evolution

    // Clone(int clone_ID, double time_occur, double fitness);

    ~Clone() = default;
    Clone(const Clone& other) = default;
    Clone& operator=(const Clone& other) = default;
    Clone& operator=(Clone&& other) = default;

    Clone(){
            clone_ID = 0;
            num_clonal_mutation= 0;
            num_novel_mutation= 0;
            // time_occur = 0;
            frequency = 0;
            fitness = 0; // neutral evolution
            tot_division = 0;
            time_end = 0;
    }


    Cell* get_cell_from_ID(vector<Cell>& cells, int cID){
        for(int i = 0; i < cells.size(); i++){
            Cell* cell = &cells[i];
            if(cell->cell_ID == cID) return cell;
        }
        return NULL;
    }

    void initialize(double death_rate, double mutation_rate, int& mut_ID, int& nu, int num_clonal_mutation, int verbose = 0){
            this->clone_ID = 0;
            this->cells.clear();
            this->curr_cells.clear();
            // this->mut_freq.clear();
            this->subclone_ID.clear();
            this->subclone_time.clear();
            // this->subclone_freq.clear();
         this->subclone_fitness.clear();
            this->num_clonal_mutation = num_clonal_mutation;
            this->num_novel_mutation = 0;

            Cell ncell(1, 0);

            ncell.mutation_rate = mutation_rate;

            ncell.death_rate = death_rate;
            // Genome germline(CHR_BIN_SIZE, NORM_PLOIDY);
            // ncell.genome = germline;

            cout << "generating a new cell with ID " << ncell.cell_ID << endl;
            if(num_clonal_mutation>0) {
                nu = num_clonal_mutation;
                ncell.generate_mutations_fixed(mut_ID, 0, num_clonal_mutation);
            }
            else if (mutation_rate>0) {
                nu = ncell.generate_mutations(mutation_rate, mut_ID, 0);
            }
            else{

            }

            if(verbose > 1) {
                this->cells.push_back(ncell);
            }
            this->curr_cells.push_back(ncell);
    }


    void initialize_with_cnv(int model, double birth_rate, double death_rate, double mutation_rate, double arm_prob, double chr_prob, double multi_prob, double fitness, int& mut_ID, int& nu, int num_clonal_mutation, string file_cmut="", int verbose = 0){
            this->clone_ID = 0;
            this->cells.clear();
            this->curr_cells.clear();
            // this->mut_freq.clear();
            this->subclone_ID.clear();
            this->subclone_time.clear();
            // this->subclone_freq.clear();
            this->subclone_fitness.clear();
            this->num_clonal_mutation = num_clonal_mutation;
            this->num_novel_mutation = 0;
            this->time_first_mut = 0;
            this->tot_division = 0;
            this->model = model;

            Cell ncell(1, 0);

            ncell.birth_rate = birth_rate;
            ncell.death_rate = death_rate;

            // ncell.mutation_rate = mutation_rate;
            ncell.arm_prob = arm_prob;
            ncell.chr_prob = chr_prob;
            ncell.multi_prob = multi_prob;
            ncell.mutation_rate = arm_prob + chr_prob + multi_prob;
            ncell.multi_nchr = MULTI_NCHR;
            ncell.fitness = fitness;

            // Genome germline(CHR_BIN_SIZE, NORM_PLOIDY);
            // ncell.genome = germline;

            if(verbose > 1) cout << "generating a new cell with ID " << ncell.cell_ID << endl;
            // TODO: read from file or fix mutation number
            if(file_cmut!="") {
                set_clonal_mut(file_cmut);
                for(auto cp : this->clonal_cn_profile){
                    ncell.cn_profile[cp.first] = cp.second;
                }
            }
            else if (num_clonal_mutation>0) {
                // ncell.generate_mutations_fixed(mut_ID, 0, num_clonal_mutation);
                nu = num_clonal_mutation;
                if(verbose > 1) cout << "Generating " << nu << " clonal mutations in cell " << ncell.cell_ID << endl;
                // mut_ID += nu;
                for (int j=0; j < num_clonal_mutation; j++) {
                    ncell.generate_CNV(mut_ID, 0, verbose);
                    // cout << mut_ID << "\t" << 0 << endl;
                    Mutation mut(mut_ID, 0);
                    ncell.mutations.push_back(mut);
                }
            }
            else{   // no clonal mutations
                // nu = ncell.generate_CNV(mut_ID, 0);

                // cout << "mutation probabilities of each chromosome on cell " << ncell.cell_ID << endl;
                // for(int i = 0; i < NUM_CHR; i++){
                //     ncell.chr_probs[i] = (double) 1/NUM_CHR;
                //     cout << "\t" << ncell.chr         // }
                // cout << endl;
            }

            if(verbose > 1) {
                if(num_clonal_mutation>0 || file_cmut!=""){
                    cout << "copy number profile after assigning clonal mutations" << endl;
                    ncell.print_cn_profile();
                }
            }
            this->cells.push_back(ncell);
            this->curr_cells.push_back(ncell);
    }


    /*
       lambda: The net growth rate of the background host population
       freq: The frequency of a subclone
     */
    double get_subclone_fitness(double lambda, double freq, double tend, double t1){
            double s = (lambda * t1 + log(freq / (1 - freq))) / (lambda * (tend - t1));
            return s;
    }

    /*
    This function computes the theoretical subclone frequency.
       lambda: The net growth rate of the background host population
       freq: The frequency of a subclone
     */
    double get_subclone_freq_exp(double lambda, double fitness, double tend, double t1){
            double numerator = exp(lambda * (1 + fitness) * (tend - t1));
            // double f = numerator / (numerator + exp(lambda * tend) + exp(lambda * (tend - t1)));
            double f = numerator / (numerator + exp(lambda * tend));
            return f;
    }


    void print_summary(string outfile) {
            ofstream out;
            out.setf(ios::fixed);
            out.setf(ios::showpoint);
            out.precision(9);
            out.open(outfile);

            map<int, int> subclone_nmut = get_subclone_nmut();
            map<int, double> subclone_freq = get_subclone_freq();
            map<int, int> subclone_ndiv = get_subclone_ndiv();
            map<int, double> subclone_adiv = get_subclone_adiv();
            double ploidy = get_avg_ploidy();

            out << "Average ploidy of the population: " << ploidy << endl;

            double lambda = log(2);
            out << "Information for host population:" << endl;
            for(auto cell : curr_cells){
                if(cell.clone_ID == 0){
                    lambda = cell.birth_rate - cell.death_rate;
                    out << "\tCell ID: " << cell.cell_ID << endl;
                    out << "\tMutation rate: " << cell.mutation_rate << endl;
                    out << "\tChr-level CNA rate: " << cell.chr_prob << endl;
                    out << "\tArm-level CNA rate: " << cell.arm_prob << endl;
                    out << "\tMultipolar spindle (hopeful monster) rate: " << cell.multi_prob << endl;
                    out << "\tBirth rate: " << cell.birth_rate << endl;
                    out << "\tDeath rate: " << cell.death_rate << endl;
                    out << "\tEffective mutation rate (μ/β): " << cell.mutation_rate / ((cell.birth_rate-cell.death_rate)/cell.birth_rate) << endl;
                    out << endl;
                    if(model==0) break;     // same rates under neutral evolution
                }
            }
            out << endl;
            out << "\tNumber of clonal mutation: "<< num_clonal_mutation << endl;
            out << "\tNumber of subclonal mutation: "<< num_novel_mutation << endl;
            out << "\tNumber of total divisions: "<< tot_division << endl;
            out << "\tEnd time of simulation: "<< time_end << endl;
            int num_subclone = subclone_ID.size();
            out << "Number of subclones: " << num_subclone << endl;
            if (num_subclone > 0) {
                    out << "Information for each subclone:" << endl;
                    for(int i = 0; i < num_subclone; i++){
                        out << "Subclone " << subclone_ID[i] << endl;
                        for(auto cell : curr_cells){
                            if(cell.clone_ID == subclone_ID[i]){
                                out << "\tMutation rate: " << cell.mutation_rate << endl;
                                out << "\tBirth rate: " << cell.birth_rate << endl;
                                out << "\tDeath rate: " << cell.death_rate << endl;
                                out << "\tEffective mutation rate (μ/β): " << cell.mutation_rate / ((cell.birth_rate-cell.death_rate)/cell.birth_rate) << endl;
                                break;
                            }
                        }
                        out << "\tFrequency: " << subclone_freq[subclone_ID[i]] << endl;
                        out << "\tNumber of mutations in subclone: " << subclone_nmut[subclone_ID[i]] << endl;
                        out << "\tFitness advantage: " << subclone_fitness[subclone_ID[i]] << endl;
                        out << "\tTime subclone emerges (simulation time): " << subclone_time[subclone_ID[i]] << endl;
                        out << "\tNumber of divisions: " << subclone_ndiv[subclone_ID[i]] << endl;
                        out << "\tAverage number of divisions per cell: " << subclone_adiv[subclone_ID[i]] << endl;
                        out << "\tPopulation size when subclone emerges: " << subclone_psize[subclone_ID[i]] << endl;
                        out << "\tTime subclone emerges (tumor doublings): " << log(subclone_psize[subclone_ID[i]])/(lambda) << endl;
                        // out << "\tParent of subclone (0 is host): " << subclone_parent[[subclone_ID[i]]];
                        out << endl;
                    }
            }
            else{
                    out << "No clones\n\n";
            }
            out.close();
    }


    double get_rmax(){
            // Find  maximum birth rate (bmax) and maximum death rate (dmax) in the population
            double bmax = 0;
            double dmax = 0;
            for(unsigned int i=0; i<curr_cells.size(); i++) {
                    Cell ci = curr_cells[i];
                    double bi = ci.birth_rate;
                    double di = ci.death_rate;
                    if (bi > bmax) {
                            bmax = bi;
                    }
                    if(di > dmax) {
                            dmax = di;
                    }
            }
            // cout << "Maximum birth rate: " << bmax << endl;
            // cout << "Maximum death rate: " << dmax << endl
            return bmax + dmax;
    }

    /*
       This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
       intput:
        Nend -- the number of cells in the final population;
        mutation_rate -- the Mutation rate per Cell division;
        time_occur -- defined in terms of population doublings
       output:
        a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
     */
    void grow(int num_subclone, int num_clonal_mutation, vector<double> fitness, vector<double> time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, int verbose = 0){
            // time is defined in terms of population doublings
            double t = 0;
            int cell_count = 0;     // To count the total number of cells in histor
            // Initialize the simulation with one Cell
            int mut_ID = 0;
            int nu = 0; // The number of new mutations

            unsigned int fitmutant = 0;
            vector<double> timeN_occur;
            double lambda = birth_rate - death_rate;
            // convert time_occur from tumor doublings time to real time
            cout << "Subclone occurring time (tumor cell number):" << endl;
            for (unsigned int i = 0; i< time_occur.size(); i++){
                timeN_occur.push_back(ceil(exp(lambda * time_occur[i])));
                cout << "\t"  << timeN_occur[i];
            }
            // timeN_occur.push_back(0);
            cout << "\n";

            initialize(death_rate, mutation_rate, mut_ID, nu, num_clonal_mutation, verbose);
            cell_count += 1;
            while(this->curr_cells.size() < Nend) {
                    if (this->curr_cells.size() == 0) {
                        t = 0;
                        mut_ID = 0;
                        nu = 0;
                        cell_count = 0;
                        initialize(death_rate, mutation_rate, mut_ID, nu, num_clonal_mutation, verbose);
                        cell_count += 1;
                        continue;
                    }
                    // Choose a random Cell from the current population
                    // uniform_int_distribution<> iunif(0, this->curr_cells.size() - 1);
                    // int rindex = iunif(eng);
                    int rindex = myrng(this->curr_cells.size());

                    Cell rcell = this->curr_cells[rindex];
                    // cout << "Selecting " << rindex+1 << "th cell" << endl;
                    // cout << "Selecting cell " << rcell.cell_ID << endl;
                    double rmax = get_rmax();

                    // increase time
                    // uniform_real_distribution<double> runifu(0,1);
                    // double tau = -log(runifu(eng)); // an exponentially distributed random variable
                    double tau = -log(runiform(r, 0, 1));
                    double deltaT = tau/(rmax * this->curr_cells.size());
                    t += deltaT;

                    // draw a random number
                    // uniform_real_distribution<double> runif(0.0,rmax);
                    // double rb = runif(eng);
                    double rb = runiform(r, 0, rmax);
                    // cout << "random number " << rb << endl;
                    // birth event if r<birthrate
                    if(rb < rcell.birth_rate) {
                            // cout << "Number of generated cells " << this->cells.size() << endl;
                            // increase one Cell
                            // cout << "birth event at time " << t << endl;
                            int parent_ID = rcell.cell_ID;
                            // cout << "   parent " << parent_ID << endl;
                            Cell dcell1 = Cell(rcell);
                            dcell1.cell_ID = cell_count + 1;
                            dcell1.parent_ID = parent_ID;
                            dcell1.num_division = rcell.num_division + 1;

                            Cell dcell2 = Cell(rcell);
                            dcell2.cell_ID = cell_count + 2;
                            dcell2.parent_ID = parent_ID;
                            dcell2.num_division = rcell.num_division + 1;
                            cell_count += 2;
                            // cout << "   children " << dcell1.cell_ID  << "\t" << dcell2.cell_ID << endl;
                            // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                            if (mutation_rate>0) {
                                    nu +=  dcell1.generate_mutations(mutation_rate, mut_ID, t);
                                    // mut_ID += nu;
                                    nu +=  dcell2.generate_mutations(mutation_rate, mut_ID, t);
                                    // mut_ID += nu;
                            }
                            // introduce a fitter mutatant
                            if(fitmutant < num_subclone && this->curr_cells.size() >= timeN_occur[fitmutant]) {
                                    double fitval = fitness[fitmutant];
                                    cout << "Introducing a fitter mutatant with fitness " << fitval << endl;
                                    dcell1.fitness = fitval;
                                    // dcell1.death_rate = runifu(eng) * rcell.death_rate;
                                    dcell1.death_rate = runiform(r, 0, 1) * rcell.death_rate;
                                    dcell1.birth_rate = (1 + dcell1.fitness) * (rcell.birth_rate - rcell.death_rate) + dcell1.death_rate;
                                    int clone_ID = this->clone_ID + fitmutant + 1;
                                    cout << "   birth rate: " << dcell1.birth_rate << endl;
                                    cout << "   death rate: " << dcell1.death_rate << endl;
                                    dcell1.clone_ID = clone_ID;
                                    this->subclone_ID.push_back(clone_ID);
                                    this->subclone_fitness[clone_ID] = fitval;
                                    this->subclone_time[clone_ID] = t;
                                    this->subclone_psize[clone_ID] = this->curr_cells.size();
                                    fitmutant += 1;
                            }
                            if(verbose > 1) {
                                    for(unsigned int i = 0; i < this->cells.size(); i++) {
                                            if (this->cells[i].cell_ID==rcell.cell_ID) {
                                                    rcell.num_division += 1;
                                                    rcell.flag = 1;
                                                    this->cells[i] = rcell;
                                                    break;
                                               }
                                    }
                            }
                            // Remove the parent cell from the list of current cells
                            // cout << "Removing cell " << this->curr_cells[rindex].cell_ID << endl;
                            this->curr_cells.erase(this->curr_cells.begin()+rindex);
                            dcell1.time_occur = t;
                            dcell2.time_occur = t;
                            this->curr_cells.push_back(dcell1);
                            this->curr_cells.push_back(dcell2);

                            if(verbose > 1) {
                                this->cells.push_back(dcell1);
                                this->cells.push_back(dcell2);
                            }
                    }
                    // death event if b<r<b+d
                    if(rb >= rcell.birth_rate && rb < rcell.birth_rate + rcell.death_rate) {
                            // cout << " death event" << endl;
                            if(verbose > 1) {
                                    for(unsigned int i = 0; i < this->cells.size(); i++) {
                                            if (this->cells[i].cell_ID==rcell.cell_ID) {
                                                    rcell.flag = -1;
                                                    this->cells[i] = rcell;
                                                    break;
                                            }
                                    }
                            }
                            this->curr_cells.erase(this->curr_cells.begin()+rindex);
                    }
                    // cout << "===========================" << endl;
            }
            // this->set_ploidy(this->curr_cells);
            // cout << "The average ploidy of this clone: " << this->ploidy << endl;
            cout << "Generated " << cell_count << " cells with " << nu << " mutations"  << " in " << tot_division << " divisions"<< endl;
            cout << "End time: " << t << endl;
    }


    /*
       This method simulates tumour growth with a rejection-kinetic Monte Carlo algorithm.
       intput:
        Nend -- the cell population size at the end
        mutation_rate -- the mutation rate per cell division;
        model -- 0: neutral, 1: gradual, 2: punctuated
       output:
        a tree-like structure. For each Cell, its children, occurence time, birth rate, death rate
     */
     void grow_with_cnv(int num_subclone, int num_clonal_mutation, int model, const vector<double>& fitness, const vector<double>& time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, double arm_prob, double chr_prob, double multi_prob, int genotype_diff = 0, double chr_weight = 1.0, int chr_sel = 0, string file_cmut = "", int verbose = 0){
             double t = 0;
             int cell_count = 0;     // To count the total number of cells in history
             int mut_ID = 0;
             int nu = 0; // The number of new mutations

             // Initialize the simulation with one cell
             initialize_with_cnv(model, birth_rate, death_rate, mutation_rate, arm_prob, chr_prob, multi_prob, 0, mut_ID, nu, num_clonal_mutation, file_cmut, verbose);
             cell_count += 1;

             if(verbose > 1) cout << "simulating tumour growth with CNA" << endl;
             unsigned int fitmutant = 0;
             int num_mut_event = 0; // count the number of times a CNA event is introduced

             while(this->curr_cells.size() < Nend) {
                 if (this->curr_cells.size() == 0) {
                     t = 0;
                     mut_ID = 0;
                     nu = 0;
                     cell_count = 0;
                     initialize_with_cnv(model, birth_rate, death_rate, mutation_rate, arm_prob, chr_prob, multi_prob, 0, mut_ID, nu, num_clonal_mutation, file_cmut, verbose);
                     cell_count += 1;
                     continue;
                 }
                 // Choose a random Cell from the current population
                 int rindex = myrng(this->curr_cells.size());

                 Cell rcell = this->curr_cells[rindex];
                 // cout << "Selecting " << rindex+1 << "th cell" << endl;
                 // cout << "Selecting cell " << rcell.cell_ID << endl;
                 double rmax = get_rmax();
                 if(verbose > 1) cout << "max sum of birth rate and death rate is " << rmax << endl;
                 // increase time
                 double tau = -log(runiform(r, 0, 1));
                 double deltaT = tau/(rmax * this->curr_cells.size());
                 t += deltaT;

                 // draw a random number
                 double rb = runiform(r, 0, rmax);
                 // cout << "random number " << rb << endl;
                 // birth event if r<birthrate
                 if(rb < rcell.birth_rate) {
                     // cout << "Number of generated cells " << this->cells.size() << endl;
                     // increase one cell
                     if(verbose > 1){
                         cout << "birth rate for cell "<< rcell.cell_ID  << " is " << rcell.birth_rate << endl;
                     }
                     // cout << "birth event at time " << t << endl;
                     int parent_ID = rcell.cell_ID;
                     // cout << "   parent " << parent_ID << endl;
                     Cell dcell1 = Cell(rcell);
                     // cout << "birth rates: " << dcell1.birth_rate << "\t" << rcell.birth_rate << endl;
                     dcell1.cell_ID = cell_count + 1;
                     dcell1.parent_ID = parent_ID;
                     dcell1.num_division = rcell.num_division + 1;
                     copy(rcell.chr_probs, rcell.chr_probs + NUM_CHR, dcell1.chr_probs);

                     Cell dcell2 = Cell(rcell);
                     dcell2.cell_ID = cell_count + 2;
                     dcell2.parent_ID = parent_ID;
                     dcell2.num_division = rcell.num_division + 1;
                     copy(rcell.chr_probs, rcell.chr_probs + NUM_CHR, dcell2.chr_probs);

                     dcell1.sibling = &dcell2;
                     dcell2.sibling = &dcell1;

                     rcell.daughters.push_back(dcell1.cell_ID);
                     rcell.daughters.push_back(dcell2.cell_ID);

                     cell_count += 2;
                     tot_division += 1;
                     // cout << "   children " << dcell1.cell_ID  << "\t" << dcell2.cell_ID << endl;
                     // daughter cells aquire nu new mutations, where nu ~ Poisson(mutation_rate)
                     if (mutation_rate > 0) {
                         // nu1: a vector of 4 numbers (#total CNAs, #chr-level CNAs, #arm-level CNAs, #MP CNAs)
                         vector<int> nu1 = dcell1.generate_CNV(mut_ID, t, verbose);
                         nu += nu1[0];
                         vector<int> nu2 = dcell2.generate_CNV(mut_ID, t, verbose);
                         nu += nu2[0];

                         if(verbose > 1){
                             cout << "Number of mutations (total, chr-level, arm-level) in cell " << dcell1.cell_ID << ": ";
                             for(int i=0; i<nu1.size(); i++){
                                 cout << "\t" << nu1[i];
                             }
                             cout << endl;

                             cout << "Number of mutations (total, chr-level, arm-level) in cell " << dcell2.cell_ID << ": ";
                             for(int i=0; i<nu2.size(); i++){
                                 cout << "\t" << nu2[i];
                             }
                             cout << endl;
                        }

                         int sentinel_nu1 = nu1[0];
                         int sentinel_nu2 = nu2[0];
                         if(chr_sel == 1){  // only assume selection on chr-level events
                             sentinel_nu1 = nu1[1];
                             sentinel_nu2 = nu2[1];
                         }
                         if(sentinel_nu1 > 0 || sentinel_nu2 > 0) // all simulated CNAs are reciprocal
                         {
                            num_mut_event += 1;     // count the order of mutations introduced
                            if(num_mut_event == 1)  time_first_mut = tot_division;
                            // Introduce selection to cells with CNAs using specified fitness
                            if(model > 0) assert(fitness.size() > 0);

                            double gdiff1 = 1;
                            double gdiff2 = 1;
                            if(genotype_diff > 0){
                                gdiff1 = dcell1.get_cn_diff(chr_weight);
                                gdiff2 = dcell2.get_cn_diff(chr_weight);
                            }

                            if(model == 2 && num_mut_event == 1){ // only introduce selection at first hit
                                 // dcell1.birth_rate = (1 + gdiff1 * fitness[0]) * birth_rate;
                                 // dcell2.birth_rate = (1 + gdiff2 * fitness[0]) * birth_rate;
                                 if(genotype_diff > 0){
                                     dcell1.birth_rate = birth_rate / (1 + gdiff1 * fitness[0]);
                                     dcell2.birth_rate = birth_rate / (1 + gdiff2 * fitness[0]);
                                 }else{
                                     dcell1.birth_rate = rcell.birth_rate * (1 + fitness[0]);
                                     dcell2.birth_rate = rcell.birth_rate * (1 + fitness[0]);
                                 }

                                 if(verbose > 1){
                                     cout << "new birth rate for cell "<< dcell1.cell_ID  << " is " << dcell1.birth_rate << " with genotype difference " << gdiff1 << endl;
                                     cout << "new birth rate for cell "<< dcell2.cell_ID  << " is " << dcell2.birth_rate << " with genotype difference " << gdiff2 << endl;
                                 }
                            }
                            if(model == 1){     // gradual selection, each CNA leads to fitness change of daughter cells
                                 if(genotype_diff > 0){ // The genotype changes are based on initial cell
                                     dcell1.birth_rate = birth_rate / (1 + gdiff1 * fitness[0]);
                                     dcell2.birth_rate = birth_rate / (1 + gdiff2 * fitness[0]);
                                 }else{
                                     dcell1.birth_rate = rcell.birth_rate * (1 + fitness[0]);
                                     dcell2.birth_rate = rcell.birth_rate * (1 + fitness[0]);
                                 }
                                 if(verbose > 1){
                                     cout << "new birth rate for cell "<< dcell1.cell_ID  << " is " << dcell1.birth_rate << " with genotype difference " << gdiff1 << endl;
                                     cout << "new birth rate for cell "<< dcell2.cell_ID  << " is " << dcell2.birth_rate << " with genotype difference " << gdiff2 << endl;
                                 }
                            }
                            // cout << "fitness " << fitness[0] << endl;
                            if(model == 3){
                                dcell1.birth_rate = birth_rate * (1 + gdiff1 * fitness[0]);
                                dcell2.birth_rate = birth_rate * (1 + gdiff2 * fitness[0]);
                                if(verbose > 1){
                                    cout << "new birth rate for cell "<< dcell1.cell_ID  << " is " << dcell1.birth_rate << " with genotype difference " << gdiff1 << endl;
                                    cout << "new birth rate for cell "<< dcell2.cell_ID  << " is " << dcell2.birth_rate << " with genotype difference " << gdiff2 << endl;
                                }
                            }
                         }
                     }

                     // if(verbose > 1) {
                         for(unsigned int i = 0; i < this->cells.size(); i++) {
                             if (this->cells[i].cell_ID==rcell.cell_ID) {
                                     rcell.num_division += 1;
                                     rcell.flag = 1;
                                     this->cells[i] = rcell;
                                     break;
                             }
                         }
                     // }
                     // Remove the parent cell from the list of current cells
                     // cout << "Removing cell " << this->curr_cells[rindex].cell_ID << endl;
                     this->curr_cells.erase(this->curr_cells.begin()+rindex);
                     dcell1.time_occur = t;
                     dcell2.time_occur = t;
                     this->curr_cells.push_back(dcell1);
                     this->curr_cells.push_back(dcell2);

                     // if(verbose > 1) {
                         this->cells.push_back(dcell1);
                         this->cells.push_back(dcell2);
                     // }
                 }
                 // death event if b<r<b+d
                 if(rb >= rcell.birth_rate && rb < rcell.birth_rate + rcell.death_rate) {
                     // cout << " death event" << endl;
                     // if(verbose > 1) {
                         for(unsigned int i = 0; i < this->cells.size(); i++) {
                             if (this->cells[i].cell_ID==rcell.cell_ID) {
                                 rcell.flag = -1;
                                 this->cells[i] = rcell;
                                 break;
                             }
                         }
                     // }
                     this->curr_cells.erase(this->curr_cells.begin()+rindex);
                 }
                 // cout << "===========================" << endl;
             }
             // this->set_ploidy(this->curr_cells);
             // cout << "The average ploidy of this clone: " << this->ploidy << endl;
             time_end = t;
             num_novel_mutation = nu;
             if(verbose > 1) cout << "Generated " << cell_count << " cells in total with " << nu << " mutations in " << tot_division << " divisions during time " << t << endl;

     }

     /*
        This method prints out the sum of ratios of branch lenghts before and after a division (skipping first k branches)
      */
      void get_treelen_ratios(vector<double>& avg, vector<Cell>& cells, int skip = 3){
         double sum1 = 0, sum2 = 0;
         int n1 = 0, n2 = 0;
         // svector<float> uniq_ratios;
         // cout << "\nbranch length ratios:";
         for(unsigned int i = 0; i < cells.size() ; i++) {
             Cell *cell = &cells[i];
             if(cell->parent_ID <= skip) continue;

             Cell* pcell = cell->get_parent(cells);
             if(pcell->parent_ID <= skip) continue;

             if(cell->daughters.size() <= 0) continue;

             Cell* ppcell = pcell->get_parent(cells);
             double branch_len1 = pcell->time_occur - ppcell->time_occur;

             Cell* d = get_cell_from_ID(cells, cell->daughters[0]);
             double branch_len2 = d->time_occur - cell->time_occur;
             double ratio = branch_len2 / branch_len1;
             // uniq_ratios.push_back(ratio);
             if(ratio > 1){
                sum1 += ratio;
                n1 += 1;
            }else{
                sum2 += ratio;
                n2 += 1;
            }
             // cout << "\t" << ratio;
         }
         // cout << endl;
         // cout << sum1 << "\t" << n1 << endl;
         // cout << sum2 << "\t" << n2 << endl;
         double avg1 = sum1 / n1;
         double avg2 = sum2 / n2;
         avg.push_back(avg1);
         avg.push_back(avg2);

         // return avg;
     }


     /*
        This method prints out unique branch lenghts (skipping first k branches)
      */
     void get_treelen_vec(vector<double>& blens, vector<Cell>& cells, int skip = 3){
         double sum = 0;
         for(unsigned int i = 0; i < cells.size() ; i++) {
             Cell cell = cells[i];
             if(cell.parent_ID <= skip) continue;
             // write the cell lineages in a file for visualization
             Cell* pcell = cell.get_parent(cells);
             double branch_len = cell.time_occur - pcell->time_occur;
             blens.push_back(branch_len);
         }
     }


     /*
        This method prints out half the sum of branch lenghts (skipping first k branches)
      */
     double get_treelen(vector<Cell>& cells, int skip = 3){
         double sum = 0;
         for(unsigned int i = 0; i < cells.size() ; i++) {
                 Cell cell = cells[i];
                 if(cell.parent_ID <= skip) continue;
                 // if (i <= skip){
                 //     continue;
                 // }
                 // write the cell lineages in a file for visualization
                 Cell* pcell = cell.get_parent(cells);
                 double branch_len = cell.time_occur - pcell->time_occur;
                 // if (i <= skip){
                 //     cout << "Skip branch " << i << "\t" << cell.parent_ID << "\t" << cell.cell_ID << "\t" << branch_len << endl;
                 //     continue;
                 // }
                 sum += branch_len;

         }
         sum = sum / 2;
         return sum;
     }


     /*
        This method prints out the lineage of cells in a clone in the format of edge list
      */
     void write_tree(vector<Cell>& cells, string fname){
         ofstream tout(fname);
         tout.precision(9);

         // string header = "parent\tchild\tbranch_len\n";
         // tout << header;

         for(unsigned int i = 0; i < cells.size() ; i++) {
                 Cell cell = cells[i];
                 if(cell.parent_ID == 0) continue;

                 // write the cell lineages in a file for visualization
                 Cell* pcell = cell.get_parent(cells);
                 double branch_len = cell.time_occur - pcell->time_occur;
                 tout << cell.parent_ID << "\t" << cell.cell_ID << "\t" << branch_len << endl;
         }
         tout.close();
     }


     // Find the preorder of nodes in the genealogy tree
     void get_nodes_preorder(vector<Cell>& cells, Cell* root, vector<Cell*>& nodes_preorder){
         nodes_preorder.push_back(root);
         for(int j=0; j < root->daughters.size();j++){
             Cell* d = get_cell_from_ID(cells, root->daughters[j]);
             get_nodes_preorder(cells, d, nodes_preorder);
         }
     }

     /*
        This method prints out the lineage of cells in a clone in the format of NEWICK
      */
     void write_newick(vector<Cell>& cells, string fname){
         ofstream tout(fname);
         int precision = 9;
         tout.precision(precision);
         string newick = "";
         const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
         const boost::format internal_node_format(boost::str(boost::format(")%%d:%%.%df") % precision));
         const boost::format root_node_format(boost::str(boost::format(")%%d")));
         stack<Cell*> node_stack;
         vector<Cell*> nodes_preorder;
         Cell* root;

         for(int i=0; i<cells.size(); ++i){
           if(cells[i].parent_ID == 0){
               root = &cells[i];
               break;
           }
         }

         // cout << "Get nodes preorder" << endl;
         get_nodes_preorder(cells, root, nodes_preorder);
         // for(int i = 0; i < nodes_preorder.size(); i++){
         //     cout << "\t" << nodes_preorder[i]->cell_ID;
         // }
         // cout << endl;

         Cell* pcell;
         double branch_len = 0;

         // cout << "Traverse nodes in preorder" << endl;
         // Traverse nodes in preorder
         for (int i = 0; i<nodes_preorder.size(); i++)
         {
             Cell* nd = nodes_preorder[i];
             // cout << nd->cell_ID << endl;
             if (nd->daughters.size()>0) // internal nodes
             {
                 newick += "(";
                 node_stack.push(nd);
             }
             else
             {
                 pcell = nd->get_parent(cells);
                 // cout << nd->cell_ID << "\t" << pcell->cell_ID << "\t" << pcell->daughters[0] << "\t" << pcell->daughters[1] << endl;
                 // assert(pcell->daughters.size()>0);
                 branch_len = nd->time_occur - pcell->time_occur;

                 newick += boost::str(boost::format(tip_node_format) % (nd->cell_ID) % branch_len);

                 if (nd->cell_ID == pcell->daughters[0])   //left child
                 {
                     newick += ",";
                 }
                 else
                 {
                     Cell* popped = (node_stack.empty() ? 0 : node_stack.top());
                     pcell = popped->get_parent(cells);
                     while (popped && popped->parent_ID > 0 && popped->cell_ID == pcell->daughters[1]) // right sibling of the previous node
                     {
                         // cout << popped->cell_ID << "\t" << pcell->cell_ID << "\t" << pcell->daughters[0] << "\t" << pcell->daughters[1] << endl;
                         branch_len = popped->time_occur - pcell->time_occur;

                         node_stack.pop();

                         newick += boost::str(boost::format(internal_node_format) % (popped->cell_ID) % branch_len);
                         popped = node_stack.top();
                         pcell = popped->get_parent(cells);
                     }

                     //  cout << popped->cell_ID << "\t" << pcell->cell_ID << "\t" << pcell->daughters[0] << "\t" << pcell->daughters[1] << endl;
                     if (popped && popped->parent_ID > 0 && popped->cell_ID == pcell->daughters[0]) // left child, with another sibling
                     {
                         branch_len = popped->time_occur - pcell->time_occur;

                         node_stack.pop();

                         newick += boost::str(boost::format(internal_node_format) % (popped->cell_ID) % branch_len);
                         newick += ",";
                     }

                     if (node_stack.empty())
                     {
                         newick += ")";
                     }
                 }
             }
             // cout << newick << endl;
         }
         newick +=  boost::str(boost::format(root_node_format) % (root->cell_ID));
         newick += ";";

         tout << newick;
         tout.close();
     }


    /*
       This method prints out the details of cells in a clone
     */
    void print_lineage(vector<Cell>& cells, string outdir, string suffix, int verbose = 0){
            string fname = outdir + "all_cells_lineage" + suffix + ".txt";
            ofstream out(fname);
            out.precision(9);

            fname= outdir + "all_cells_cn" + suffix + ".txt";
            ofstream fcn(fname);

            fname= outdir + "all_cells_mut" + suffix + ".txt";
            ofstream fmut(fname);

            int num_cell = cells.size();
            if(verbose > 1) cout << "There are " << num_cell << " cells in the history of tumor growth" << endl;

            // string header = "id\tparent_ID\tflag\tbirth_rate\tdeath_rate\tmutation_rate\tploidy\tnum_division\ttime_occur\n";
            string header = "id\tparent_ID\tflag\tnum_mut\tclone_ID\ttime_occur\tbranch_len\n";
            out << header;

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = cells[i];
                    if(cell.parent_ID == 0) continue;
                    // if(cell.daughters.size()==0){
                    //     cout << cell.cell_ID << "\t" << cell.parent_ID << "\t" << -1 << "\t" << -1 << "\t" << cell.flag << "\t" << cell.birth_rate << "\t" << cell.death_rate << "\t" << cell.mutation_rate  << "\t" << cell.ploidy  << "\t" << cell.num_division << "\t" << cell.time_occur << endl;
                    // }else{
                    //     cout << cell.cell_ID << "\t" << cell.parent_ID << "\t" << cell.daughters[0] << "\t" << cell.daughters[1] << "\t" << cell.flag << "\t" << cell.birth_rate << "\t" << cell.death_rate << "\t" << cell.mutation_rate  << "\t" << cell.ploidy  << "\t" << cell.num_division << "\t" << cell.time_occur << endl;
                    // }
                    int num_mut = cell.get_num_mut();
                    // int num_mut = cell.mutations.size();

                    // write the cell lineages in a file for visualization
                    Cell* pcell = cell.get_parent(cells);
                    double branch_len = cell.time_occur - pcell->time_occur;
                    out << cell.cell_ID << "\t" << cell.parent_ID << "\t" << cell.flag << "\t" << num_mut << "\t" << cell.clone_ID  << "\t" << cell.time_occur << "\t" << branch_len << endl;

                    cell.write_obs_cn(fcn, 1);

                    if(verbose > 1) cout << "\t" << num_mut << " mutations in cell " << cell.cell_ID << endl;
                    for(unsigned int j = 0; j < num_mut; j++){
                        Mutation mut = cell.mutations[j];
                        fmut << cell.cell_ID << "\t" << mut.mut_ID << "\t" << mut.time_occur << "\t" << mut.chr + 1 << "\t" << mut.arm << "\t" << mut.type << "\t" << mut.reciprocal << endl;
                    }
            }
            out.close();
            fcn.close();
            fmut.close();
    }

    /*
       This method computes the allele frequency of each mutation
     */
    map<int, double> get_allele_freq(){
            map<int, double> mut_freq;

            int num_cell = this->curr_cells.size();
            // Collect mutations
            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    for(auto mut : cell.mutations) {
                            mut_freq[mut.mut_ID] += mut.number;
                    }
            }
            // Compute frequency
            double ploidy = get_avg_ploidy();
            // cout << "Average ploidy of the population: " << ploidy << endl;
            for(auto it : mut_freq) {
                    double vaf = it.second/num_cell;
                    vaf = vaf / ploidy;
                    mut_freq[it.first] = vaf;
            }
            return mut_freq;
    }


    /*
       This method computes the average ploidy of a tumor population.
     */
    double get_avg_ploidy(){
            double ploidy = 0;
            int num_cell = this->curr_cells.size();
            // Collect mutations
            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    ploidy += cell.ploidy;
            }
            ploidy = ploidy / num_cell;
            return ploidy;
    }

    /*
       This method computes the number of unique mutations in each subclone.
     */
    map<int, int> get_subclone_nmut(){
            map<int, set<double>> subclone_muts;
            map<int, int> subclone_nmut;

            int num_cell = this->curr_cells.size();

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    for (auto mut : cell.mutations){
                        subclone_muts[cell.clone_ID].insert(mut.mut_ID);
                    }
            }

            for(auto muts : subclone_muts) {
                    subclone_nmut[muts.first] = muts.second.size();
            }

            return subclone_nmut;
    }

    /*
       This method computes the maximum number of cell division in each subclone.
     */
    map<int, int> get_subclone_ndiv(){
            map<int, set<int>> subclone_divs;
            map<int, int> subclone_ndiv;

            int num_cell = this->curr_cells.size();

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    subclone_divs[cell.clone_ID].insert(cell.num_division);
            }

            // cout << "Number of unique divisions in subclones: " << endl;
            for(auto divs : subclone_divs) {
                // for(auto num : divs.second){
                //     cout << num << "\t";
                // }
                // cout << endl;
                // set<int>::iterator min = divs.second.begin();
                set<int>::reverse_iterator max = divs.second.rbegin();
                subclone_ndiv[divs.first] = *max;
                // cout << divs.first << "\t" << *min << "\t" << *max << "\n";
            }
            return subclone_ndiv;
    }


    /*
       This method computes the average number of cell division in each subclone.
     */
    map<int, double> get_subclone_adiv(){
            map<int, vector<int>> subclone_divs;
            map<int, double> subclone_adiv;

            int num_cell = this->curr_cells.size();

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    // cout << cell.cell_ID << "\t" << cell.clone_ID << "\t" << cell.num_division << endl;
                    subclone_divs[cell.clone_ID].push_back(cell.num_division);
            }

            // cout << "Number of divisions in subclones: " << endl;
            for(auto divs : subclone_divs) {
                int sum = 0;
                for (auto num : divs.second)
                    sum += num;
                subclone_adiv[divs.first] = sum / (divs.second).size();
                // cout << "\t" << divs.first << "\t"<< sum << "\t" << (divs.second).size() << "\n";
            }

            return subclone_adiv;
    }


    /*
       This method computes the subclone frequencies of a tumor population.
     */
    map<int, double> get_subclone_freq(){
            map<int, double> subclone_freq;
            int num_cell = this->curr_cells.size();

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = this->curr_cells[i];
                    subclone_freq[cell.clone_ID] += 1;
            }

            for(auto freq : subclone_freq) {
                    subclone_freq[freq.first] = freq.second / num_cell;
            }

            return subclone_freq;
    }

    // void initialize(double death_rate, double mutation_rate, int& mut_ID, int& nu, int num_clonal_mutation, int verbose);
    // void initialize_with_cnv(double birth_rate, double death_rate, double mutation_rate, double arm_prob, double chr_prob, int& mut_ID, int& nu, int num_clonal_mutation, string file_cmut, int verbose);
    //
    // double get_subclone_fitness(double lambda, double freq, double tend, double t1);
    // double get_subclone_freq_exp(double lambda, double fitness, double tend, double t1);
    // double get_rmax();
    // double get_avg_ploidy();
    // map<int, double> get_allele_freq();
    // map<int, double> get_subclone_freq();
    // map<int, int> get_subclone_nmut();
    // map<int, int> get_subclone_ndiv();
    // map<int, double> get_subclone_adiv();
    //
    // void grow(int num_subclone, int num_clonal_mutation, vector<double> fitness, vector<double> time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, int verbose);
    // void grow_with_cnv(int num_subclone, int num_clonal_mutation, vector<double> fitness, vector<double> time_occur, int Nend, double birth_rate, double death_rate, double mutation_rate, double arm_prob, double chr_prob, string file_cmut, int verbose);
    //
    //
    // void print_lineage(vector<Cell> cells, string outdir, string suffix, int verbose);
    // void print_summary(string outfile);

    template <typename T> void print_map(map<int, T> m, string outfile){
        ofstream out;
        out.setf(ios::fixed);
        out.setf(ios::showpoint);
        out.precision(9);
        out.open(outfile);

        for(auto it : m) {
                out << it.first << "\t" << it.second << endl;
        }

        out.close();
    }

    /*
       This method prints out the copy numbers of final cells in a clone
     */
    void print_obs_cn(vector<Cell> cells, string fname, int verbose = 0){
            ofstream fcn(fname);

            int num_cell = cells.size();
            if(verbose > 1) cout << "Printing copy numbers of " << num_cell << " cells" << endl;

            for(unsigned int i = 0; i < num_cell; i++) {
                    Cell cell = cells[i];
                    cell.write_obs_cn(fcn, 0);
            }
            fcn.close();
    }

    /*
       This method read clonal mutations from file
     */
    void set_clonal_mut(string fname){
        ifstream infile(fname.c_str());
        if (infile.is_open()){
          std::string line;
          while(!getline(infile,line).eof()){
            if(line.empty()) continue;

            std::vector<std::string> split;
            std::string buf;
            stringstream ss(line);
            while (ss >> buf) split.push_back(buf);
            assert(split.size()==3);

            int chr = atoi(split[0].c_str()) - 1;
            int arm = atoi(split[1].c_str());
            int cn = atoi(split[2].c_str());
            // cout << "read chr: " << chr << " arm: " << arm << " cn: " << cn << endl;
            pair<int, int> pos(chr, arm);
            clonal_cn_profile[pos] = cn;
          }
      }
    }

    /*
       This method summarize subclonal mutations in terms of frequencies of events
     */
     void write_prop(vector<Cell> cells, string fname){
         // compute the frequencies of arm-level events by chromosome
         map<int, double> parm_freq;
         map<int, double> qarm_freq;
         for(int i = 0; i < NUM_CHR; i++){
             parm_freq[i] = 0.0;
             qarm_freq[i] = 0.0;
         }
         int num_parm_cnv = 0;
         int num_qarm_cnv = 0;

         for(unsigned int i = 0; i < cells.size(); i++) {
             Cell cell = cells[i];
             cell.set_obs_cn();
             // cout << "cell " << cell.cell_ID << endl;
             for(auto cp : cell.obs_cn_profile){
                 // cout << cp.first.first << "\t" << cp.first.second << "\t" << cp.second << endl;
                 if(cp.first.second == 1 && cp.second != 0){
                     num_parm_cnv += 1;
                     parm_freq[cp.first.first] += 1;
                 }
                 if(cp.first.second == 2 && cp.second != 0){
                     num_qarm_cnv += 1;
                     qarm_freq[cp.first.first] += 1;
                 }
             }
         }
         // cout << "There are " << num_parm_cnv << " parm-level events" << endl;
         if(num_parm_cnv > 0){
             for(int i = 0; i < NUM_CHR; i++){
                 parm_freq[i] = (double) parm_freq[i] / num_parm_cnv;
                 // cout << "\tchr " << i+1 << " with arm-level freq " << arm_freq[i] << endl;
             }
         }
         if(num_qarm_cnv > 0){
             for(int i = 0; i < NUM_CHR; i++){
                 qarm_freq[i] = (double) qarm_freq[i] / num_qarm_cnv;
                 // cout << "\tchr " << i+1 << " with arm-level freq " << arm_freq[i] << endl;
             }
         }
         // compute the frequencies of chr-level events by chromosome
         map<int, double> chr_freq;
         for(int i = 0; i < NUM_CHR; i++){
             chr_freq[i] = 0.0;
         }
         int num_chr_cnv = 0;
         for(unsigned int i = 0; i < cells.size(); i++) {
             Cell cell = cells[i];
             cell.set_obs_cn();
             // cout << "cell " << cell.cell_ID << endl;
             for(auto cp : cell.obs_cn_profile){
                 // cout << cp.first.first << "\t" << cp.first.second << "\t" << cp.second << endl;
                 if(cp.first.second == 0 && cp.second != 0){
                     num_chr_cnv += 1;
                     chr_freq[cp.first.first] += 1;
                 }
             }
         }
         // cout << "There are " << num_chr_cnv << " chr-level events" << endl;
         if(num_chr_cnv > 0){
             for(int i = 0; i < NUM_CHR; i++){
                 chr_freq[i] = (double) chr_freq[i] / num_chr_cnv;
                 // cout << "\tchr " << i+1 << " with chr-level freq " << chr_freq[i] << endl;
             }
         }

         ofstream fout(fname);
         for(int i = 0; i < NUM_CHR; i++){
             fout << chr_freq[i] << endl;
         }
         for(int i = 0; i < NUM_CHR; i++){
             fout << parm_freq[i] << endl;
         }
         for(int i = 0; i < NUM_CHR; i++){
             fout << qarm_freq[i] << endl;
         }
         fout.close();
     }

     // get unique CNVs
     void set_uniq_cn(vector<Cell> cells){
         for(int c = 0; c < cells.size(); c++){
             cells[c].set_obs_cn();
             map<pair<int, int>, int> obs_cns = cells[c].obs_cn_profile;
             for (auto cp : obs_cns){
                 tuple<int, int, int> cnvec(cp.first.first, cp.first.second, cp.second);
                 uniq_cns.insert(cnvec);
             }
         }

         // cout << "There are " << uniq_cns.size() << " unique CNA events in the clone " << endl;
         // for(auto cv : uniq_cns){
         //    cout << get<0>(cv) << "\t"  << get<1>(cv) << "\t"  << get<2>(cv) << endl;
         // }
     }

     // print average of absolute relative changes
     void write_avg_cn(vector<Cell> cells, string fname){
         vector<float> avg_cn(NUM_CHR * 3, 0);
         // int n = 0;

         for(int c = 0; c < cells.size(); c++){
             cells[c].set_cn_all();
         }
         for(int i = 0; i < cells.size(); i++){
             // for(int j = i+1; j < cells.size(); j++){
             //     n++;
             //     // assert(cells[i].cn_all.size() == cells[j].cn_all.size());
             //     for(int k = 0; k < cells[i].cn_all.size(); k++){
             //         sum_cn[k] += cells[i].cn_all[k] + cells[j].cn_all[k];   // reciprocal events will be cancelled out
             //     }
             // }
             for(int k = 0; k < cells[i].cn_all.size(); k++){
                 avg_cn[k] += abs(cells[i].cn_all[k]);
             }
         }

         // cout << "There are " << n << " combinations of cells" << endl;
         for(int i = 0; i < avg_cn.size(); i++){
             // avg_cn[i] = avg_cn[i] / n;
             avg_cn[i] = avg_cn[i] / cells.size();
         }

         ofstream fout(fname);
         for(int i = 0; i < avg_cn.size(); i++){
             fout << avg_cn[i] << endl;
         }
         fout.close();
     }

     // print average of unique absolute relative changes
     void write_avg_uniq_cn(vector<Cell> cells, string fname){
         map<pair<int, int>, int> cp;
         vector<float> avg_cn(NUM_CHR * 3, 0);

         set_uniq_cn(cells);

         for(int i = 0; i < NUM_CHR; i++){
             for (int j = 0; j < 3; j++){
                 pair<int, int> pos(i,j);
                 cp[pos] = 0;
                 avg_cn[i * 3 + j] = 0;
             }
         }

         for(auto cv: uniq_cns){
             pair<int, int> pos(get<0>(cv), get<1>(cv));
             cp[pos] += abs(get<2>(cv));
             avg_cn[get<0>(cv) * 3 + get<1>(cv)] += abs(get<2>(cv));
         }

         // cout << "number of CNA changes at each position" << endl;
         // for(auto c : cp){
         //     if(c.second!=0)
         //        cout << c.first.first << "\t"  << c.first.second << "\t" << c.second << endl;
         // }

         for(int i = 0; i < avg_cn.size(); i++){
             avg_cn[i] = avg_cn[i] / cells.size();
         }

         ofstream fout(fname);
         for(int i = 0; i < avg_cn.size(); i++){
             fout << avg_cn[i] << endl;
         }
         fout.close();
     }

     // print average of unique absolute relative changes by ignoring chromosome information
     void write_aggregated_uniq_cn(vector<Cell> cells, string fname){
         vector<float> avg_cn(2, 0);    // only count chr and arm level

         set_uniq_cn(cells);

         for (int j = 0; j < 2; j++){
             avg_cn[j] = 0;
         }

         for(auto cv: uniq_cns){
             if(get<1>(cv)==0){
                 avg_cn[0] += abs(get<2>(cv));
             }else{
                 avg_cn[1] += abs(get<2>(cv));
             }
         }

         // cout << "There are " << n << " combinations of cells" << endl;
         for(int i = 0; i < avg_cn.size(); i++){
             avg_cn[i] = avg_cn[i] / cells.size();
         }

         ofstream fout(fname);
         for(int i = 0; i < avg_cn.size(); i++){
             fout << avg_cn[i] << endl;
         }
         fout.close();
     }


     void get_avg_reciprocal_cn(vector<Cell> cells, vector<int>& avg_cn, bool adjust_bias=false){
         map<pair<int, int>, int> cp;
         // vector<int> avg_cn(2, 0);    // only count chr and arm level

         set_uniq_cn(cells);

         for(int i = 0; i < NUM_CHR; i++){
             for (int j = 0; j < 3; j++){
                 pair<int, int> pos(i,j);
                 cp[pos] = 0;
             }
         }

         for (int j = 0; j < 2; j++){
             avg_cn[j] = 0;
         }

         for(auto cv: uniq_cns){
             pair<int, int> pos(get<0>(cv), get<1>(cv));
             cp[pos] += abs(get<2>(cv));
         }

         // sum of reciprocal events will have sum being even
         // cout << "number of absolute CNA changes at each position (adjusting for odd arm-level events)" << endl;
         for(auto cv: cp){
             if(adjust_bias){
                 if(cv.second % 2 !=0)  // probably caused by further arm-level events
                 {
                     // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\t" << cv.second << endl;
                     // arm- and chr- level have different bias
                     if(cv.first.second==0){
                         cv.second = cv.second + 1;
                     }
                     else{
                        cv.second = cv.second - 1;
                     }
                     // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\tafter" << cv.second << endl;
                 }
             }
             // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\t" << cv.second << endl;
             if(cv.first.second==0){
                 avg_cn[0] += cv.second;
             }else{
                 avg_cn[1] += cv.second;
             }
         }

         // return avg_cn;
     }


     // print average of reciprocal CN changes by ignoring chromosome information
     void write_avg_reciprocal_cn(vector<Cell> cells, string fname, bool adjust_bias=false){
         vector<int> avg_cn(2, 0);
         get_avg_reciprocal_cn(cells, avg_cn, adjust_bias);

         ofstream fout(fname);
         for(int i = 0; i < avg_cn.size(); i++){
             fout << (double) avg_cn[i] / cells.size() << endl;
         }
         fout.close();
     }

     // print odd ratio of reciprocal CN changes by ignoring chromosome information
     void write_relative_reciprocal_cn(vector<Cell> cells, string fname, bool adjust_bias=false){
         map<pair<int, int>, int> cp;
         vector<int> avg_cn(2, 0);    // only count chr and arm level

         set_uniq_cn(cells);

         for(int i = 0; i < NUM_CHR; i++){
             for (int j = 0; j < 3; j++){
                 pair<int, int> pos(i,j);
                 cp[pos] = 0;
             }
         }

         for (int j = 0; j < 2; j++){
             avg_cn[j] = 0;
         }

         for(auto cv: uniq_cns){
             pair<int, int> pos(get<0>(cv), get<1>(cv));
             cp[pos] += abs(get<2>(cv));
         }

         // sum of reciprocal events will have sum being even
         // cout << "number of absolute CNA changes at each position (adjusting for odd arm-level events)" << endl;
         for(auto cv: cp){
             if(adjust_bias){
                 if(cv.second % 2 !=0)  // probably caused by further arm-level events
                 {
                     // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\t" << cv.second << endl;
                     // arm- and chr- level have different bias
                     if(cv.first.second==0){
                         cv.second = cv.second + 1;
                     }
                     else{
                        cv.second = cv.second - 1;
                     }
                     // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\tafter" << cv.second << endl;
                 }
             }
             // cout << cv.first.first + 1 << "\t"  << cv.first.second << "\t" << cv.second << endl;
             if(cv.first.second==0){
                 avg_cn[0] += cv.second;
             }else{
                 avg_cn[1] += cv.second;
             }
         }

         double odd = (double) avg_cn[0] / (avg_cn[0] + avg_cn[1]);
         ofstream fout(fname);
         // for(int i = 0; i < avg_cn.size(); i++){
         //     fout << (double) avg_cn[i] / cells.size() << endl;
         // }
         fout << odd << endl;
         fout.close();
     }

     void write_summary_stats(vector<double>& sum_stats, string fname){
         ofstream fout(fname);
         fout.precision(9);
         for(int i = 0; i < sum_stats.size(); i++){
             fout << sum_stats[i] << endl;
         }
         fout.close();
     }

};

#endif
