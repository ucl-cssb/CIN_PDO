#ifndef STATS_HPP
#define STATS_HPP


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_multimin.h>

#include <cstdio>
#include <algorithm>
#include <vector>
#include <map>
#include <cmath>
#include <set>
#include <climits>
#include <cstdlib>

#include <unistd.h>

gsl_rng * r;
std::mt19937 eng;

unsigned long setup_rng(unsigned long set_seed){
  gsl_rng_env_setup();

  const gsl_rng_type* T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  if( set_seed != 0 ){
    gsl_rng_set(r, set_seed);
    return(set_seed);
  }else{
    int t = time(NULL);
    int pid = getpid();
    long s = t*pid;
    //cout << "pid:" << "\t" << getpid() << endl;
    // cout << "seed:" << "\t" << t << "\t" << pid << "\t" << abs(s) << endl;
    std::cout << "seed:" << "\t" << abs(s) << std::endl;
    gsl_rng_set (r, abs(s));
    return(abs(s));
  }
}

// unary function and pointer to unary function
// allows use of gsl rng for standard template algorithms
long unsigned myrng(long unsigned n)
{
  return gsl_rng_uniform_int(r, n);
}

long unsigned (*fp)(long unsigned) = myrng;

// factorial for choose(n,k)
int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}

// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b-a)*gsl_rng_uniform (r);

  while(myrandom==0){
    myrandom = a + (b-a)*gsl_rng_uniform (r);
  }
  return myrandom;
}


// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const std::vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i=0; i<rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  std::vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i=0; i<rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  std::vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i=0; i<rates.size(); ++i){
    if( u < psum[i] ){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}


#endif
