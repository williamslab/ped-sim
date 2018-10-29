// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

// This code was ported from the xoi R package written by Karl Broman and
// implemented in C (https://github.com/kbroman/xoi)

#include <math.h>
#include <random>
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>
#include "cointerfere.h"

// Note: we don't use GSL version of gamma cdf to reduce library dependencies,
// but it is almost 2x faster than boost
//#include <gsl/gsl_sf_gamma.h>

uniform_real_distribution<double> unif_prob(0.0, 1.0);
extern uniform_int_distribution<int> coinFlip;


// approximate distribution to first (start) crossover using <N_BINS4START>
// bins
// store away for computational efficiency
void COInterfere::initStartProb() {
  for(int s = 0; s < 2; s++) {
    // original code used scale = 1 / rate, but rate is more convenient for the
    // libraries we use
    double rate = (2.0 * nu[s] * (1.0 - p[s]));

    double step = length[s] / N_BINS4START;

    // Using R library:
//    startProb[s][0] = 2.0 * (1.0 - p[s]) * pgamma(0.5*step, nu[s], 1/rate, 0, 0) * step;
    // Using Boost (and GSL commented out):
    startProb[s][0] = 2.0 * (1.0 - p[s]) *
		  boost::math::gamma_q(/*shape=*/nu[s], 0.5*step * rate) * step;
//		  gsl_sf_gamma_inc_Q(/*shape=*/nu[s], 0.5*step * rate) * step;
    for(int i = 1; i < N_BINS4START; i++) {
      startProb[s][i] = startProb[s][i-1] +
//	2.0*(1.0 - p[s]) * pgamma((i + 0.5)*step, nu[s], 1/rate, 0, 0) * step;
	2.0*(1.0 - p[s]) *
	    boost::math::gamma_q(/*shape=*/nu[s], (i + 0.5)*step * rate) * step;
//	    gsl_sf_gamma_inc_Q(/*shape=*/nu[s], (i+0.5)*step * rate) * step;
    }
  }
}

// locations: stores sampled crossover locations (assumed empty initially)
// sex: 0 or 1 for male or female meiosis, respectively
// randomGen: random number generator
void COInterfere::simStahl(vector<double> &locations, int sex,
			   mt19937 &randomGen) {
  if (fabs(nu[sex] - 1.0) < 1e-8) { // looks like a Poisson model
    poisson_distribution<int> standard(length[sex]);
    int nxo = standard(randomGen);

    for(int j = 0; j < nxo; j++)
      locations.push_back( unif_prob(randomGen) * length[sex] );

    sort(locations.begin(), locations.end());

    return; // done; below code is for mixture model
  }

  /////////////////////////////////////////////////////////////////////////////
  // using mixture model

  // original code used scale = 1 / rate, but rate is more convenient for the
  // libraries we use
  double rate = (2.0 * nu[sex] * (1.0 - p[sex]));

  double step = length[sex] / N_BINS4START;

  // sample location of current crossover -- initially none: so chr start
  double curloc = 0.0;

  // locations of chiasmata from the gamma model
  // shape = nu, rate = 2*nu*(1-p) [scale = 1/{2*nu*(1-p)}]

  double u = unif_prob(randomGen);
  if ( u > startProb[sex][ N_BINS4START - 1 ] )
    curloc = length[sex]+1; // no crossovers: at end of chromosome
  else {
    // faster binary search:
    if (u <= startProb[sex][0])
      // corner case that doesn't work well with binary search
      curloc = 0.5 * step;
    else {
      int low = 0, high = N_BINS4START;
      while (high - low > 1) {
	int mid = (high - low) / 2 + low;
	if (u <= startProb[sex][mid])
	  high = mid;
	else if (u > startProb[sex][mid])
	  low = mid;
      }
      curloc = ((double) high + 0.5) * step;
    }
    if(coinFlip(randomGen)) // on this chromatid? coin toss
      locations.push_back(curloc);

    // original linear search:
//    for(int j = 0; j < N_BINS4START; j++) {
//      if(u <= startProb[sex][j]) {
//	curloc = ((double) j + 0.5) * step;
//	if(coinFlip(randomGen)) // on this chromatid? coin toss
//	  locations.push_back(curloc);
//
//	break;
//      }
//    }
  }

  gamma_distribution<double> gammaRand(nu[sex], /*scale=*/ 1 / rate);
  while(curloc < length[sex]) {
    curloc += gammaRand(randomGen); // location of next chiasmata
    // is it before the end of the chromosome and on this chromatid?
    if(curloc < length[sex] && coinFlip(randomGen))
      // coin toss decides chromatid
      locations.push_back(curloc);
  }

  // locations of crossovers from the no interference mechanism
  if(p[sex] > 0) {
    poisson_distribution<int> no_interfere(length[sex] * p[sex]);
    int n_nixo = no_interfere(randomGen);

    for(int j = 0; j < n_nixo; j++)
      // sample position of the non-interference derived crossovers uniformly
      locations.push_back( unif_prob(randomGen) * length[sex] );

    // Same as above but drawing distance to next event from exponential:
//    double curCOPos = 0.0;
//    while (curCOPos < length[sex]) {
//      double distToNext = crossoverDist(randomGen);
//      curCOPos += distToNext;
//      if (curCOPos < length[sex] && unif_prob(randomGen) < p[sex]) {
//	locations.push_back(curCOPos);
//      }
//    }
  }

  sort(locations.begin(), locations.end());
}
