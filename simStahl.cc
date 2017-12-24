// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

// This code was ported from the xoi R package written by Karl Broman and
// implemented in C (https://github.com/kbroman/xoi)

#include <math.h>
#include <random>
#include <algorithm>
#include <boost/math/special_functions/gamma.hpp>

// Note: we don't use GSL version of gamma cdf to reduce library dependencies,
// but it is almost 2x faster than boost
//#include <gsl/gsl_sf_gamma.h>


using namespace std;

uniform_real_distribution<double> unif_prob(0.0, 1.0);
extern uniform_int_distribution<int> coinFlip;

// locations: stores sampled crossover locations (assumed empty initially)
// sex: 0 or 1 for male or female meiosis
// chr: chromosome number subtracted by 1 so 1-22 maps to 0-21
// randomGen: random number generator
// n_bins4start: number of bins for approximating cdf of first crossover
void simStahl(vector<double> &locations, double nu, double p, double length,
	      mt19937 &randomGen, int n_bins4start) {
  if (fabs(nu - 1.0) < 1e-8) { // looks like a Poisson model
    poisson_distribution<int> standard(length);
    int nxo = standard(randomGen);

    for(int j = 0; j < nxo; j++)
      locations.push_back( unif_prob(randomGen) * length );

    sort(locations.begin(), locations.end());

    return; // done; below code is for mixture model
  }

  /////////////////////////////////////////////////////////////////////////////
  // using mixture model

  // original code used scale = 1 / rate, but rate is more convenient for the
  // libraries we use
  double rate = (2.0 * nu * (1.0 - p));

  // approximate distribution to first (start) crossover using
  // <n_bins4start> bins
  double *startprob = new double[n_bins4start];
  double step = length / n_bins4start;

  // Using R library:
//  startprob[0] = 2.0 * (1.0 - p) * pgamma(0.5*step, nu, 1/rate, 0, 0) * step;
  // Using Boost (and GSL commented out):
  startprob[0] = 2.0 * (1.0 - p) *
		  boost::math::gamma_q(/*shape=*/ nu, 0.5*step * rate) * step;
//		  gsl_sf_gamma_inc_Q(/*shape=*/nu, 0.5*step * rate) * step;
  for(int i = 1; i < n_bins4start; i++) {
    startprob[i] = startprob[i-1] +
//	2.0*(1.0 - p) * pgamma((i + 0.5)*step, nu, 1/rate, 0, 0) * step;
	2.0*(1.0 - p) *
	  boost::math::gamma_q(/*shape=*/ nu, (i + 0.5)*step * rate) * step;
//	  gsl_sf_gamma_inc_Q(/*shape=*/ nu, (i+0.5)*step * rate) * step;
  }

  // now sample location of current crossover -- initially none: so chr start
  double curloc = 0.0;

  // locations of chiasmata from the gamma model
  // shape = nu, rate = 2*nu*(1-p) [scale = 1/{2*nu*(1-p)}]

  double u = unif_prob(randomGen);
  if ( u > startprob[ n_bins4start - 1 ] )
    curloc = length+1; // no crossovers: at end of chromosome
  else {
    // TODO: should switch to binary search
    for(int j = 0; j < n_bins4start; j++) {
      if(u <= startprob[j]) {
	curloc = ((double) j + 0.5) * step;
	if(coinFlip(randomGen)) // on this chromatid? coin toss
	  locations.push_back(curloc);

	break;
      }
    }
  }

  gamma_distribution<double> gammaRand(nu, /*scale=*/ 1 / rate);
  while(curloc < length) {
    curloc += gammaRand(randomGen); // location of next chiasmata
    // is it before the end of the chromosome and on this chromatid?
    if(curloc < length && coinFlip(randomGen)) // coin toss to decide chromatid
      locations.push_back(curloc);
  }

  // locations of crossovers from the no interference mechanism
  if(p > 0) {
    poisson_distribution<int> no_interfere(length * p);
    int n_nixo = no_interfere(randomGen);

    for(int j = 0; j < n_nixo; j++)
      // sample position of the non-interference derived crossovers uniformly
      locations.push_back( unif_prob(randomGen) * length );

    // Same as above but drawing distance to next event from exponential:
//    double curCOPos = 0.0;
//    while (curCOPos < length) {
//      double distToNext = crossoverDist(randomGen) / 100;
//      curCOPos += distToNext;
//      if (curCOPos < length && unif_prob(randomGen) < p) {
//	locations.push_back(curCOPos);
//      }
//    }
  }

  sort(locations.begin(), locations.end());

  delete [] startprob;
}
