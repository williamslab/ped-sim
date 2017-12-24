// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>

#ifndef SIMSTAHL_H
#define SIMSTAHL_H

using namespace std;

void simStahl(vector<double> &locations, int sex, int chr, mt19937 &randomGen,
	      int n_bins4start = 10000);

#endif // SIMSTAHL_H
