// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include "geneticmap.h"

#ifndef FIXEDCOS_H
#define FIXEDCOS_H

#ifndef NOFIXEDCO // allow disabling of fixed CO functionality at compile time
struct FixedCOs {
  static void read(const char *fixedCOfile, GeneticMap &map);

  static vector<int> & getCOs(int sex, unsigned int idx, unsigned int chrIdx) {
    return theCOs[sex][idx][chrIdx];
  }

  // Array of size 2: one set for males, one for females
  // Outer vector is for individuals, next level for chromosomes, and the last
  // stores a vector of positions of COs
  static vector< vector< vector<int> > > theCOs[2];
};
#endif // NOFIXEDCO


#endif // FIXEDCOS_H
