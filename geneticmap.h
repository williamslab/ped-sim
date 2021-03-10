// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include "cmdlineopts.h"

#ifndef GENETICMAP_H
#define GENETICMAP_H

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Used to store the genetic map
struct PhysGeneticPos {
  PhysGeneticPos(int p, double m1, double m2) {
    physPos = p; mapPos[0] = m1; mapPos[1] = m2;
  }
  int physPos; double mapPos[2];
};


class GeneticMap {
  public:
    GeneticMap(char *mapFile, bool &sexSpecificMaps);

    size_t size() { return map.size(); }
    const char * chromName(int chrIdx) { return map[chrIdx].first; }
    bool isX(int chrIdx) {
      if (strcmp( chromName(chrIdx), CmdLineOpts::chrX ) == 0)
	return true;
      else
	return false;
    }
    size_t chromNumPos(int chrIdx) { return map[chrIdx].second->size(); }

    int chromStartPhys(int chrIdx) {
      return map[chrIdx].second->front().physPos;
    }
    int chromEndPhys(int chrIdx) {
      return map[chrIdx].second->back().physPos;
    }
    int chromPhysPos(int chrIdx, int entry) {
      return (*map[chrIdx].second)[entry].physPos;
    }

    double chromStartGenet(int chrIdx, int sex) {
      return map[chrIdx].second->front().mapPos[sex];
    }
    double chromEndGenet(int chrIdx, int sex) {
      return map[chrIdx].second->back().mapPos[sex];
    }
    double chromGenetPos(int chrIdx, int sex, int entry) {
      return (*map[chrIdx].second)[entry].mapPos[sex];
    }
    double chromGenetLength(int chrIdx, int sex) {
      return chromEndGenet(chrIdx, sex) - chromStartGenet(chrIdx, sex);
    }

    bool haveXmap() {
      for(size_t i = 0; i < map.size(); i++) {
	if (isX(i)) {
	  return true;
	}
      }
      return false;
    }

  private:
    GeneticMap() { }; // disallow default constructor

    vector< pair<char*, vector<PhysGeneticPos>* > > map;
};

#endif // GENETICMAP_H
