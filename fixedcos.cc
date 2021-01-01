// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <assert.h>
#include "fixedcos.h"

#ifndef NOFIXEDCO // allow disabling of fixed CO functionality at compile time
vector< vector< vector<int> > > FixedCOs::theCOs[2];

// Given a file <fixedCOfile> containing sets of crossovers labeled paternal and
// maternal, reads in the crossovers for use in simulating
void FixedCOs::read(const char *fixedCOfile, GeneticMap &map) {
  size_t bytesRead = 1024, bytesReadOther = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  char *bufferOther = (char *) malloc(bytesRead + 1);
  if (buffer == NULL || bufferOther == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  const char *delim = " \t\n";

  FILE *in = fopen(fixedCOfile, "r");
  if (!in) {
    fprintf(stderr, "ERROR: could not open fixed crossover file %s!\n", fixedCOfile);
    exit(1);
  }

  char *lastId = NULL;
  int chrIdx = -1; // index of current chromosome
  int lastPM_idx = -1;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    char *id, *pat_mat, *chrom, *posStr;
    char *saveptr, *endptr;

    // get sample id, paternal/maternal type
    id = strtok_r(buffer, delim, &saveptr);
    pat_mat = strtok_r(NULL, delim, &saveptr);

    int pm_idx = (pat_mat[0] == 'M') ? 1 : 0; // 1 for maternal, 0 for paternal

    if (pm_idx != lastPM_idx || lastId == NULL || strcmp(lastId, id) != 0) {
      if (lastId) {
	while (chrIdx < (int) map.size() - 1) {
	  // new (empty) chr:
	  theCOs[lastPM_idx].back().emplace_back();
	  chrIdx++;
	}
      }
      // new person:
      theCOs[pm_idx].emplace_back();
      assert(theCOs[pm_idx].size() < UINT_MAX);
      chrIdx = -1; // reset
    }

    // get chromosome
    chrom = strtok_r(NULL, delim, &saveptr);

    while (chrIdx == -1 || strcmp(map.chromName(chrIdx), chrom) != 0) {
      // new chr:
      theCOs[pm_idx].back().emplace_back();
      chrIdx++;
      if (chrIdx > (int) map.size()) {
	fprintf(stderr, "ERROR: chromosome %s in %s does not match any chromosome in genetic map\n",
		chrom, fixedCOfile);
      }
    }

    // get position
    posStr = strtok_r(NULL, delim, &saveptr);
    errno = 0; // initially
    int pos = strtol(posStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: column three of %s contains %s, which cannot be converted to an integer\n",
	  fixedCOfile, posStr);
      exit(2);
    }

    if (pos >= map.chromStartPhys(chrIdx) && pos <= map.chromEndPhys(chrIdx))
      // add to COs so long as the position is within the range of the genetic
      // map
      theCOs[pm_idx].back().back().push_back( pos );

    // Swap values for next round:
    lastId = id;
    lastPM_idx = pm_idx;

    char *tmpBuf = buffer;
    buffer = bufferOther;
    bufferOther = tmpBuf;
    size_t tmpBytesRead = bytesRead;
    bytesRead = bytesReadOther;
    bytesReadOther = tmpBytesRead;
  }

  fclose(in);
}

#endif // NOFIXEDCO
