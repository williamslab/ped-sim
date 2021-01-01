// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <vector>
#include "geneticmap.h"

// Read in genetic map from <mapFile>. Also determines whether there are male
// and female maps present and sets <sexSpecificMaps> to true if so. If only
// one map is present, this is assumed to be the sex-averaged map and in that
// case, sets <sexSpecificMap> to false.
GeneticMap::GeneticMap(char *mapFile, bool &sexSpecificMaps) {
  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  const char *delim = " \t\n";

  FILE *in = fopen(mapFile, "r");
  if (!in) {
    printf("ERROR: could not open map file %s!\n", mapFile);
    perror("open");
    exit(1);
  }

  char *curChr = NULL, *endptr;
  int prevPhysPos = -1;
  vector<PhysGeneticPos> *curMap = NULL;
  sexSpecificMaps = false; // will be updated on first pass below

  while (getline(&buffer, &bytesRead, in) >= 0) {
    char *chrom, *physPosStr, *mapPos1Str, *mapPos2Str;
    char *saveptr;

    if (buffer[0] == '#')
      continue; // comment

    // get all the tokens:
    chrom = strtok_r(buffer, delim, &saveptr);
    physPosStr = strtok_r(NULL, delim, &saveptr);
    mapPos1Str = strtok_r(NULL, delim, &saveptr);
    mapPos2Str = strtok_r(NULL, delim, &saveptr);

    if (curChr == NULL && mapPos2Str != NULL)
      sexSpecificMaps = true;

    // need a new entry in <map> for a new chrom?
    if (curChr == NULL || strcmp(chrom, curChr) != 0) {
      curChr = new char[ strlen(chrom) + 1 ];
      if (curChr == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
      strcpy(curChr, chrom);
      curMap = new vector<PhysGeneticPos>;
      if (curMap == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
      map.emplace_back(curChr, curMap);
      prevPhysPos = -1;
    }

    int physPos;
    double mapPos1, mapPos2 = 0.0;
    errno = 0; // initially
    physPos = strtol(physPosStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: could not parse column 2 of map file as integer\n");
      if (errno != 0)
	perror("strtol");
      exit(2);
    }
    if (physPos <= prevPhysPos) {
      fprintf(stderr, "ERROR: column 2 of map file is not sorted\n");
      exit(2);
    }
    mapPos1 = strtod(mapPos1Str, &endptr);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: could not parse column 3 of map file as floating point\n");
      if (errno != 0)
	perror("strtod");
      exit(2);
    }
    if (sexSpecificMaps) {
      mapPos2 = strtod(mapPos2Str, &endptr);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: could not parse column 4 of map file as floating point\n");
	if (errno != 0)
	  perror("strtod");
	exit(2);
      }
    }
    else if (mapPos2Str != NULL) {
      fprintf(stderr, "ERROR: expected three columns on all lines in map file but more seen\n");
      exit(2);
    }

    curMap->emplace_back(physPos, mapPos1, mapPos2);
    prevPhysPos = physPos;
  }

  free(buffer);
  fclose(in);
}
