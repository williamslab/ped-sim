#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <random>
#include <chrono>
#include <assert.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Used to store necessary details about the source and length of each segment:
//
// note: start marker is implicit
struct Segment { int foundHapNum, endPos; };
typedef vector<Segment> Haplotype;
struct Person {
  Person() { sex = 0; } // by default assume using sex-averaged map: all 0
  int sex;
  vector<Haplotype> haps[2]; // haplotype pair for <this>
};

////////////////////////////////////////////////////////////////////////////////
// Used to store the genetic map
struct PhysGeneticPos {
  PhysGeneticPos(int p, double m1, double m2) {
    physPos = p; mapPos[0] = m1; mapPos[1] = m2;
  }
  int physPos; double mapPos[2];
};

////////////////////////////////////////////////////////////////////////////////
// Function decls
void readDat(int *numSampsToRetain, int numGen, char *datFile);
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps);
void simulate(Person ***generations[2], int *numSampsToRetain, int numGen,
	      int numFam,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps);
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap);
void printBPs(Person ***generations[2], int *numSampsToRetain, int numGen,
	      int numFam,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *outFile);
void makeVCF(Person ***generations[2], int *numSampsToRetain, int numGen,
	     int numFam, char *inVCFfile, char *outVCFfile,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap);
template<typename T>
void pop_front(vector<T> &vec);
void printUsage(char **argv);

mt19937 randomGen;
uniform_int_distribution<int> coinFlip(0,1);
exponential_distribution<double> crossoverDist(1.0 / 100); // in cM units

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 8) {
    printUsage(argv);
  }

  // seed random number generator:
  unsigned int seed = random_device().entropy();
  if (seed == 0 && random_device().entropy() == 0) {
    // random_device is not a real random number generator: fall back on
    // using time to generate a seed:
    seed = chrono::system_clock::now().time_since_epoch().count();
  }
//  seed = 695558636u; // for testing
  printf("Using random seed: %u\n\n", seed);
  randomGen.seed(seed);

  int numFam = atoi(argv[1]);

  int numGen = atoi(argv[2]);
  // We'll save generated haplotypes for the indicated number of samples in the
  // indicated generations. This information is indicated in the dat file read
  // below. But we initially assume that we'll save 0
  int *numSampsToRetain = new int[numGen];
  for(int i = 0; i < numGen; i++)
    numSampsToRetain[i] = 0;

  readDat(numSampsToRetain, numGen, /*datFile=*/ argv[3]);

  vector< pair<char*, vector<PhysGeneticPos>* > > geneticMap;
  bool sexSpecificMaps;
  readMap(geneticMap, /*mapFile=*/ argv[4], sexSpecificMaps);

  // Besides the top-most generation, there are two sides to the pedigree, with
  // element 0 in <generation> containing data from one side and element 1 the
  // other.
  // The second index (e.g., <g> in <generations[side][g]>) is the generation
  // number.
  // The third index is the sample number in the given generation
  // The fourth index is the haplotype 0 or 1 of the given sample
  Person ***generations[2];

  printf("Simulating... ");
  fflush(stdout);
  simulate(generations, numSampsToRetain, numGen, numFam, geneticMap,
	   sexSpecificMaps);
  printf("done.\n");

  printf("Printing break points... ");
  fflush(stdout);
  printBPs(generations, numSampsToRetain, numGen, numFam, geneticMap,
	   /*outFile=*/ argv[7]);
  printf("done.\n");

  printf("Generating VCF file... ");
  fflush(stdout);
  makeVCF(generations, numSampsToRetain, numGen, numFam, /*inVCFfile=*/ argv[5],
	  /*outVCFfile=*/ argv[6], geneticMap);
  printf("done.\n");

  return 0;
}

// Reads the number of samples to retain in each generation from the dat file
void readDat(int *numSampsToRetain, int numGen, char *datFile) {
  // open dat file:
  FILE *in = fopen(datFile, "r");
  if (!in) {
    printf("ERROR: could not open dat file %s!\n", datFile);
    exit(1);
  }

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *delim = " \t\n";

  bool modifiedNumToRetain = false;
  int line = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    line++;

    char *genNumStr, *numSampsStr, *saveptr;
    genNumStr = strtok_r(buffer, delim, &saveptr);
    numSampsStr = strtok_r(NULL, delim, &saveptr);

    if (genNumStr == NULL || numSampsStr == NULL) {
      printf("ERROR: line number %d in dat file does not contain two fields\n",
	      line);
    }

    int generation = atoi(genNumStr);
    int numSamps = atoi(numSampsStr);

    if (generation < 2 || generation > numGen) { // TODO: document
      fprintf(stderr, "ERROR: generation %d in dat file below 2 or above %d (max number of generations)\n",
	      generation, numGen);
      exit(1);
    }
    if (numSamps <= 0) {
      fprintf(stderr, "ERROR: in generation %d, number of samples to simulate below 0\n",
	      generation);
      exit(2);
    }

    // subtract 1 because array is 0 based
    if (numSampsToRetain[generation - 1] != 0) {
      fprintf(stderr, "ERROR: multiple entries for generation %d in dat file\n",
	      generation);
    }
    numSampsToRetain[generation - 1] = numSamps;
    modifiedNumToRetain = true;
  }

  if (!modifiedNumToRetain) {
    fprintf(stderr, "ERROR: need dat file with information about numbers of samples to retain\n");
    exit(3);
  }

  if (numSampsToRetain[numGen-1] == 0) { // TODO: document
    fprintf(stderr, "ERROR: request to simulate %d generations; must print data from last generation\n",
	    numGen);
    fprintf(stderr, "       modify either number of generations or dat file\n");
    exit(4);
  }

  fclose(in);
}

// Read in genetic map from <mapFile> into <geneticMap>. Also determines whether
// there are male and female maps present and sets <sexSpecificMaps> to true if
// so. If only one map is present, this is assumed to be the sex-averaged map
// and in that case, sets <sexSpecificMap> to false.
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps) {
  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *delim = " \t\n";

  FILE *in = fopen(mapFile, "r");
  if (!in) {
    printf("ERROR: could not open map file %s!\n", mapFile);
    exit(1);
  }

  char *curChr = NULL;
  vector<PhysGeneticPos> *curMap = NULL;
  sexSpecificMaps = false; // will be updated on first pass below

  while (getline(&buffer, &bytesRead, in) >= 0) {
    char *chrom, *physPosStr, *mapPos1Str, *mapPos2Str;
    char *saveptr;

    // get all the tokens:
    chrom = strtok_r(buffer, delim, &saveptr);
    if (chrom[0] == '#')
      continue; // comment
    physPosStr = strtok_r(NULL, delim, &saveptr);
    mapPos1Str = strtok_r(NULL, delim, &saveptr);
    mapPos2Str = strtok_r(NULL, delim, &saveptr);

    if (curChr == NULL && mapPos2Str != NULL)
      sexSpecificMaps = true;

    // need a new entry in <geneticMap> for a new chrom?
    // TODO: document that maps must be in order in terms of chroms and
    // physical position
    if (curChr == NULL || strcmp(chrom, curChr) != 0) {
      curChr = (char *) malloc(strlen(chrom) + 1);
      strcpy(curChr, chrom);
      curMap = new vector<PhysGeneticPos>;
      geneticMap.emplace_back(curChr, curMap);
    }

    int physPos;
    double mapPos1, mapPos2 = 0.0;
    physPos = atoi(physPosStr);
    mapPos1 = atof(mapPos1Str);
    if (sexSpecificMaps) {
      mapPos2 = atof(mapPos2Str);
    }
    else {
      assert(mapPos2Str == NULL); // TODO: do something smarter here
    }

    curMap->emplace_back(physPos, mapPos1, mapPos2);
  }

  free(buffer);
  fclose(in);
}

// Simulate data for each generation
void simulate(Person ***generations[2], int *numSampsToRetain, int numGen,
	      int numFam,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps) {
  // Note: throughout we use 0-based values for generations though the input
  // dat file is 1-based

  // Make Person objects for all those we will simulate, assigning sex if
  // <sexSpecificMaps> is true.
  // Top-most generation:
  Person dad, mom;
  if (sexSpecificMaps)
    mom.sex = 1;
  // this loop is setup to improve memory locality for each family
  for (int fam = 0; fam < numFam; fam++) {
    for(int side = 0; side < 2; side++) {
      if (fam == 0)
	generations[side] = new Person**[numFam];
      generations[side][fam] = new Person*[numGen];
    }
  }
  for(int side = 0; side < 2; side++) {
    for (int fam = 0; fam < numFam; fam++) {
      for(int curGen = 1; curGen < numGen; curGen++) {
	// Determine how many samples we need data for in <curGen>:
	int numPersons = numSampsToRetain[curGen];
	if (numPersons == 0) // not saving, but need parent of next generation
	  numPersons = 1;
	// additional person that is the other parent of next generation
	numPersons++;

	generations[side][fam][curGen] = new Person[numPersons];

	if (sexSpecificMaps) {
	  // the two individuals that reproduce are index 0 (a founder) and 1
	  // randomly decide which one to make female
	  int theFemale = coinFlip(randomGen);
	  generations[side][fam][curGen][theFemale].sex = 1;
	}
      }
    }
  }

  // have exactly 2 founders / 4 founder haplotypes per generation
  int numFounderHapsPerFam = numGen * 4;

  for(int fam = 0; fam < numFam; fam++) { // simulate each family

    for (int h = 0; h < 2; h++) { // clear out the old haplotypes for <ad,mom
      dad.haps[h].clear();
      mom.haps[h].clear();
    }

    // simulate haplotypes for each chrom:
    for(auto it = geneticMap.begin(); it != geneticMap.end(); it++) {
      vector<PhysGeneticPos> *curMap = it->second;

      // Make trivial haplotypes for generation 0 founders:
      Segment trivialSeg;
      // no recombinations in founders
      trivialSeg.endPos = curMap->back().physPos;
      for(int h = 0; h < 2; h++) {
	trivialSeg.foundHapNum = fam * numFounderHapsPerFam + h;
	// makes a copy of <trivialSeg>, can reuse
	dad.haps[h].emplace_back();
	dad.haps[h].back().push_back(trivialSeg);
	trivialSeg.foundHapNum = fam * numFounderHapsPerFam + h + 2;
	mom.haps[h].emplace_back();
	mom.haps[h].back().push_back(trivialSeg);
      }

      // do the simulation for each generation
      for(int curGen = 1; curGen < numGen; curGen++) {

	// First make trivial haplotypes for the two founders in <curGen>; use
	// convention that sample 0 is the founder on both sides:
	if (curGen != numGen -1) { //no need for founders in the last generation
	  for(int side = 0; side < 2; side++) {
	    for(int h = 0; h < 2; h++) {
	      // 4 founder haplotypes per generation
	      trivialSeg.foundHapNum = fam * numFounderHapsPerFam +
		curGen * 4 + side * 2 + h;
	      // the following makes a copy of <trivialSeg>, so we can reuse it
	      generations[side][fam][curGen][0].haps[h].emplace_back();
	      generations[side][fam][curGen][0].haps[h].back().push_back(
		  trivialSeg);
	    }
	  }
	}

	int numPersons = numSampsToRetain[curGen];
	if (numPersons == 0)
	  // not saving any, but need parent of next generation
	  numPersons = 1;
	// additional person that is the other parent of next generation
	numPersons++;

	// Now simulate the non-founders in <curGen>
	for(int side = 0; side < 2; side++) {
	  for(int ind = 1; ind < numPersons; ind++) {
	    if (curGen == 1) {
	      // no "sides" for creating generation 1; use <dad> and <mom>
	      generations[side][fam][curGen][ind].haps[0].emplace_back();
	      generations[side][fam][curGen][ind].haps[1].emplace_back();
	      generateHaplotype(
			     generations[side][fam][curGen][ind].haps[0].back(),
			     dad, curMap);
	    generateHaplotype(
			     generations[side][fam][curGen][ind].haps[1].back(),
			     mom, curMap);
	    }
	    else {
	      if (sexSpecificMaps) {
		assert(generations[side][fam][curGen-1][0].sex !=
				      generations[side][fam][curGen-1][1].sex);
	      }
	      for(int parIdx = 0; parIdx < 2; parIdx++) {
		int hapIdx = parIdx;
		if (sexSpecificMaps)
		  hapIdx = generations[side][fam][curGen-1][parIdx].sex;

		generations[side][fam][curGen][ind].haps[hapIdx].emplace_back();
		generateHaplotype(
		    generations[side][fam][curGen][ind].haps[hapIdx].back(),
		    /*parent=*/ generations[side][fam][curGen-1][parIdx],
		    curMap);
	      }
	    }
	  }
	}

      } // <curGen>
    } // <geneticMap> (chroms)
  } // <fam>
}

// Simulate one haplotype <toGenerate> by sampling crossovers and switching
// between the two haplotypes stored in <parent>.
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap) {
  // For the two haplotypes in <parent>, which index is the current
  // <switchMarker> position contained in?
  unsigned int curIndex[2] = { 0, 0 };

  int mapIdx = parent.sex; // either 0 or 1 for sex-averaged/male or female

  // Pick haplotype for the beginning of the transmitted one:
  int curHap = coinFlip(randomGen);
  // centiMorgan position of next crossover -- exponentially random distance
  // from first physical position in the map
  double cMPosNextCO = curMap->front().mapPos[mapIdx] +
						      crossoverDist(randomGen);
  // initially assume we'll recombine between first and positions with map info
  int switchIdx = 1;

  int mapSize = curMap->size();
  double lastcMPos = curMap->back().mapPos[mapIdx];

  while (cMPosNextCO < lastcMPos) {
    // TODO: slow linear search for switching index, but probably fast enough
    for( ; switchIdx < mapSize; switchIdx++) {
      if ((*curMap)[switchIdx].mapPos[mapIdx] > cMPosNextCO)
	break;
    }
    if (switchIdx == mapSize)
      break; // let code below this while loop insert the final segments

    // segment ends between map indexes 1 less than the current <switchIdx>, so:
    switchIdx--;

    // get physical position using linear interpolation:
    double frac = (cMPosNextCO - (*curMap)[switchIdx].mapPos[mapIdx]) /
	      ((*curMap)[switchIdx+1].mapPos[mapIdx] -
					   (*curMap)[switchIdx].mapPos[mapIdx]);
    assert(frac >= 0.0 && frac <= 1.0);
    int switchPos = (*curMap)[switchIdx].physPos +
	frac * ((*curMap)[switchIdx+1].physPos - (*curMap)[switchIdx].physPos);

    // copy Segments from <curHap>
    // TODO: comment about back()
    for( ; curIndex[curHap] < parent.haps[curHap].back().size();
							  curIndex[curHap]++) {
      Segment &seg = parent.haps[curHap].back()[ curIndex[curHap] ];
      if (seg.endPos >= switchPos) {
	// last segment to copy, and we will break it at <switchPos>
	Segment copy = seg; // don't modify <seg.endPos> directly
	copy.endPos = switchPos;
	toGenerate.push_back(copy);
	if (seg.endPos == switchPos)
	  curIndex[curHap]++;
	break; // done copying
      }
      else {
	toGenerate.push_back(seg);
      }
    }
    assert(curIndex[curHap] < parent.haps[curHap].back().size());

    // swap haplotypes
    curHap ^= 1;
    // must update <curIndex[curHap]>
    for( ; curIndex[curHap] < parent.haps[curHap].back().size();
							  curIndex[curHap]++) {
      Segment &seg = parent.haps[curHap].back()[ curIndex[curHap] ];
      if (seg.endPos > switchPos)
	// current segment spans from just after <switchMarker> to <endMarker>
	break;
    }
    assert(curIndex[curHap] < parent.haps[curHap].back().size());

    cMPosNextCO += crossoverDist(randomGen);

    // TODO: implement? document if so
    // find location of next crossover; we use this loop to ensure that the next
    // <switchPos> value is strictly greater than the current one:
    // (Presumably this will get triggered next to never, but just in case. Note
    // that this changes the distribution of the crossover locations, but since
    // there's CO interference in real data, this seems fine.)
//    do {
//      cMPosNextCO += crossoverDist(randomGen);
//    } while (cMPosNextCO < curMap[switchIdx].mapPos[mapIdx]);
  }

  // copy through to the end of the chromosome:
  for( ; curIndex[curHap] < parent.haps[curHap].back().size();
							  curIndex[curHap]++) {
    Segment &seg = parent.haps[curHap].back()[ curIndex[curHap] ];
    toGenerate.push_back(seg);
  }
}

// Print the break points to <outFile> along with meta data for each haplotype
// (sample id with generation number, side of the pedigree, and sample number,
// as well as which haplotype is being printed)
void printBPs(Person ***generations[2], int *numSampsToRetain, int numGen,
	      int numFam,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *outFile) {
  FILE *out = fopen(outFile, "w");
  if (!out) {
    printf("ERROR: could not open output file %s!\n", outFile);
    exit(1);
  }

  assert(numSampsToRetain[0] == 0);
  for(int fam = 0; fam < numFam; fam++) {
    for(int side = 0; side < 2; side++) {
      for(int gen = 1; gen < numGen; gen++) {
	if (numSampsToRetain[gen] > 0) {
	  for(int ind = 1; ind < numSampsToRetain[gen] + 1; ind++) {
	    for(int h = 0; h < 2; h++) {
	      fprintf(out, "f%d_s%d_g%d_i%d h%d", fam+1, side, gen+1, ind, h);

	      for(unsigned int chr = 0; chr < geneticMap.size(); chr++) {
		// print chrom name and starting position
		fprintf(out, " %s|%d", geneticMap[chr].first,
			geneticMap[chr].second->front().physPos);
		Haplotype &curHap = generations[side][fam][gen][ind].
								   haps[h][chr];
		for(unsigned int s = 0; s < curHap.size(); s++) {
		  Segment &seg = curHap[s];
		  fprintf(out, " %d:%d", seg.foundHapNum, seg.endPos);
		}
	      }
	      fprintf(out, "\n");
	    }
	  }
	}
      }
    }
  }

  fclose(out);
}

// Given the simulated break points for individuals in each generation stored in
// <generations> and other necessary information, reads input VCF format data
// from the file named <inVCFfile> and prints the simulated haplotypes for each
// sample to <outVCFfile> in VCF format.
void makeVCF(Person ***generations[2], int *numSampsToRetain, int numGen,
	     int numFam, char *inVCFfile, char *outVCFfile,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap) {
  // open input VCF file:
  FILE *in = fopen(inVCFfile, "r");
  if (!in) {
    printf("ERROR: could not open input VCF file %s!\n", inVCFfile);
    exit(1);
  }

  // open output VCF file:
  FILE *out = fopen(outVCFfile, "w");
  if (!out) {
    printf("ERROR: could not open output VCF file %s!\n", outVCFfile);
    exit(1);
  }

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *tab = "\t";
  const char *bar = "|";
  // Below when we print the VCF output, we alternate printing tab
  // and | between successive alleles. Make this simpler with:
  const char *betweenAlleles[2] = { tab, bar };

  // have exactly 2 founders / 4 founder haplotypes per generation
  int numFounderHaps = numGen * 4 * numFam;
  int *founderHaps = new int[numFounderHaps];

  // iterate over chromosomes in the genetic map
  unsigned int chrIdx = 0; // index of current chromosome number;
  char *chrName = geneticMap[chrIdx].first;
  int chrBegin = geneticMap[chrIdx].second->front().physPos;
  int chrEnd = geneticMap[chrIdx].second->back().physPos;

  bool gotSomeData = false;

  while (getline(&buffer, &bytesRead, in) >= 0) { // read each line of input VCF
    if (buffer[0] == '#' && buffer[1] == '#') {
      // header line: print to output
      fprintf(out, "%s", buffer);
      continue;
    }

    if (buffer[0] == '#') {
      // header line indicating fields and sample ids -- print version for the
      // simulated individuals
      fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

      // sample ids:
      for(int fam = 0; fam < numFam; fam++)
	for(int side = 0; side < 2; side++)
	  for(int gen = 1; gen < numGen; gen++)
	    if (numSampsToRetain[gen] > 0)
	      for(int ind = 1; ind < numSampsToRetain[gen] + 1; ind++)
		fprintf(out, "\tf%d_s%d_g%d_i%d", fam+1, side, gen+1, ind);
      fprintf(out, "\n");
      continue;
    }

    char *saveptr;
    char *chrom = strtok_r(buffer, tab, &saveptr);

    if (strcmp(chrom, chrName) != 0) {
      if (gotSomeData) {
	chrIdx++;
	if (chrIdx == geneticMap.size())
	  // no more chromosomes to process; will ignore remainder of VCF
	  break;
	chrName = geneticMap[chrIdx].first;
      }
      if (!gotSomeData || strcmp(chrom, chrName) != 0) {
	printf("ERROR: chromosome %s in VCF file either out of order or not present\n",
	       chrom);
	printf("       in genetic map\n");
	exit(5);
      }

      // update beginning / end positions for this chromosome
      chrBegin = geneticMap[chrIdx].second->front().physPos;
      chrEnd = geneticMap[chrIdx].second->back().physPos;
    }

    gotSomeData = true;

    char *posStr = strtok_r(NULL, tab, &saveptr);
    int pos = atoi(posStr);
    if (pos < chrBegin || pos > chrEnd)
      continue; // no genetic map information for this position: skip

    // read/save the ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT fields:
    char *otherFields[7];
    for(int i = 0; i < 7; i++)
      otherFields[i] = strtok_r(NULL, tab, &saveptr);

    // read in/store the haplotypes
    int numRead = 0;
    char *token;
    while(numRead < numFounderHaps && (token = strtok_r(NULL, tab, &saveptr))) {
      char *alleles[2];
      char *saveptr2;
      alleles[0] = strtok_r(token, bar, &saveptr2);
      alleles[1] = strtok_r(NULL, bar, &saveptr2);
      if (alleles[1] == NULL) {
	printf("ERROR: VCF contains data field %s, which is not phased\n",
		token);
	exit(5);
      }
      if (strtok_r(NULL, bar, &saveptr2) != NULL) {
	printf("ERROR: multiple '|' chacters in data field\n");
	exit(5);
      }

      for(int h = 0; h < 2; h++) {
	founderHaps[ numRead + h ] = atoi( alleles[h] );
      }
      numRead += 2;
    }

    if (numRead < numFounderHaps) {
      printf("ERROR: fewer than the needed %d haplotypes found in VCF\n",
	      numFounderHaps);
      exit(6);
    }

    // Print this line to the output file
    fprintf(out, "%s\t%s", chrom, posStr);
    for(int i = 0; i < 7; i++)
      fprintf(out, "\t%s", otherFields[i]);

    for(int fam = 0; fam < numFam; fam++)
      for(int side = 0; side < 2; side++)
	for(int gen = 1; gen < numGen; gen++)
	  if (numSampsToRetain[gen] > 0)
	    for(int ind = 1; ind < numSampsToRetain[gen] + 1; ind++) {
	      for(int h = 0; h < 2; h++) {
		Haplotype &curHap = generations[side][fam][gen][ind].
								haps[h][chrIdx];
		while (curHap.front().endPos < pos) {
		  pop_front(curHap);
		}
		assert(curHap.front().endPos >= pos);
		int foundHapNum = curHap.front().foundHapNum;
		fprintf(out, "%s%d", betweenAlleles[h],
			founderHaps[foundHapNum]);
	      }
	    }
    fprintf(out, "\n");
  }

  free(buffer);
  fclose(out);
  fclose(in);
}

// removes the first element from <vec>. This has time that is linear in the
// number of elements. We could use a list<> instead of a vector<> to avoid this
// linear time, but having random access ability on the vectors is very
// convenient, so we'll live with this.
template<typename T>
void pop_front(vector<T> &vec) {
  vec.erase(vec.begin());
}

void printUsage(char **argv) {
  printf("Usage:\n");
  printf("  %s [numFam] [numGen] [in.dat] [map file] [in.vcf] [out.vcf] [out.bp]\n\n", argv[0]);
  printf("Where:\n");
  printf("  [numFam] is an integer indicating the number of families to simulate\n");
  printf("  [numGen] is an integer specifying the number of generations in the pedigree\n");
  printf("  [map file] contains either a sex averaged genetic map or both male and\n");
  printf("             female maps\n\n");
  printf("  The genetic map file should be formatted with three or four columns:\n\n");
  printf("    [chrom name] [physical position] [map position1 (cM)] <map position2 (cM)>\n\n");
  printf("  [chrom name] must match with the names in the VCF file of phased samples\n");
  printf("  [map position1] gives the sex-averaged map if there are only three columns\n");
  printf("                  or it gives the male map\n");
  printf("  [map position2] gives the female map (if using sex-specific maps)\n");
  exit(1);
}
