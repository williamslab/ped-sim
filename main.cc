#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <random>
#include <chrono>
#include <algorithm>
#include <assert.h>

// TODO: command line options for printing phased VCF and for retaining a
//       specified number of unused samples

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Used to store details about each simulation
struct SimDetails {
  SimDetails(char t, int nFam, int nGen, int *retain, int *mult) {
    type = t;
    numFam = nFam;
    numGen = nGen;
    numSampsToRetain = retain;
    multiplicity = mult;
  }
  // type: either 'f' for full sibs/cousins, 'h' for half sibs/cousins, or
  // 'd' for double cousins
  char type;
  int numFam;
  int numGen;
  int *numSampsToRetain;
  int *multiplicity;
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
// Function decls
void readDat(vector<SimDetails> &simDetails, char *datFile);
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps);
int simulate(vector <SimDetails> &simDetails, Person *****&theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps);
void makeParents(Person *parents[2], char pedType, bool sexSpecificMaps);
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap);
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *outFile);
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, char *inVCFfile, char *outVCFfile,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap);
template<typename T>
void pop_front(vector<T> &vec);
void printUsage(char **argv);

mt19937 randomGen;
uniform_int_distribution<int> coinFlip(0,1);
exponential_distribution<double> crossoverDist(1.0 / 100); // in cM units

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  if (argc != 6) {
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

  vector<SimDetails> simDetails;
  readDat(simDetails, /*datFile=*/ argv[1]);

  vector< pair<char*, vector<PhysGeneticPos>* > > geneticMap;
  bool sexSpecificMaps;
  readMap(geneticMap, /*mapFile=*/ argv[2], sexSpecificMaps);

  // The first index is the pedigree number corresponding to the description of
  // the pedigree to be simulated in the dat file
  // The second index is the family: we replicate the same pedigree structure
  // some number of times as specified in the dat file
  // The third index is the generation number (0-based)
  // The fourth index is the branch of the pedigree; the number of branches
  // is determined by the multiplicity
  // The fifth index is the individual number
  Person *****theSamples;

  printf("Simulating... ");
  fflush(stdout);
  int totalFounderHaps = simulate(simDetails, theSamples, geneticMap,
				  sexSpecificMaps);
  printf("done.\n");

  printf("Printing break points... ");
  fflush(stdout);
  printBPs(simDetails, theSamples, geneticMap, /*outFile=*/ argv[5]);
  printf("done.\n");

  printf("Generating VCF file... ");
  fflush(stdout);
  makeVCF(simDetails, theSamples, totalFounderHaps, /*inVCFfile=*/ argv[3],
	  /*outVCFfile=*/ argv[4], geneticMap);
  printf("done.\n");

  return 0;
}

// Reads in the pedigree formats from the dat file, including the type of the
// pedigree (full, half, or double) and the number of samples to produce in
// every generation
void readDat(vector<SimDetails> &simDetails, char *datFile) {
  // open dat file:
  FILE *in = fopen(datFile, "r");
  if (!in) {
    printf("ERROR: could not open dat file %s!\n", datFile);
    exit(1);
  }

  // dat file contains information about how many samples to generate / store
  // information for in some of the generations; we store this in an array
  // with length equal to the number of generations to be simulated
  int *curNumSampsToRetain = NULL;
  // Have variable number of branches in each generation that is determined by
  // the multiplicity factor stored in the dat file
  int *curMultiplicity = NULL;
  // Have we seen an entry in the dat file for the corresponding generation?
  bool *seen = NULL;
  int curNumGen = 0;
  char curType = '\0';

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *delim = " \t\n";

  int line = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    line++;

    char *token, *saveptr;
    token = strtok_r(buffer, delim, &saveptr);

    if (token == NULL || token[0] == '#') {
      // blank line or comment -- skip
      continue;
    }

    if (strcmp(token, "full") == 0 || strcmp(token, "half") == 0 ||
	strcmp(token, "double") == 0) {
      // new pedigree description
      curType = token[0];
      char *numFamStr = strtok_r(NULL, delim, &saveptr);
      char *numGenStr = strtok_r(NULL, delim, &saveptr);
      if (numFamStr == NULL || numGenStr == NULL ||
				      strtok_r(NULL, delim, &saveptr) != NULL) {
	fprintf(stderr, "ERROR: line %d in dat: expect three fields for pedigree declaration:\n",
		line);
	fprintf(stderr, "       [full/half/double] [numFam] [numGen]\n");
	exit(5);
      }
      int curNumFam = atoi(numFamStr);
      curNumGen = atoi(numGenStr);

      if (curType == 'd' && curNumGen < 3) {
	fprintf(stderr, "ERROR: line %d in dat: request to simulate double cousins with fewer\n",
		line);
	fprintf(stderr, "       than 3 generations\n");
	exit(5);
      }

      curNumSampsToRetain = new int[curNumGen];
      curMultiplicity = new int[curNumGen];
      if (seen)
	delete [] seen;
      seen = new bool[curNumGen];
      for(int i = 0; i < curNumGen; i++) {
	// initially
	curNumSampsToRetain[i] = 0;
	curMultiplicity[i] = 1;
	seen[i] = false;
      }
      curMultiplicity[1] = 2; // default
      simDetails.emplace_back(curType, curNumFam, curNumGen,
			      curNumSampsToRetain, curMultiplicity);
      continue;
    }

    // line contains information about sample storage for the current pedigree
    char *genNumStr = token;
    char *numSampsStr = strtok_r(NULL, delim, &saveptr);
    char *multStr = strtok_r(NULL, delim, &saveptr);

    if (numSampsStr == NULL || strtok_r(NULL, delim, &saveptr) != NULL) {
      printf("ERROR: improper line number %d in dat file: expected two or three fields\n",
	      line);
    }

    int generation = atoi(genNumStr);
    int numSamps = atoi(numSampsStr);

    if (generation < 1 || generation > curNumGen) { // TODO: document
      fprintf(stderr, "ERROR: line %d in dat: generation %d below 1 or above %d (max number\n",
	      line, generation, curNumGen);
      fprintf(stderr, "       of generations)\n");
      exit(1);
    }
    if (numSamps < 0) {
      fprintf(stderr, "ERROR: line %d in dat: in generation %d, number of samples to print\n",
	      line, generation);
      fprintf(stderr, "       below 0\n");
      exit(2);
    }
    if (generation == 1 && numSamps > 1) {
      fprintf(stderr, "ERROR: line %d in dat: in generation 1, if founders are to be printed must\n",
	      line);
      fprintf(stderr, "       list 1 as the number to be printed (others invalid)\n");
      exit(2);
    }

    // subtract 1 because array is 0 based
    if (seen[generation - 1]) {
      fprintf(stderr, "ERROR: line %d in dat: multiple entries for generation %d\n",
	      line, generation);
      exit(2);
    }
    curNumSampsToRetain[generation - 1] = numSamps;
    seen[generation - 1] = true;

    if (multStr != NULL) {
      int mult = atoi(multStr);

      if (generation == curNumGen && mult != 0) {
	fprintf(stderr, "ERROR: line %d in dat: last generation %d multiplicity factor non-zero\n",
	    line, generation);
	exit(2);
      }
      else if (generation != curNumGen && mult <= 0) {
	fprintf(stderr, "ERROR: line %d in dat: in generation %d, multiplicity factor zero or below\n",
		line, generation);
	exit(2);
      }
      else if (generation == 1 && curType != 'f' && mult != 2) {
	fprintf(stderr, "ERROR: line %d in dat: for half and double type pedigrees, generation 1\n",
		line);
	fprintf(stderr, "       multiplicity factor must be 2\n");
	exit(2);
      }
      // TODO: remove this, I think -- any other rules?
      // TODO: one rule we could consider making is that multiplicity of
      // curNumGen - 1 is always 1. If it's > 1, you end up with distinct
      // branches that are in fact full siblings of each other; it's probably
      // simpler if these are forced to be in the same branch
//      else if (generation == 2 && type[0] == 'd' && mult != 1) {
//	fprintf(stderr, "ERROR: line %d in dat: for double type pedigrees, in generation 2,\n",
//		line);
//	fprintf(stderr, "       multiplicity must be 1\n");
//	exit(2);
//      }
      else if (generation < curNumGen)
	// Note: the arrays are 0-based vs. the dat file 1-based in terms of
	// generations. The fact that we assign the index <generation> and
	// not <generation>-1 reflects the fact that the multiplicity in one
	// generation affects the number of branches in the next generation, not
	// the one listed.
	curMultiplicity[generation] = mult;
    }
  }

  for(auto it = simDetails.begin(); it != simDetails.end(); it++) {
    if (it->numSampsToRetain[ it->numGen - 1 ] == 0) {
      // TODO: document
      const char *typeName;
      if (it->type == 'f')
	typeName = "full";
      else if (it->type == 'h')
	typeName = "half";
      else if (it->type == 'd')
	typeName = "double";
      else
	typeName = "error";

      fprintf(stderr, "ERROR: request to simulate '%s' type pedigree, %d families, %d generations\n",
	      typeName, it->numFam, it->numGen);
      fprintf(stderr, "       but no request print any samples from last generation (number %d)\n",
	      it->numGen);
      exit(4);
    }
  }

  if (simDetails.size() == 0) {
    fprintf(stderr, "ERROR: dat file does not contain pedigree descriptions;\n");
    fprintf(stderr, "       nothing to simulate\n");
    exit(3);
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
    else if (mapPos2Str != NULL) {
      fprintf(stderr, "ERROR: expected three columns on all lines in map file but more seen\n");
    }

    curMap->emplace_back(physPos, mapPos1, mapPos2);
  }

  free(buffer);
  fclose(in);
}

// Simulate data for each specified pedigree type for the number of requested
// families. Returns the number of founder haplotypes used to produce these
// simulated samples.
int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     bool sexSpecificMaps) {
  // Note: throughout we use 0-based values for generations though the input
  // dat file is 1-based

  int totalFounderHaps = 0;

  theSamples = new Person****[simDetails.size()];
  for(unsigned int ped = 0; ped < simDetails.size(); ped++) { // for each ped
    char pedType = simDetails[ped].type;
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numSampsToRetain = simDetails[ped].numSampsToRetain;
    // Note: we modify this such that it stores the number of branches in each
    // generation below. To make the code more intuitive, we'll use different
    // names for the <multiplicity> and <numBranches> arrays, but since we
    // never reuse the multiplicity values, they're both in fact the same array
    int *multiplicity = simDetails[ped].multiplicity;
    int *numBranches = multiplicity; // to be updated below

    ////////////////////////////////////////////////////////////////////////////
    // Allocate space and make Person objects for all those we will simulate,
    // assigning sex if <sexSpecificMaps> is true
    theSamples[ped] = new Person***[numFam];
    for (int fam = 0; fam < numFam; fam++) {

      theSamples[ped][fam] = new Person**[numGen];

      // Allocate space for first generation and make top-most generation
      // parents. For full siblings, have only one set of parents, for other
      // types, have two sets:
      if (fam == 0) // only need assign once
	numBranches[0] = (pedType == 'f') ? 1 : 2;
      theSamples[ped][fam][0] = new Person*[ numBranches[0] ];
      // We have an array of size two since, for the half-sib and double cousin
      // simulation, the multiplicity in the first generation is required to
      // be 2. We only use index 0 for full siblings.
      Person *parents[2];
      makeParents(parents, pedType, sexSpecificMaps);

      for(int branch = 0; branch < numBranches[0]; branch++) {
	theSamples[ped][fam][0][branch] = parents[branch];
      }

      // Now all other generations:
      for(int curGen = 1; curGen < numGen; curGen++) {
	// first work out what the number of branches is
	int curMult;
	if (fam == 0) { // only do the below updates once
	  curMult = multiplicity[curGen];
	  if (curGen == 1 && (pedType == 'h' || pedType == 'd')) {
	    // In this case, we should have that curMult == 2, but we don't want
	    // four branches in the second generation but only two. Since the
	    // first generation had numBranches[0] == 2, must change the
	    // multiplication factor here in order to get the desired
	    // numBranches[1] == 2
	    assert(curMult == 2 && numBranches[0] == 2);
	    curMult = 1;
	  }
	  numBranches[curGen] = numBranches[ curGen - 1 ] * curMult;
	}
	else
	  curMult = numBranches[curGen] / numBranches[curGen-1];

	theSamples[ped][fam][curGen] = new Person*[ numBranches[curGen] ];

	// Determine how many samples we need data for in each branch in
	// <curGen>:
	int numPersons = numSampsToRetain[curGen];
	if (numPersons == 0) // not saving, but need parent of next generation
	  numPersons = 1;
	// additional person that is the other parent of next generation
	numPersons++;


	// allocate Persons for each branch of <curGen> and assign their sex
	// ... and do the simulation for these allocated individuals
	for(int branch = 0; branch < numBranches[curGen]; branch++) {
	  theSamples[ped][fam][curGen][branch] = new Person[numPersons];

	  if (sexSpecificMaps) {
	    if (pedType == 'd' && curGen == 1) {
	      // won't randomize the sex assignment for the second generation
	      // when simulating double cousins. Want to use the convention
	      // that ind 1 from one branch has children with ind 0 from the
	      // other branch. Therefore, we don't want the sexes to be
	      // randomly assigned
	      theSamples[ped][fam][curGen][branch][1].sex = 1;
	    }
	    else {
	      // the two individuals that reproduce are index 0 (a founder) and
	      // 1 randomly decide which one to make female
	      int femaleIdx = coinFlip(randomGen);
	      theSamples[ped][fam][curGen][branch][femaleIdx].sex = 1;
	    }
	  }

	  /////////////////////////////////////////////////////////////////////
	  // All samples allocated for this pedigree/family/generation/branch:
	  // simulate the actual samples

	  // for each chromosome:
	  for(unsigned int chrIdx = 0; chrIdx < geneticMap.size(); chrIdx++) {
	    vector<PhysGeneticPos> *curMap = geneticMap[chrIdx].second;

	    // Make trivial haplotypes for generation 0 founders:
	    // no recombinations in founders
	    Segment trivialSeg;
	    trivialSeg.endPos = curMap->back().physPos;
	    if (curGen == 1) {
	      for(int par = 0; par < 2; par++) { // each parent
		if (pedType == 'f' && branch > 0)
		  // already initialized parents for branch == 0, and for
		  // pedType == 'f' they're the same objects so no need to
		  // change anything
		  break;
		for(int h = 0; h < 2; h++) {
		  if (pedType == 'h' && branch == 1 && par == 0)
		    // shared parent -- make identical to branch 0
		    trivialSeg.foundHapNum =
		      theSamples[ped][fam][0][0][0].haps[h].back().
							    back().foundHapNum;
		  else {
		    if (chrIdx == 0)
		      // new sample: new haplotype index:
		      trivialSeg.foundHapNum = totalFounderHaps++;
		    else
		      // want the same founder on all chromosomes, so access
		      // the haplotype number assigned to the previous
		      // chromosome for this person:
		      trivialSeg.foundHapNum =
			theSamples[ped][fam][0][branch][par].haps[h].back().
							    back().foundHapNum;
		  }
		  // Note: pedType == 'd' will create two separate sets of
		  // parents for the two branches/sides automatically with this
		  // code. This is what we want.

		  // makes a copy of <trivialSeg>, can reuse
		  theSamples[ped][fam][0][branch][par].haps[h].emplace_back();
		  theSamples[ped][fam][0][branch][par].haps[h].back().push_back(
								    trivialSeg);
		}
	      }
	    }

	    // First make trivial haplotypes for the two founders in <curGen>;
	    // use convention that sample 0 is the founder in each branch:
	    //
	    // no founders in the last generation and
	    // TODO: document this
	    // no founders in second generation when pedType == 'd' -- the two
	    // full sibs on both sides reproduce to create the next generation
	    if (curGen != numGen - 1 && !(pedType == 'd' && curGen == 1)) {
	      for(int h = 0; h < 2; h++) {
		// 4 founder haplotypes per generation
		if (chrIdx == 0)
		  trivialSeg.foundHapNum = totalFounderHaps++;
		else
		  // want the same founder on all chromosomes, so access
		  // the haplotype number assigned to the previous
		  // chromosome for this person:
		  trivialSeg.foundHapNum =
			theSamples[ped][fam][curGen][branch][0].haps[h].back().
							    back().foundHapNum;

		// the following copies <trivialSeg>, so we can reuse it
		theSamples[ped][fam][curGen][branch][0].haps[h].emplace_back();
		theSamples[ped][fam][curGen][branch][0].haps[h].back().
							  push_back(trivialSeg);
	      }
	    }

	    // Now simulate the non-founders in <branch> of <curGen>
	    // We take parents from the current branch number divided by the
	    // multiplicity factor. Thus, each branch in the previous generation
	    // produces <curMult> new branches in the current one:
	    int prvBrch = branch / curMult;
	    int startInd = 1;
	    if (pedType == 'd' && curGen == 1)
	      // As per above, individual 0 is not a founder in the second
	      // generation for double cousins simulations, so want to simulate
	      // this person (as a sibling of individual 1)
	      startInd = 0;
	    for(int ind = startInd; ind < numPersons; ind++) {
	      // If we're using sex-specific maps, the two parents' sexes should
	      // differ:
	      assert(!sexSpecificMaps ||
		      theSamples[ped][fam][curGen-1][prvBrch][0].sex !=
				theSamples[ped][fam][curGen-1][prvBrch][1].sex);
	      for(int parIdx = 0; parIdx < 2; parIdx++) {
		// haplotype index for the simulated sample
		int hapIdx = parIdx;
		if (sexSpecificMaps)
		  // match the sex of the parent if using sex-specific maps
		  hapIdx = theSamples[ped][fam][curGen-1][prvBrch][parIdx].sex;

		theSamples[ped][fam][curGen][branch][ind].haps[hapIdx].
								 emplace_back();
		Person *theParent;
		if (pedType == 'd' && curGen == 2 && parIdx == 0) {
		  // for double cousins, in the third generation, have ind 1
		  // on one side/branch have children with ind 0 on the other
		  // side/branch
		  assert(curMult == 1 && numBranches[1] == 2); // sanity check
		  int othrBrnch = branch ^ 1;
		  theParent =&theSamples[ped][fam][curGen-1][othrBrnch][parIdx];
		}
		else {
		  theParent = &theSamples[ped][fam][curGen-1][prvBrch][parIdx];
		}

		generateHaplotype(
		  theSamples[ped][fam][curGen][branch][ind].haps[hapIdx].back(),
		  *theParent, curMap);
	      } // <parIdx> (simulate each transmitted haplotype)
	    } // <ind>
	  } // <geneticMap> (chroms)
	} // <branch>
      } // <curGen>
    } // <fam>

  } // <ped>

  return totalFounderHaps;
}

// Allocate Person objects for the parents in the top-most generation. The
// way this is setup depends on the pedigree type: full, half, or double
void makeParents(Person *parents[2], char pedType, bool sexSpecificMaps) {
  if (pedType == 'f') {
    // full siblings in second generation: same parents for all branches
    parents[0] = new Person[2];
    if (sexSpecificMaps)
      parents[0][1].sex = 1;
    parents[1] = NULL;
  }
  else if (pedType == 'h') {
    // half-siblings in second generation: one shared parent because we're not
    // using pointers here but actual Person objects, we'll make two copies of
    // the same founder; this is done in simulate()
    parents[0] = new Person[2];
    parents[1] = new Person[2];

    // decide whether the shared parent -- index 0 in each parents array --
    // is male or female
    int sharedParSex = coinFlip(randomGen);
    parents[0][0].sex = parents[1][0].sex = sharedParSex;
    // opposite for other parent:
    parents[0][1].sex = parents[1][1].sex = 1 ^ sharedParSex;
  }
  else if (pedType == 'd') {
    // double cousins: two completely separate sets of parents for the two
    // sides/branches, and two full siblings from each side/branch reproduce to
    // create the double cousins in the third generation
    parents[0] = new Person[2];
    parents[1] = new Person[2];
    if (sexSpecificMaps)
      parents[0][1].sex = parents[1][1].sex = 1;
  }
  else {
    fprintf(stderr, "ERROR: unsupported pedigree type %c\n", pedType);
    exit(5);
  }
}

// Simulate one haplotype <toGenerate> by sampling crossovers and switching
// between the two haplotypes stored in <parent>. Uses the genetic map stored
// in <curMap> which is expected to correspond to the chromosome being
// simulated and may contain either one sex-averaged map or a male and female
// map.
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

// Print the break points to <outFile>
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *outFile) {
  FILE *out = fopen(outFile, "w");
  if (!out) {
    printf("ERROR: could not open output file %s!\n", outFile);
    exit(1);
  }

  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    char pedType = simDetails[ped].type;
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numSampsToRetain = simDetails[ped].numSampsToRetain;
    // The simulate() method updated the multiplicity array to store the
    // number of branches per generation
    int *numBranches = simDetails[ped].multiplicity;

    for(int fam = 0; fam < numFam; fam++) {
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  if (numSampsToRetain[gen] > 0) {
	    for(int ind = 0; ind < numSampsToRetain[gen] + 1; ind++) {
	      if (gen == 0) {
		// should only be one set of parents / branch for pedType == 'f'
		assert(pedType != 'f' || branch == 0);
		if (pedType == 'h' && branch > 0 && ind == 0)
		  // for pedType == 'h', ind 0 in the top-most generation is
		  // the same person; need only print once, so skip
		  continue;
	      }
	      if (ind == 0 && gen == numGen-1)
		// no founders (by convention stored as ind == 0) in the
		// last generation
		continue;

	      for(int h = 0; h < 2; h++) {
		int sex = theSamples[ped][fam][gen][branch][ind].sex;
		fprintf(out, "%c%d_f%d_b%d_g%d_i%d s%d h%d", pedType, ped+1,
			fam+1, branch, gen+1, ind, sex, h);

		for(unsigned int chr = 0; chr < geneticMap.size(); chr++) {
		  // print chrom name and starting position
		  fprintf(out, " %s|%d", geneticMap[chr].first,
			  geneticMap[chr].second->front().physPos);
		  Haplotype &curHap = theSamples[ped][fam][gen][branch][ind].
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
  }

  fclose(out);
}

// Given the simulated break points for individuals in each pedigree/family
// stored in <theSamples> and other necessary information, reads input VCF
// format data from the file named <inVCFfile> and prints the simulated
// haplotypes for each sample to <outVCFfile> in VCF format.
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, char *inVCFfile, char *outVCFfile,
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
  // technically tab and newline; we want the latter so that the last sample id
  // on the header line doesn't include the newline character in it
  const char *tab = "\t\n";
  const char *bar = "|";
  // Below when we print the VCF output, we alternate printing tab
  // and / between successive alleles. Make this simpler with:
  const char betweenAlleles[2] = { '\t', '/' };

  char **hapAlleles = NULL; // stores all alleles from input sample
  char **founderHaps = new char*[totalFounderHaps]; // alleles for founder haps

  // iterate over chromosomes in the genetic map
  unsigned int chrIdx = 0; // index of current chromosome number;
  char *chrName = geneticMap[chrIdx].first;
  int chrBegin = geneticMap[chrIdx].second->front().physPos;
  int chrEnd = geneticMap[chrIdx].second->back().physPos;

  bool gotSomeData = false;

  int numInputSamples = 0;
  vector<int> shuffHaps; // For randomizing the assigned haplotypes

  while (getline(&buffer, &bytesRead, in) >= 0) { // read each line of input VCF
    if (buffer[0] == '#' && buffer[1] == '#') {
      // header line: print to output
      fprintf(out, "%s", buffer);
      continue;
    }

    if (buffer[0] == '#') {
      // header line with sample ids

      // skip all the header fields relating to meta-data:
      char *saveptr;
      char *token = strtok_r(buffer, tab, &saveptr);
      for(int i = 1; i < 9; i++)
	token = strtok_r(NULL, tab, &saveptr);

      // now parse / store the sample ids:
      vector<char*> sampleIds;
      while ((token = strtok_r(NULL, tab, &saveptr))) {
	sampleIds.push_back(token);
      }
      numInputSamples = sampleIds.size();
      hapAlleles = new char*[numInputSamples * 2]; // 2 for diploid samples

      // Next generate an ordered list of haplotype indexes (2 * sample_index)
      // and randomly shuffle it. The index of shuffHaps is the sample index
      // and the value stored at the index is the (randomly assigned) haplotype
      // index for its first allele. Ultimately we only care about the haplotype
      // indexes that are < totalFounderHaps; those >= totalFounderHaps will
      // not be used for the simulated samples and we will print their data
      // (see below)
      for(int i = 0; i < numInputSamples; i++)
	shuffHaps.push_back(2 * i);
      shuffle(shuffHaps.begin(), shuffHaps.end(), randomGen);

      if (2 * numInputSamples < totalFounderHaps) {
	fprintf(stderr, "ERROR: need %d founder haplotypes, but input only contains %d\n",
		totalFounderHaps, 2 * numInputSamples);
	exit(5);
      }

      // Now print the header line indicating fields and sample ids for the
      // output VCF
      fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

      // print sample ids:
      for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
	char pedType = simDetails[ped].type;
	int numFam = simDetails[ped].numFam;
	int numGen = simDetails[ped].numGen;
	int *numSampsToRetain = simDetails[ped].numSampsToRetain;
	// The simulate() method updated the multiplicity array to store the
	// number of branches per generation
	int *numBranches = simDetails[ped].multiplicity;

	for(int fam = 0; fam < numFam; fam++)
	  for(int gen = 0; gen < numGen; gen++)
	    for(int branch = 0; branch < numBranches[gen]; branch++)
	      if (numSampsToRetain[gen] > 0)
		for(int ind = 0; ind < numSampsToRetain[gen] + 1; ind++) {
		  if (gen == 0) {
		    // should only be one set of parents / branch for
		    // pedType == 'f'
		    assert(pedType != 'f' || branch == 0);
		    if (pedType == 'h' && branch > 0 && ind == 0)
		      // for pedType == 'h', ind 0 in the top-most generation is
		      // the same person; need only print once, so skip
		      continue;
		  }
		  if (ind == 0 && gen == numGen-1)
		    // no founders (by convention stored as ind == 0) in the
		    // last generation
		    continue;
		  fprintf(out, "\t%c%d_f%d_b%d_g%d_i%d", pedType, ped+1, fam+1,
			  branch, gen+1, ind);
		}
      }

      // print ids for the samples that will not be used in generating the
      // simulated data
      for(int i = 0; i < numInputSamples; i++) {
	if (shuffHaps[i] >= totalFounderHaps)
	  fprintf(out, "\t%s", sampleIds[i]);
      }

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
    int inputIndex = 0;
    int numStored = 0;
    char *token;
    while((token = strtok_r(NULL, tab, &saveptr)) &&
					      numStored < numInputSamples * 2) {
      char *alleles[2];
      char *saveptr2;
      alleles[0] = strtok_r(token, bar, &saveptr2);
      alleles[1] = strtok_r(NULL, bar, &saveptr2);

      for(int h = 0; h < 2; h++)
	hapAlleles[numStored++] = alleles[h];

      int founderIndex = shuffHaps[ inputIndex ];
      inputIndex++;
      if (founderIndex < totalFounderHaps) {
	for(int h = 0; h < 2; h++)
	  founderHaps[founderIndex + h] = alleles[h];
      }

      // error check:
      if (alleles[1] == NULL) {
	printf("ERROR: VCF contains data field %s, which is not phased\n",
		token);
	exit(5);
      }
      if (strtok_r(NULL, bar, &saveptr2) != NULL) {
	printf("ERROR: multiple '|' chacters in data field\n");
	exit(5);
      }
    }

    bool fewer = numStored < numInputSamples * 2;
    bool more = token != NULL;
    if (fewer || more) {
      printf("ERROR: line in VCF file has data for %s than the indicated %d samples\n",
	     (more) ? "more" : "fewer", numInputSamples);
      exit(6);
    }

    // Print this line to the output file
    fprintf(out, "%s\t%s", chrom, posStr);
    for(int i = 0; i < 7; i++)
      fprintf(out, "\t%s", otherFields[i]);

    for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
      char pedType = simDetails[ped].type;
      int numFam = simDetails[ped].numFam;
      int numGen = simDetails[ped].numGen;
      int *numSampsToRetain = simDetails[ped].numSampsToRetain;
      // The simulate() method updated the multiplicity array to store the
      // number of branches per generation
      int *numBranches = simDetails[ped].multiplicity;

      for(int fam = 0; fam < numFam; fam++)
	for(int gen = 0; gen < numGen; gen++)
	  for(int branch = 0; branch < numBranches[gen]; branch++)
	    if (numSampsToRetain[gen] > 0)
	      for(int ind = 0; ind < numSampsToRetain[gen] + 1; ind++) {
		if (gen == 0) {
		  // should only be one set of parents / branch for
		  // pedType == 'f'
		  assert(pedType != 'f' || branch == 0);
		  if (pedType == 'h' && branch > 0 && ind == 0)
		    // for pedType == 'h', ind 0 in the top-most generation is
		    // the same person; need only print once, so skip
		    continue;
		}
		if (ind == 0 && gen == numGen-1)
		  // no founders (by convention stored as ind == 0) in the
		  // last generation
		  continue;
		for(int h = 0; h < 2; h++) {
		  Haplotype &curHap = theSamples[ped][fam][gen][branch][ind].
								haps[h][chrIdx];
		  while (curHap.front().endPos < pos) {
		    pop_front(curHap);
		  }
		  assert(curHap.front().endPos >= pos);
		  int foundHapNum = curHap.front().foundHapNum;
		  fprintf(out, "%c%s", betweenAlleles[h],
			  founderHaps[foundHapNum]);
		}
	      }
    }
    // print data for samples not used in generating the simulated data
    for(int i = 0; i < numInputSamples; i++) {
      if (shuffHaps[i] >= totalFounderHaps)
	for(int h = 0; h < 2; h++)
	  fprintf(out, "%c%s", betweenAlleles[h],
		  hapAlleles[ shuffHaps[i] + h ]);
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
  printf("  %s [in.dat] [map file] [in.vcf] [out.vcf] [out.bp]\n\n", argv[0]);
  printf("Where:\n");
//  printf("  [numFam] is an integer indicating the number of families to simulate\n");
//  printf("  [numGen] is an integer specifying the number of generations in the pedigree\n");
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
