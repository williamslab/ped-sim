// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <random>
#include <sys/time.h>
#include <algorithm>
#include <assert.h>
#include <boost/dynamic_bitset.hpp>
#include <zlib.h>
#include "cmdlineopts.h"
#include "cointerfere.h"

// TODO! only use sexConstraints array when there are sex-specific maps?
// TODO! make branchNumSpouses positive

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Used to store details about each simulation
struct SimDetails {
  SimDetails(int nFam, int nGen, int *retain, int *branches, int **parents,
	     int **sexes, int **spouses, char *theName) {
    numFam = nFam;
    numGen = nGen;
    numSampsToRetain = retain;
    numBranches = branches;
    branchParents = parents;
    sexConstraints = sexes;
    branchNumSpouses = spouses;
    name = new char[ strlen(theName) + 1 ];
    strcpy(name, theName);
  }
  int numFam;
  int numGen;
  int *numSampsToRetain;
  int *numBranches;
  int **branchParents;
  int **sexConstraints;
  int **branchNumSpouses;
  char *name;
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
void readDef(vector<SimDetails> &simDetails, char *defFile);
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				int *&thisGenBranchParents,
				int *prevGenSpouseNum = NULL,
				vector<bool> *branchParentsAssigned = NULL);
void readBranchParents(int prevGenNumBranches, int thisGenNumBranches,
		       int *&thisGenBranchParents, int *&prevGenSexConstraints,
		       int *&prevGenSpouseNum,
		       vector<bool> &branchParentsAssigned,
		       vector< boost::dynamic_bitset<>* > &spouseDependencies,
		       const char *delim, char *&saveptr, char *&endptr,
		       int line);
void updateSexConstraints(int *&prevGenSexConstraints, int parIdx[2],
			  int prevGenNumBranches,
			  vector< boost::dynamic_bitset<>*> &spouseDependencies,
			  int line);
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps);
void readInterfere(vector<COInterfere> &coIntf, char *interfereFile,
		   vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		   bool &sexSpecificMaps);
int simulate(vector <SimDetails> &simDetails, Person *****&theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps, vector<COInterfere> &coIntf);
void getPersonCounts(int curGen, int numGen, int branch, int *numSampsToRetain,
		     int **branchParents, int **branchNumSpouses,
		     int &numFounders, int &numNonFounders);
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap,
		       vector<COInterfere> &coIntf, unsigned int chrIdx);
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *bpFile);
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, char *inVCFfile, char *outFileBuf,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     FILE *outs[2]);
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
	      char *famFile);
template<typename T>
void pop_front(vector<T> &vec);

mt19937 randomGen;
uniform_int_distribution<int> coinFlip(0,1);
exponential_distribution<double> crossoverDist(1.0);

struct fileOrGZ {
  FILE *fp;
  gzFile gfp;
  bool isGZ;
};

bool fileOrGZ_open(fileOrGZ &fgz, const char *filename, const char *mode,
		   bool isGZ);
int fileOrGZ_getline(char **lineptr, size_t *n, fileOrGZ &fgz);
int fileOrGZ_printf(fileOrGZ &fgz, const char *format, ...);
int fileOrGZ_close(fileOrGZ &fgz);


////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  int outPrefixLen = strlen(CmdLineOpts::outPrefix);
  char *outFile = new char[ outPrefixLen + 7 + 1 ]; //+7 for .vcf.gz, + 1 for \0

  // open the log file
  sprintf(outFile, "%s.log", CmdLineOpts::outPrefix);
  FILE *log = fopen(outFile, "w");
  if (!log) {
    printf("ERROR: could not open log file %s!\n", outFile);
    perror("open");
    exit(1);
  }

  // seed random number generator if needed
  if (CmdLineOpts::autoSeed) {
    CmdLineOpts::randSeed = random_device().entropy();
    if (CmdLineOpts::randSeed == 0 && random_device().entropy() == 0) {
      // random_device is not a real random number generator: fall back on
      // using time to generate a seed:
      timeval tv;
      gettimeofday(&tv, NULL);
      CmdLineOpts::randSeed = tv.tv_sec * tv.tv_usec;
    }
  }
  randomGen.seed(CmdLineOpts::randSeed);

  FILE *outs[2] = { stdout, log };

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Pedigree simulator!  v%s    (Released %s)\n\n",
	    VERSION_NUMBER, RELEASE_DATE);

    fprintf(outs[o], "  Def file:\t\t%s\n", CmdLineOpts::defFile);
    fprintf(outs[o], "  Map file:\t\t%s\n", CmdLineOpts::mapFile);
    fprintf(outs[o], "  Input VCF:\t\t%s\n", CmdLineOpts::inVCFfile);
    fprintf(outs[o], "  Output prefix:\t%s\n\n", CmdLineOpts::outPrefix);

    fprintf(outs[o], "  Random seed:\t\t%u\n\n", CmdLineOpts::randSeed);

    fprintf(outs[o], "  Interference file:\t%s\n\n",
	    CmdLineOpts::interfereFile == NULL ? "[none]" :
						    CmdLineOpts::interfereFile);

    fprintf(outs[o], "  Genotype error rate:\t%.1le\n",
	    CmdLineOpts::genoErrRate);
    fprintf(outs[o], "  Opposite homozygous error rate:\t%.2lf\n",
	    CmdLineOpts::homErrRate);
    fprintf(outs[o], "  Missingness rate:\t%.1le\n\n",
	    CmdLineOpts::missRate);

    if (CmdLineOpts::retainExtra < 0) {
      fprintf(outs[o], "  Retaining all unused samples in output VCF\n");
    }
    else if (CmdLineOpts::retainExtra == 0) {
      fprintf(outs[o], "  Not retaining extra samples in output VCF (printing only simulated samples)\n");
    }
    else {
      fprintf(outs[o], "  Retaining %d unused samples in output VCF\n",
	      CmdLineOpts::retainExtra);
    }
    if (CmdLineOpts::keepPhase) {
      fprintf(outs[o], "  Output VCF will contain phased data\n\n");
    }
    else {
      fprintf(outs[o], "  Output VCF will contain unphased data\n\n");
    }
  }

  vector<SimDetails> simDetails;
  readDef(simDetails, CmdLineOpts::defFile);

  vector< pair<char*, vector<PhysGeneticPos>* > > geneticMap;
  bool sexSpecificMaps;
  readMap(geneticMap, CmdLineOpts::mapFile, sexSpecificMaps);

  vector<COInterfere> coIntf;
  if (CmdLineOpts::interfereFile) {
    readInterfere(coIntf, CmdLineOpts::interfereFile, geneticMap,
		  sexSpecificMaps);
  }

  // The first index is the pedigree number corresponding to the description of
  // the pedigree to be simulated in the def file
  // The second index is the family: we replicate the same pedigree structure
  // some number of times as specified in the def file
  // The third index is the generation number (0-based)
  // The fourth index is the branch of the pedigree
  // The fifth index is the individual number
  Person *****theSamples;

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Simulating haplotype transmissions... ");
    fflush(outs[o]);
  }
  int totalFounderHaps = simulate(simDetails, theSamples, geneticMap,
				  sexSpecificMaps, coIntf);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Printing break points... ");
    fflush(outs[o]);
  }
  sprintf(outFile, "%s.bp", CmdLineOpts::outPrefix);
  printBPs(simDetails, theSamples, geneticMap, /*bpFile=*/ outFile);
  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "done.\n");
  }

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Generating output VCF... ");
    fflush(outs[o]);
  }
  makeVCF(simDetails, theSamples, totalFounderHaps, CmdLineOpts::inVCFfile,
	  /*outVCFfile=*/ outFile, geneticMap, outs);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "Printing fam file... ");
  fflush(stdout);
  sprintf(outFile, "%s.fam", CmdLineOpts::outPrefix);
  printFam(simDetails, theSamples, /*famFile=*/ outFile);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  fclose(log);

  return 0;
}

// Reads in the pedigree formats from the def file, including the type of the
// pedigree (full, half, or double) and the number of samples to produce in
// every generation
void readDef(vector<SimDetails> &simDetails, char *defFile) {
  // open def file:
  FILE *in = fopen(defFile, "r");
  if (!in) {
    printf("ERROR: could not open def file %s!\n", defFile);
    perror("open");
    exit(1);
  }

  // def file contains information about how many samples to generate / store
  // information for in some of the generations; we store this in an array
  // with length equal to the number of generations to be simulated
  int *curNumSampsToRetain = NULL;
  // Have variable number of branches in each generation
  int *curNumBranches = NULL;
  // Who are the parents of each branch in each generation?
  // Contains <curNumGen> rows, and <2*curNumBranches[gen]> columns on each row.
  // Stores the previous branch numbers that contain the two parents.
  // Negative values correspond to founders that are stored in the same
  // branch number as the other parent.
  int **curBranchParents = NULL;
  // Gives numerical values indicating dependencies of sex assignments for each
  // branch. For example, if the person in branch 1 has children with the
  // individual in branch 2 and 3, branch 2 and 3 must have the same sex.
  int **curSexConstraints = NULL;
  // Counts of number of non-founder spouses for each generation/branch
  int **curBranchNumSpouses = NULL;
  int curNumGen = 0;
  // for ensuring generations are in increasing order. This requirement arises
  // from the fact that we assign the number of branches in each generation to
  // be equal to the previous generation (except generation 2), so we need to
  // know which generation we've assigned the generation numbers to and update
  // branch counts for any generations that aren't explicitly listed.
  int lastReadGen = -1;

  // Tracks whether there has been an explicit assignment of the parents of
  // each branch to avoid double assignments and giving default assignments.
  vector<bool> branchParentsAssigned;
  // Stores sets of individuals that are required to have the same and/or
  // opposite sex assignments by virtue of their being spouses.
  vector< boost::dynamic_bitset<>* > spouseDependencies;

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *delim = " \t\n";

  int line = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    line++;

    char *token, *saveptr, *endptr;
    token = strtok_r(buffer, delim, &saveptr);

    if (token == NULL || token[0] == '#') {
      // blank line or comment -- skip
      continue;
    }

    if (strcmp(token, "def") == 0) {

      /////////////////////////////////////////////////////////////////////////
      // parse new pedigree description
      char *name = strtok_r(NULL, delim, &saveptr);
      char *numFamStr = strtok_r(NULL, delim, &saveptr);
      char *numGenStr = strtok_r(NULL, delim, &saveptr);
      if (name == NULL || numFamStr == NULL || numGenStr == NULL ||
				      strtok_r(NULL, delim, &saveptr) != NULL) {
	fprintf(stderr, "ERROR: line %d in def: expect four fields for pedigree definition:\n",
		line);
	fprintf(stderr, "       def [name] [numFam] [numGen]\n");
	exit(5);
      }
      int curNumFam = strtol(numFamStr, &endptr, 10);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: line %d in def: expected number of families to simulate as second token\n",
		line);
	if (errno != 0)
	  perror("strtol");
	exit(2);
      }
      curNumGen = strtol(numGenStr, &endptr, 10);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: line %d in def: expected number of generations to simulate as third",
		line);
	fprintf(stderr, "      token\n");
	if (errno != 0)
	  perror("strtol");
	exit(2);
      }

      // TODO: slow linear search to ensure lack of repetition of the pedigree
      // names; probably fast enough
      for(auto it = simDetails.begin(); it != simDetails.end(); it++) {
	if (strcmp(it->name, name) == 0) {
	  fprintf(stderr, "ERROR: line %d in def: name of pedigree is same as previous pedigree\n",
		  line);
	  exit(5);
	}
      }

      curNumSampsToRetain = new int[curNumGen];
      curNumBranches = new int[curNumGen];
      curBranchParents = new int*[curNumGen];
      curSexConstraints = new int*[curNumGen];
      curBranchNumSpouses = new int*[curNumGen];
      if (lastReadGen >= 0)
	lastReadGen = -1; // reset

      for(int gen = 0; gen < curNumGen; gen++) {
	// initially
	curNumSampsToRetain[gen] = 0;
	// set to -1 initially so we know these are unassigned; will update
	// later
	curNumBranches[gen] = -1;
	curBranchParents[gen] = NULL;
	curSexConstraints[gen] = NULL;
	curBranchNumSpouses[gen] = NULL;
      }
      simDetails.emplace_back(curNumFam, curNumGen, curNumSampsToRetain,
			      curNumBranches, curBranchParents,
			      curSexConstraints, curBranchNumSpouses, name);
      continue;
    }

    ///////////////////////////////////////////////////////////////////////////
    // parse line with information about a generation in the current pedigree

    // is there a current pedigree?
    if (curNumSampsToRetain == NULL) {
      fprintf(stderr, "ERROR: line %d in def: expect four fields for pedigree definition:\n",
	      line);
      fprintf(stderr, "       def [name] [numFam] [numGen]\n");
      exit(5);
    }

    char *genNumStr = token;
    char *numSampsStr = strtok_r(NULL, delim, &saveptr);
    char *branchStr = strtok_r(NULL, delim, &saveptr);

    int generation = strtol(genNumStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: line %d in def: expected generation number or \"def\" as first token\n",
	  line);
      if (errno != 0)
	perror("strtol");
      exit(2);
    }

    if (numSampsStr == NULL) {
      printf("ERROR: improper line number %d in def file: expected at least two fields\n",
	      line);
    }
    int numSamps = strtol(numSampsStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: line %d in def: expected number of samples to print as second token\n",
	  line);
      if (errno != 0)
	perror("strtol");
      exit(2);
    }

    if (generation < 1 || generation > curNumGen) {
      fprintf(stderr, "ERROR: line %d in def: generation %d below 1 or above %d (max number\n",
	      line, generation, curNumGen);
      fprintf(stderr, "       of generations)\n");
      exit(1);
    }
    if (numSamps < 0) {
      fprintf(stderr, "ERROR: line %d in def: in generation %d, number of samples to print\n",
	      line, generation);
      fprintf(stderr, "       below 0\n");
      exit(2);
    }
    if (generation == 1 && numSamps > 1) {
      fprintf(stderr, "ERROR: line %d in def: in generation 1, if founders are to be printed must\n",
	      line);
      fprintf(stderr, "       list 1 as the number to be printed (others invalid)\n");
      exit(2);
    }

    if (generation <= lastReadGen) {
      fprintf(stderr, "ERROR: line %d in def: generation numbers must be in increasing order\n",
	      line);
      exit(7);
    }

    // if <curNumBranches> != -1, have prior definition for generation.
    // subtract 1 from <generation> because array is 0 based
    if (curNumBranches[generation - 1] != -1) {
      fprintf(stderr, "ERROR: line %d in def: multiple entries for generation %d\n",
	      line, generation);
      exit(2);
    }
    curNumSampsToRetain[generation - 1] = numSamps;

    // Assign number of branches (and parents) for generations that are not
    // explicitly listed. In general the number of branches is equal to the
    // number in the previous generation. The exceptions are generation 1 which
    // defaults to 1 branch, and generation 2 which defaults to 2 branches when
    // generation 1 has only 1 branch (otherwise it's assigned the same as the
    // previous generation)
    for(int i = lastReadGen + 1; i < generation - 1; i++) {
      if (i == 0)
	curNumBranches[0] = 1;
      else if (i == 1 && curNumBranches[0] == 1)
	curNumBranches[1] = 2;
      else
	curNumBranches[i] = curNumBranches[i-1];

      // assign default parents for each branch:
      if (i > 0)
	assignDefaultBranchParents(curNumBranches[i-1], curNumBranches[i],
				   curBranchParents[i]);
    }

    int thisGenNumBranches;
    if (branchStr != NULL) {
      thisGenNumBranches = strtol(branchStr, &endptr, 10);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: line %d in def: optional third token must be numerical value giving\n",
		line);
	fprintf(stderr, "      number of branches\n");
	if (errno != 0)
	  perror("strtol");
	exit(2);
      }

      if (thisGenNumBranches <= 0) {
	fprintf(stderr, "ERROR: line %d in def: in generation %d, branch number zero or below\n",
		line, generation);
	exit(2);
      }
      else {
	curNumBranches[generation - 1] = thisGenNumBranches;
      }
    }
    else {
      if (generation - 1 == 0)
	thisGenNumBranches = 1;
      else if (generation - 1 == 1 && curNumBranches[0] == 1)
	thisGenNumBranches = 2;
      else
	thisGenNumBranches = curNumBranches[generation-2];

      curNumBranches[generation - 1] = thisGenNumBranches;
    }

    lastReadGen = generation - 1;

    // now read in and assign (if only using the defaults) the branch parents
    // for this generation. Note that in the first generation, all individuals
    // are necessarily founders so there should not be any specification.
    if (generation - 1 > 0)
      readBranchParents(/*prevGenBranches=*/curNumBranches[generation - 2],
			thisGenNumBranches, curBranchParents[generation - 1],
			/*prevGenSexConst=*/curSexConstraints[generation - 2],
			/*prevSpouseNum=*/curBranchNumSpouses[generation - 2],
			branchParentsAssigned, spouseDependencies, delim,
			saveptr, endptr, line);
    else if (strtok_r(NULL, delim, &saveptr) != NULL) {
      fprintf(stderr, "ERROR: line %d in def: first generation cannot have parent specifications\n",
	      line);
      exit(8);
    }
  }


  for(auto it = simDetails.begin(); it != simDetails.end(); it++) {
    if (it->numSampsToRetain[ it->numGen - 1 ] == 0) {
      fprintf(stderr, "ERROR: request to simulate pedigree \"%s\" with %d generations\n",
	      it->name, it->numGen);
      fprintf(stderr, "       but no request to print any samples from last generation (number %d)\n",
	      it->numGen);
      exit(4);
    }
  }

  if (simDetails.size() == 0) {
    fprintf(stderr, "ERROR: def file does not contain pedigree definitions;\n");
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
    perror("open");
    exit(1);
  }

  char *curChr = NULL, *endptr;
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

    // need a new entry in <geneticMap> for a new chrom?
    if (curChr == NULL || strcmp(chrom, curChr) != 0) {
      curChr = (char *) malloc(strlen(chrom) + 1);
      strcpy(curChr, chrom);
      curMap = new vector<PhysGeneticPos>;
      geneticMap.emplace_back(curChr, curMap);
    }

    int physPos;
    double mapPos1, mapPos2 = 0.0;
    physPos = strtol(physPosStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: could not parse column 2 of map file as integer\n");
      if (errno != 0)
	perror("strtol");
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
    }

    curMap->emplace_back(physPos, mapPos1, mapPos2);
  }

  free(buffer);
  fclose(in);
}

// Reads in interference parameters nu and p for males and females from
// <interfereFile> and stores them in <intfParams>.
void readInterfere(vector<COInterfere> &coIntf, char *interfereFile,
		   vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		   bool &sexSpecificMaps) {
  if (!sexSpecificMaps) {
    fprintf(stderr, "ERROR: Must use sex specific genetic maps in order to simulate with interference\n");
    exit(6);
  }

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  const char *delim = " \t\n";

  FILE *in = fopen(interfereFile, "r");
  if (!in) {
    printf("ERROR: could not open interference file %s!\n", interfereFile);
    exit(1);
  }

  // Which chromosome index (into <geneticMap>) are we on? This allows us to
  // ensure the names of the chromosomes listed in the interference file match
  // those in <geneticMap>
  unsigned int chrIdx = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    char *chrom, *nuStr[2], *pStr[2];
    char *saveptr, *endptr;
    double nu[2], p[2];

    if (buffer[0] == '#')
      continue; // comment

    // get chromosome tokens:
    chrom = strtok_r(buffer, delim, &saveptr);

    if (chrIdx >= geneticMap.size()) {
      fprintf(stderr, "ERROR: read chrom %s from interference file, but last genetic map chromosome\n",
	      chrom);
      fprintf(stderr, "       is %s\n", geneticMap[chrIdx-1].first);
      exit(5);
    }

    // read remaining tokens:
    for(int i = 0; i < 2; i++) {
      nuStr[i] = strtok_r(NULL, delim, &saveptr);
      pStr[i] = strtok_r(NULL, delim, &saveptr);

      nu[i] = strtod(nuStr[i], &endptr);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: chrom %s, could not parse %s interference nu parameter\n",
		chrom, (i == 0) ? "male" : "female");
	if (errno != 0)
	  perror("strtod");
	exit(5);
      }
      p[i] = strtod(pStr[i], &endptr);
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: chrom %s, could not parse %s interference p parameter\n",
		chrom, (i == 0) ? "male" : "female");
	if (errno != 0)
	  perror("strtod");
	exit(5);
      }
    }

    char *tok;
    if ((tok = strtok_r(NULL, delim, &saveptr)) != NULL) {
      fprintf(stderr, "ERROR: read extra token %s in interference file (chrom %s)\n",
	      tok, chrom);
      exit(5);
    }

    if (strcmp(chrom, geneticMap[chrIdx].first) != 0) {
      fprintf(stderr, "ERROR: order of interference chromosomes different from genetic map:\n");
      fprintf(stderr, "       expected chromosome %s in interference file, read %s\n",
	      geneticMap[chrIdx].first, chrom);
      exit(10);
    }

    // Get the genetic lengths of the male and female maps for this chromosome
    // This is the last element in the genetic maps for the corresponding
    // chromosome:
    double len[2];
    for(int i = 0; i < 2; i++)
      len[i] = geneticMap[chrIdx].second->back().mapPos[i] / 100; // in Morgans
    coIntf.emplace_back(nu, p, len);
    chrIdx++;
  }

  if (chrIdx != geneticMap.size()) {
    fprintf(stderr, "ERROR: read %u chromosomes from interference file, but genetic map has %lu\n",
	    chrIdx, geneticMap.size());
    exit(5);
  }

  free(buffer);
  fclose(in);
}

// Gives the default parent assignment for any branches that have not had
// their parents explicitly specified.
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				int *&thisGenBranchParents,
				int *prevGenSpouseNum,
				vector<bool> *branchParentsAssigned) {
  // how many new branches is each previous branch the parent of?
  int multFactor = thisGenNumBranches / prevGenNumBranches;
  if (multFactor == 0)
    multFactor = 1; // for branches that survive, map prev branch i to cur i

  // allocate space to store the parents of each branch
  if (thisGenBranchParents == NULL)
    thisGenBranchParents = new int[2 * thisGenNumBranches];

  for(int prevB = 0; prevB < prevGenNumBranches; prevB++) {
    if (prevB >= thisGenNumBranches)
      break; // defined all the branches for this generation
    // same founder spouse for all <multFactor> branches that <prevB>
    // is the parent of
    int spouseNum = -1;
    bool spouseNumDefined = false;

    for(int multB = 0; multB < multFactor; multB++) {
      int curBranch = prevB * multFactor + multB;
      if (branchParentsAssigned && (*branchParentsAssigned)[curBranch])
	continue; // skip assignment of branches that were assigned previously
      thisGenBranchParents[curBranch*2] = prevB;
      if (!spouseNumDefined) {
	if (prevGenSpouseNum) {
	  prevGenSpouseNum[ prevB ]--;
	  spouseNum = prevGenSpouseNum[ prevB ];
	}
	else
	  spouseNum = -1;
	spouseNumDefined = true;
      }
      thisGenBranchParents[curBranch*2 + 1] = spouseNum;
    }
  }
  // For any branches in this generation that are not an exact multiple of
  // the number of branches in the previous generation, make them brand new
  // founders. They will contain exactly one person regardless of the number
  // of samples to print
  for(int newB = prevGenNumBranches * multFactor; newB < thisGenNumBranches;
								      newB++) {
    if (branchParentsAssigned && (*branchParentsAssigned)[newB])
      continue; // skip assignment of branches that were assigned previously
    // undefined parents for excess branches: new founders
    thisGenBranchParents[newB*2] = thisGenBranchParents[newB*2 + 1] = -1;
  }
}

// Reads in branch parent specification and makes the parent assignments
void readBranchParents(int prevGenNumBranches, int thisGenNumBranches,
		       int *&thisGenBranchParents, int *&prevGenSexConstraints,
		       int *&prevGenSpouseNum,
		       vector<bool> &branchParentsAssigned,
		       vector< boost::dynamic_bitset<>* > &spouseDependencies,
		       const char *delim, char *&saveptr, char *&endptr,
		       int line) {
  assert(prevGenSpouseNum == NULL);
  prevGenSpouseNum = new int[prevGenNumBranches];
  for(int b = 0; b < prevGenNumBranches; b++)
    // What number have we assigned through for founder spouses of individuals
    // in the previous generation? Note that founders have an id (in the code
    // for the purposes of <curBranchParents>) that are always negative and
    // that by default we assign one founder spouse to marry one person in
    // each branch in the previous generation (see assignDefaultBranchParents())
    prevGenSpouseNum[b] = 0;

  thisGenBranchParents = new int[2 * thisGenNumBranches];

  spouseDependencies.clear();
  prevGenSexConstraints = new int[prevGenNumBranches];
  for(int i = 0; i < prevGenNumBranches; i++)
    prevGenSexConstraints[i] = -1;

  // so far, all the branches in the current generation are assigned default
  // parents; track which branches get explicitly assigned and throw an error
  // if the same branch is assigned more than once
  branchParentsAssigned.clear();
  for(int i = 0; i < thisGenNumBranches; i++)
    branchParentsAssigned.push_back(false);

  while (char *assignToken = strtok_r(NULL, delim, &saveptr)) {
    char *assignBranches = assignToken; // will add '\0' at ':'
    int i;

    // split on ':' to get the parent assignments on the right and the
    // branches on the left
    for(i = 0; assignToken[i] != ':' && assignToken[i] != '\0'; i++);
    if (assignToken[i] != ':') {
      fprintf(stderr, "ERROR: line %d in def: improperly formatted parent assignment field %s\n",
	  line, assignToken);
      exit(8);
    }
    assignToken[i] = '\0';

    // Get the one or two parents
    char *assignPar[2];
    assignPar[0] = &(assignToken[i+1]); // will add '\0' at '_' if present
    assignPar[1] = NULL; // initially; updated just below

    // Find the second parent if present
    for(i = 0; assignPar[0][i] != '_' && assignPar[0][i] != '\0'; i++);
    if (assignPar[0][i] == '_') {
      assignPar[0][i] = '\0';
      assignPar[1] = &(assignPar[0][i+1]);
    }

    int parIdx[2] = { -1, -1 };
    for(int p = 0; p < 2 && assignPar[p] != NULL && assignPar[p][0] != '\0';
									  p++) {
      parIdx[p] = strtol(assignPar[p], &endptr, 10) - 1; // 0 indexed => -1
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: line %d in def: unable to parse parent assignment for branches %s\n",
	    line, assignBranches);
	if (errno != 0)
	  perror("strtol");
	exit(2);
      }
      if (parIdx[p] < 0) {
	fprintf(stderr, "ERROR: line %d in def: parent assignments must be of positive branch numbers\n",
	    line);
	exit(8);
      }
      else if (parIdx[p] >= prevGenNumBranches) {
	fprintf(stderr, "ERROR: line %d in def: parent branch number %d is more than the number of\n",
		line, parIdx[p]+1);
	fprintf(stderr, "       branches (%d) in the previous generation\n",
		prevGenNumBranches);
	exit(8);
      }
    }
    if (parIdx[0] == -1) {
      // new founder
      assert(parIdx[1] == -1);
    }
    else if (parIdx[1] == -1) {
      // Have not yet assigned the numerical id of parent 1. Because the def
      // file doesn't specify this, it is a founder, and one that hasn't been
      // assigned before. As such, we'll get a unique number associated with a
      // spouse of parIdx[0]. Negative values correspond to founders, so we
      // decrement <prevGenSpouseNum>. (It is initialized to -1 above)
      prevGenSpouseNum[ parIdx[0] ]--;
      parIdx[1] = prevGenSpouseNum[ parIdx[0] ];
    }
    else {
      if (parIdx[0] == parIdx[1]) {
	fprintf(stderr, "ERROR: line %d in def: cannot have both parents be from same branch\n",
		line);
	exit(8);
      }
      updateSexConstraints(prevGenSexConstraints, parIdx, prevGenNumBranches,
			   spouseDependencies, line);
    }

    // so that we can print the parent assignment in case of errors below
    if (assignPar[1] != NULL)
      assignPar[1][-1] = '_';
    char *fullAssignPar = assignPar[0];

    // process the branches to be assigned <parIdx> as parents
    bool done = false;
    // the starting branch for a range (delimited by '-'); see below
    char *startBranch = NULL;
    while (!done) {
      for(i = 0; assignBranches[i] != ',' && assignBranches[i] != '-' &&
		 assignBranches[i] != '\0'; i++);
      if (assignBranches[i] == '-') { // have a range; just passed over start:
	assignBranches[i] = '\0';
	if (startBranch != NULL) {
	  fprintf(stderr, "ERROR: line %d in def: improperly formatted branch range \"%s-%s-\"\n",
		  line, startBranch, assignBranches);
	  exit(5);
	}
	startBranch = assignBranches;
	assignBranches = &(assignBranches[i+1]); // go through next loop
      }
      if (assignBranches[i] == ',' || assignBranches[i] == '\0') {
	if (assignBranches[i] == '\0')
	  done = true;
	else
	  assignBranches[i] = '\0';

	int curBranch = strtol(assignBranches, &endptr, 10) - 1; // 0 indexed
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: line %d in def: unable to parse branch %s to assign parent %s to\n",
		  line, assignBranches, fullAssignPar);
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}

	if (startBranch) {
	  int rangeEnd = curBranch;
	  int rangeStart = strtol(startBranch, &endptr, 10) - 1; // 0 indexed
	  if (errno != 0 || *endptr != '\0') {
	    fprintf(stderr, "ERROR: line %d in def: unable to parse branch %s to assign parent %s to\n",
		    line, startBranch, fullAssignPar);
	    if (errno != 0)
	      perror("strtol");
	    exit(2);
	  }
	  startBranch = NULL; // parsed: reset this variable

	  if (rangeStart >= rangeEnd) {
	    fprintf(stderr, "ERROR: line %d in def: assigning parents to non-increasing branch range %d-%d\n",
		    line, rangeStart, rangeEnd);
	    exit(8);
	  }

	  for(int branch = rangeStart; branch <= rangeEnd; branch++) {
	    if (branchParentsAssigned[branch]) {
	      fprintf(stderr, "ERROR: line %d in def: parents of branch number %d assigned multiple times\n",
		      line, branch+1);
	      exit(8);
	    }
	    branchParentsAssigned[branch] = true;
	    for(int p = 0; p < 2; p++)
	      thisGenBranchParents[branch*2 + p] = parIdx[p];
	  }
	}
	else {
	  if (branchParentsAssigned[curBranch]) {
	    fprintf(stderr, "ERROR: line %d in def: parents of branch number %d assigned multiple times\n",
		    line, curBranch+1);
	    exit(8);
	  }
	  branchParentsAssigned[curBranch] = true;
	  for(int p = 0; p < 2; p++)
	    thisGenBranchParents[curBranch*2 + p] = parIdx[p];
	}

	assignBranches = &(assignBranches[i+1]); // go through next loop
      }
    }

    if (startBranch != NULL) {
      fprintf(stderr, "ERROR: line %d in def: range of branches to assign parents does not terminate\n",
	      line);
      exit(8);
    }
  }

  // Now update <prevGenSexConstraints> array to have integer values
  // such that all individuals with the same integer will be assigned the same
  // sex. Even numbers will be assigned a random sex and that number+1 will be
  // assigned the opposite sex.
  //
  // Which individuals have we already updated the index of?
  boost::dynamic_bitset<> updated(prevGenNumBranches);
  assert(spouseDependencies.size() % 2 == 0);
  int nextIndex = 0;
  for(unsigned int i = 0; i < spouseDependencies.size(); i += 2) {
    if ((updated & *spouseDependencies[i]).none()) {
      // should also hold for the spousal set:
      assert((updated & *spouseDependencies[i+1]).none());

      // haven't yet updated <prevGenSexAssignment> for individuals in
      // <spouseDependencies[i,i+1]> -- do the update for each:
      for(int p = 0; p < 2; p++) { // for the two pairs of samples
	for(unsigned int samp = spouseDependencies[i+p]->find_first();
		samp < spouseDependencies[i+p]->size();
		samp = spouseDependencies[i+p]->find_next(samp)) {
	  prevGenSexConstraints[samp] = nextIndex;
	}

	updated |= *spouseDependencies[i+p]; // have now updated these samples
	nextIndex++; // different index for the next set
      }
    }
    else {
      // Ensure consistency:
      for(int p = 0; p < 2; p++)
	assert((updated & *spouseDependencies[i+p]) ==
						      *spouseDependencies[i+p]);
      // remove the pair of elements from spouseDependencies
      spouseDependencies.erase(spouseDependencies.begin() + i);
      // note: might think we should erasa(begin() + i+1), but because of the
      // the erase() command just above, that element is now at begin() + i
      spouseDependencies.erase(spouseDependencies.begin() + i);
      // so that we don't skip over the elements that are now at positions i,i+1
      i -= 2;
    }
  }
  
  // don't need the sets anymore:
  for(unsigned int i = 0; i < spouseDependencies.size(); i++)
    delete spouseDependencies[i];

  assignDefaultBranchParents(prevGenNumBranches, thisGenNumBranches,
			     thisGenBranchParents, prevGenSpouseNum,
			     &branchParentsAssigned);
}

// Given the branch indexes of two parents, adds constraints and error checks
// to ensure that this couple does not violate the requirement that parents
// must have opposite sex.
void updateSexConstraints(int *&prevGenSexConstraints, int parIdx[2],
			  int prevGenNumBranches,
			  vector< boost::dynamic_bitset<>*> &spouseDependencies,
			  int line) {
  // we check these things in the caller, but just to be sure:
  assert(parIdx[0] >= 0 && parIdx[1] >= 0);
  assert(parIdx[0] < prevGenNumBranches && parIdx[1] < prevGenNumBranches);

  if (prevGenSexConstraints[ parIdx[0] ] == -1 &&
      prevGenSexConstraints[ parIdx[1] ] == -1) {
    // neither is a member of a spouse set: create and add to
    // <spouseDependencies>
    boost::dynamic_bitset<> *sets[2];
    for(int p = 0; p < 2; p++) {
      sets[p] = new boost::dynamic_bitset<>(prevGenNumBranches);
      sets[p]->set( parIdx[p] );
      // which index in spouseDependencies is the set corresponding to this
      // parent stored in? As we're about to add it, just below this, the
      // current size will be the index
      prevGenSexConstraints[ parIdx[p] ] = spouseDependencies.size();
      spouseDependencies.push_back(sets[p]);
    }
  }
  else if (prevGenSexConstraints[ parIdx[0] ] == -1 ||
	   prevGenSexConstraints[ parIdx[1] ] == -1) {
    // one spouse is a member of a spouse set and the other is not: add the
    // other parent to the opposite spouse set, error check, and update state
    int assignedPar = -1;
    if (prevGenSexConstraints[ parIdx[0] ] >= 0)
      assignedPar = 0;
    else
      assignedPar = 1;

    int otherPar = assignedPar ^ 1;
    int assignedSetIdx = prevGenSexConstraints[ parIdx[ assignedPar ] ];

    // since <otherPar> isn't in a spouse set yet, it definitely shouldn't be in
    // the same spouse set as <assignedPar>
    assert(!spouseDependencies[ assignedSetIdx ]->test( parIdx[otherPar] ));

    // NOTE: The following is a trick that relies on the fact that sets are
    // stored as sequential pairs in <spouseDependencies>. Ex: indexes 0 and 1
    // are pairs of spouses that are dependent upon each other (those in index
    // 0 must have the same sex and must be opposite those in index 1).  To get
    // an even number that is 1 minus an odd value or the odd number that is 1
    // plus an even value, it suffices to flip the lowest order bit:
    int otherSetIdx = assignedSetIdx ^ 1;
    // add <parIdx[otherPar]> to the set containing spouses of
    // <parIdx[assignedPar]>
    spouseDependencies[ otherSetIdx ]->set( parIdx[otherPar] );
    prevGenSexConstraints[ parIdx[otherPar] ] = otherSetIdx;
  }
  else {
    assert(prevGenSexConstraints[ parIdx[0] ] >= 0 &&
				      prevGenSexConstraints[ parIdx[1] ] >= 0);
    // both are members of spouse sets

    if (prevGenSexConstraints[ parIdx[0] ] / 2 ==
				      prevGenSexConstraints[ parIdx[1] ] / 2) {
      // spouse sets from the same pair of sets
      if (prevGenSexConstraints[ parIdx[0] ] ==
					  prevGenSexConstraints[ parIdx[1] ]) {
	fprintf(stderr, "ERROR: line %d in def: assigning %d and %d as parents is impossible due to\n",
	    line, parIdx[0]+1, parIdx[1]+1);
	fprintf(stderr, "       other parent assignments: they necessarily have same sex\n");
	exit(3);
      }
      // otherwise done: are already members of sets that are pairs and
      // therefore will be constrained to be opposite
    }
    else {
      // two distinct sets of spouse sets: must generate two sets that are the
      // unions of the appropriate sets -- all the spouses of both must be
      // constrained to be opposite sex
      int setIdxes[2][2]; // current set indexes
      for(int p = 0; p < 2; p++) {
	setIdxes[p][0] = prevGenSexConstraints[ parIdx[p] ];
	// see comment denoted with NOTE a little ways above for why this works:
	setIdxes[p][1] = setIdxes[p][0] ^ 1;
      }

      // ensure the values are different pointers (if this process already
      // happened they could be the same pointers)
      if (spouseDependencies[ setIdxes[0][0] ] !=
					spouseDependencies[ setIdxes[1][1] ]) {
	assert(spouseDependencies[ setIdxes[0][1] ] !=
					  spouseDependencies[ setIdxes[1][0] ]);
	// take the union such that spouses of each are contained in both
	*spouseDependencies[ setIdxes[0][0] ] |=
					  *spouseDependencies[ setIdxes[1][1] ];
	*spouseDependencies[ setIdxes[0][1] ] |=
					  *spouseDependencies[ setIdxes[1][0] ];

	// intersection of these should be empty, otherwise there's
	// inconsistencies and individuals that are meant to be different sexes
	// simultaneously
	if ((*spouseDependencies[ setIdxes[0][0] ] &
				*spouseDependencies[ setIdxes[0][1] ]).any()) {
	  fprintf(stderr, "ERROR: line %d in def: assigning %d and %d as parents is impossible due to\n",
		  line, parIdx[0]+1, parIdx[1]+1);
	  fprintf(stderr, "       other parent assignments: they necessarily have same sex\n");
	  exit(4);
	}

	// replace pointers stored in setIdxes[1][ 0,1 ] to those in
	// setIdxes[0][ 0,1 ]. First free the objects already stored there, then
	// update
	for(int p = 0; p < 2; p++) {
	  // replace all pointers to <spouseDependencies[ setIdxes[1][p] ]> to
	  // point to <spouseDependencies[ setIdxes[0][ p^1 ] ]:
	  boost::dynamic_bitset<> *indexesToReplace =
					  spouseDependencies[ setIdxes[1][p] ];
	  for(unsigned int samp = indexesToReplace->find_first();
		  samp < indexesToReplace->size();
		  samp = indexesToReplace->find_next(samp)) {
	    spouseDependencies[ samp ] = spouseDependencies[ setIdxes[0][p^1] ];
	  }

	  delete indexesToReplace;
	}
      }
      else {
	// pointers already equal; ensure consistency among all sets:
	// everything is already done:
	assert(spouseDependencies[ setIdxes[0][1] ] ==
					  spouseDependencies[ setIdxes[1][0] ]);
	assert((*spouseDependencies[ setIdxes[0][0] ] &
				*spouseDependencies[ setIdxes[0][1] ]).none());
      }
    }
  }
}

// Simulate data for each specified pedigree type for the number of requested
// families. Returns the number of founder haplotypes used to produce these
// simulated samples.
int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     bool sexSpecificMaps, vector<COInterfere> &coIntf) {
  // Note: throughout we use 0-based values for generations though the input
  // def file is 1-based

  int totalFounderHaps = 0;
  // Stores the random assignments of sex for the parents in a given generation
  // This relates to the constraints on sexes of individuals. Individuals with
  // the same index -- an index into the <sexAssignments> vector -- will be
  // assigned the same sex. And any pair of indexes i and i^1 will have
  // opposite sex.
  vector<int> sexAssignments;

  theSamples = new Person****[simDetails.size()];
  for(unsigned int ped = 0; ped < simDetails.size(); ped++) { // for each ped
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numSampsToRetain = simDetails[ped].numSampsToRetain;
    int *numBranches = simDetails[ped].numBranches;
    int **branchParents = simDetails[ped].branchParents;
    int **sexConstraints = simDetails[ped].sexConstraints;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    ////////////////////////////////////////////////////////////////////////////
    // Allocate space and make Person objects for all those we will simulate,
    // assigning sex if <sexSpecificMaps> is true
    theSamples[ped] = new Person***[numFam];
    for (int fam = 0; fam < numFam; fam++) {

      theSamples[ped][fam] = new Person**[numGen];
      for(int curGen = 0; curGen < numGen; curGen++) {

	theSamples[ped][fam][curGen] = new Person*[ numBranches[curGen] ];

	// ready to make assignments for this generaiton
	sexAssignments.clear();

	// allocate Persons for each branch of <curGen>, assign their sex,
	// and do the simulation for these allocated individuals
	for(int branch = 0; branch < numBranches[curGen]; branch++) {
	  // Determine how many founders and non-founders we need data for in
	  // <branch>:
	  int numFounders, numNonFounders;
	  getPersonCounts(curGen, numGen, branch, numSampsToRetain,
			  branchParents, branchNumSpouses, numFounders,
			  numNonFounders);

	  // Will use convention that all founder spouses are stored as indexes
	  // 0 through <branchNumSpouses>, the "primary" samples is next and
	  // all other non-founders follow.
	  // Note that when the branch is new and has no parents, all
	  // individuals are founders.
	  int numPersons = numFounders + numNonFounders;
	  theSamples[ped][fam][curGen][branch] = new Person[numPersons];

	  if (sexSpecificMaps) {
	    int branchAssign;
	    if (sexConstraints[curGen] == NULL ||
					sexConstraints[curGen][branch] == -1) {
	      // no dependencies, just pick randomly
	      branchAssign = coinFlip(randomGen);
	    }
	    else {
	      while (sexConstraints[curGen][branch] >=
						  (int) sexAssignments.size()) {
		int rand = coinFlip(randomGen);
		sexAssignments.push_back(rand);
		sexAssignments.push_back(rand ^ 1);
	      }
	      branchAssign = sexAssignments[ sexConstraints[curGen][branch] ];
	    }
	    // How many spouses for this branch?
	    int thisBranchNumSpouses;
	    if (branchNumSpouses[curGen])
	      thisBranchNumSpouses = -branchNumSpouses[curGen][branch];
	    else if (curGen == numGen - 1) // no spouses in last generation
	      thisBranchNumSpouses = 0;
	    else                           // one spouse by default
	      thisBranchNumSpouses = 1;
	    for(int ind = 0; ind < thisBranchNumSpouses; ind++) {
	      theSamples[ped][fam][curGen][branch][ind].sex = 1 ^ branchAssign;
	    }
	    // sex of "primary" person -- who each of the above individuals have
	    // children with -- index just after the spouses
	    int primaryIdx = thisBranchNumSpouses;
	    theSamples[ped][fam][curGen][branch][primaryIdx].sex = branchAssign;
	  }

	  /////////////////////////////////////////////////////////////////////
	  // All samples allocated for this pedigree/family/generation/branch:
	  // simulate the actual samples

	  // for each chromosome:
	  for(unsigned int chrIdx = 0; chrIdx < geneticMap.size(); chrIdx++) {
	    vector<PhysGeneticPos> *curMap = geneticMap[chrIdx].second;

	    // Make trivial haplotypes for founders in current generation:
	    // no recombinations in founders
	    Segment trivialSeg;
	    trivialSeg.endPos = curMap->back().physPos;

	    // Simulate the founders for this chromosome:
	    if (curGen != numGen - 1) { // no founders in the last generation
	      for(int ind = 0; ind < numFounders; ind++) {
		for(int h = 0; h < 2; h++) { // 2 founder haplotypes per founder
		  if (chrIdx == 0)
		    trivialSeg.foundHapNum = totalFounderHaps++;
		  else
		    // want the same founder on all chromosomes, so access
		    // the haplotype number assigned to the previous
		    // chromosome for this person:
		    trivialSeg.foundHapNum =
			theSamples[ped][fam][curGen][branch][ind].haps[h].
						      back().back().foundHapNum;

		  // the following copies <trivialSeg>, so we can reuse it
		  theSamples[ped][fam][curGen][branch][ind].haps[h].
								 emplace_back();
		  theSamples[ped][fam][curGen][branch][ind].haps[h].back().
							  push_back(trivialSeg);
		}
	      }
	    }

	    if (numNonFounders == 0) {
	      assert(curGen == 0 || branchParents[curGen][branch*2] == -1);
	      continue; // no non-founders in first generation
	    }

	    assert(curGen > 0 && branchParents[curGen][branch*2] >= 0);

	    // Now simulate the non-founders in <branch> of <curGen>
	    //
	    // First figure out who the parents are:
	    int parBrch[2]; // which branch are the two parents in
	    int parIdx[2];  // index of the Person in the branch?
	    for(int p = 0; p < 2; p++) {
	      parBrch[p] = branchParents[curGen][branch*2 + p];
	      if (parBrch[p] < 0) {
		assert(p == 1);
		// founders have negative indexes that start from -1, so we
		// add 1 to get it to be 0 based and negate to get the index
		parIdx[p] = -(parBrch[p] + 1);
		parBrch[1] = parBrch[0];
	      }
	      else {
		// use "primary" person: immediately after all the spouses
		int thisBranchNumSpouses;
		if (branchNumSpouses[curGen-1])
		  thisBranchNumSpouses =-branchNumSpouses[curGen-1][parBrch[p]];
		else if (curGen == numGen - 1) // no spouses in last generation
		  thisBranchNumSpouses = 0;
		else                           // one spouse by default
		  thisBranchNumSpouses = 1;
		parIdx[p] = thisBranchNumSpouses;
	      }
	    }

	    Person **prevGenSamps = theSamples[ped][fam][curGen-1];
	    Person **curGenSamps = theSamples[ped][fam][curGen];
	    // the non-founders are stored just after the founders
	    for(int ind = numFounders; ind < numPersons; ind++) {
	      // If we're using sex-specific maps, the two parents' sexes should
	      // differ
	      if (sexSpecificMaps) {
		assert(prevGenSamps[ parBrch[0] ][ parIdx[0] ].sex !=
				   prevGenSamps[ parBrch[1] ][ parIdx[1] ].sex);
	      }

	      for(int p = 0; p < 2; p++) { // meioses from each parent index <p>
		Person &theParent = prevGenSamps[ parBrch[p] ][ parIdx[p] ];

		int hapIdx; // haplotype index for the simulated sample
		if (sexSpecificMaps)
		  // match the sex of the parent if using sex-specific maps
		  hapIdx = theParent.sex;
		else
		  hapIdx = p;

		curGenSamps[branch][ind].haps[hapIdx].emplace_back();
		generateHaplotype(curGenSamps[branch][ind].haps[hapIdx].back(),
				  theParent, curMap, coIntf, chrIdx);
	      } // <parIdx> (simulate each transmitted haplotype)
	    } // <ind>
	  } // <geneticMap> (chroms)
	} // <branch>
      } // <curGen>
    } // <fam>

  } // <ped>

  return totalFounderHaps;
}

// Returns (via parameters) the number of founders and non-founders in the given
// generation and branch.
void getPersonCounts(int curGen, int numGen, int branch, int *numSampsToRetain,
		     int **branchParents, int **branchNumSpouses,
		     int &numFounders, int &numNonFounders) {
  if (curGen > 0 && branchParents[curGen][branch*2] >= 0) { // have parent(s)?
    numNonFounders = numSampsToRetain[curGen];
    if (numNonFounders == 0)
      // not saving, but need parent of next generation: the "primary" person
      numNonFounders = 1;

    // set default number of founders (spouses):
    if (curGen == numGen - 1) // none in last generation
      numFounders = 0;
    else
      numFounders = 1;
    // more founders depending on the number of spouses in the current branch
    if (branchNumSpouses[curGen])
      // We store the number of branch spouses as negative: must negate
      numFounders = -branchNumSpouses[curGen][branch];
  }
  else {
    assert(curGen == 0 || branchParents[curGen][branch*2 + 1] == -1);

    // no parents, so only founders in this branch
    numNonFounders = 0;
    numFounders = 1; // one founder initialy: the "primary" person
    if (branchNumSpouses[curGen])
      // We store the number of branch spouses as negative: must negate
      // + 1 because the "main" person in this branch is also a founder:
      // that person is married to the <-prevGenSpouseNum[curGen][branch]>
      // other founders
      numFounders = -branchNumSpouses[curGen][branch] + 1;
    else if (curGen != numGen - 1)
      // though we haven't initialized <branchNumSpouses>,
      // will necessarily have 1 spouse (except for last generation)
      numFounders++;
  }
}

// Simulate one haplotype <toGenerate> by sampling crossovers and switching
// between the two haplotypes stored in <parent>. Uses the genetic map stored
// in <curMap> which is expected to correspond to the chromosome being
// simulated and may contain either one sex-averaged map or a male and female
// map.
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap,
		       vector<COInterfere> &coIntf, unsigned int chrIdx) {
  // For the two haplotypes in <parent>, which segment index (in
  // parent.haps[].back()) is the current <switchMarker> position contained in?
  unsigned int curSegIdx[2] = { 0, 0 };

  int mapIdx = parent.sex; // either 0 or 1 for sex-averaged/male or female

  // Pick haplotype for the beginning of the transmitted one:
  int curHap = coinFlip(randomGen);

  double firstcMPos = curMap->front().mapPos[mapIdx];
  double lastcMPos = curMap->back().mapPos[mapIdx];
  // Get the genetic length of this chromosome (for the appropriate sex).
  // This is the last element in the genetic maps for the corresponding
  // chromosome:
  // In centiMorgans:
  double chrLength = lastcMPos - firstcMPos;
  // Locations of the crossovers
  vector<double> coLocations; // in Morgans

  if (CmdLineOpts::interfereFile) {
    coIntf[chrIdx].simStahl(coLocations, parent.sex, randomGen);
  }
  else {
    double lastPos = 0.0; // position of last crossover
    while (true) { // simulate until crossover is past chromosome end
      double curPos = lastPos + crossoverDist(randomGen);
      if (curPos >= chrLength / 100)
	break;
      coLocations.push_back(curPos);
      lastPos = curPos;
    }
  }

  // initially assume we'll recombine between first and positions with map info
  int switchIdx = 1;

  int mapNumPos = curMap->size();
  for(auto it = coLocations.begin(); it != coLocations.end(); it++) {
    // Multiply by 100 to get cM:
    double cMPosNextCO = firstcMPos + (*it * 100);

    // TODO: slow linear search for switching index, but probably fast enough
    for( ; switchIdx < mapNumPos; switchIdx++) {
      if ((*curMap)[switchIdx].mapPos[mapIdx] > cMPosNextCO)
	break;
    }
    if (switchIdx == mapNumPos)
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
    for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
      Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
      if (seg.endPos >= switchPos) {
	// last segment to copy, and we will break it at <switchPos>
	Segment copy = seg; // don't modify <seg.endPos> directly
	copy.endPos = switchPos;
	toGenerate.push_back(copy);
	if (seg.endPos == switchPos)
	  curSegIdx[curHap]++;
	break; // done copying
      }
      else {
	toGenerate.push_back(seg);
      }
    }
    assert(curSegIdx[curHap] < parent.haps[curHap][chrIdx].size());

    // swap haplotypes
    curHap ^= 1;
    // must update <curSegIdx[curHap]>
    for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
      Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
      if (seg.endPos > switchPos)
	// current segment spans from just after <switchMarker> to <endMarker>
	break;
    }
    assert(curSegIdx[curHap] < parent.haps[curHap][chrIdx].size());
  }

  // copy through to the end of the chromosome:
  for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
    Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
    toGenerate.push_back(seg);
  }
}

// Print the break points to <outFile>
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *bpFile) {
  FILE *out = fopen(bpFile, "w");
  if (!out) {
    printf("ERROR: could not open output file %s!\n", bpFile);
    perror("open");
    exit(1);
  }

  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numSampsToRetain = simDetails[ped].numSampsToRetain;
    int *numBranches = simDetails[ped].numBranches;
    int **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;
    char *pedName = simDetails[ped].name;

    for(int fam = 0; fam < numFam; fam++) {
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  if (numSampsToRetain[gen] > 0) {
	    int numNonFounders, numFounders;
	    getPersonCounts(gen, numGen, branch, numSampsToRetain,
			    branchParents, branchNumSpouses, numFounders,
			    numNonFounders);
	    int numPersons = numNonFounders + numFounders;

	    for(int ind = 0; ind < numPersons; ind++) {
	      for(int h = 0; h < 2; h++) {
		int sex = theSamples[ped][fam][gen][branch][ind].sex;
		int thisBranchNumSpouses;
		if (branchNumSpouses[gen])
		  thisBranchNumSpouses = -branchNumSpouses[gen][branch];
		else if (gen == numGen - 1) // no spouses in last generation
		  thisBranchNumSpouses = 0;
		else                        // one spouse by default
		  thisBranchNumSpouses = 1;
		if (ind < thisBranchNumSpouses)
		  fprintf(out, "%s%d_g%d-b%d-s%d s%d h%d", pedName, fam+1,
			  gen+1, branch+1, ind+1, sex, h);
		else
		  fprintf(out, "%s%d_g%d-b%d-i%d s%d h%d", pedName, fam+1,
			  gen+1, branch+1, ind - thisBranchNumSpouses + 1,
			  sex, h);

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

// open <filename> either using standard I/O or as a gzipped file if <isGZ>
bool fileOrGZ_open(fileOrGZ &fgz, const char *filename, const char *mode,
		   bool isGZ) {
  fgz.isGZ = isGZ;
  if (isGZ) { // gzipped
    fgz.fp = NULL;
    fgz.gfp = gzopen(filename, mode);
    if (!fgz.gfp)
      return false;
    else
      return true;
  }
  else { // standard
    fgz.gfp = NULL;
    fgz.fp = fopen(filename, mode);
    if (!fgz.fp)
      return false;
    else
      return true;
  }
}

// NOTE: unlike getline() this assumes that *lineptr is non-NULL, and that
// n > 0.
int fileOrGZ_getline(char **lineptr, size_t *n, fileOrGZ &fgz) {
  if (fgz.isGZ) {
    char *buf = *lineptr;
    int n_read = 0;
    int c;

    while ((c = gzgetc(fgz.gfp)) != EOF) {
      // About to have read one more, so n_read + 1 needs to be less than *n.
      // Note that we use >= not > since we need one more space for '\0'
      if (n_read + 1 >= (int) *n) {
	const size_t GROW = 1024;
	buf = (char *) realloc(*lineptr, *n + GROW);
	if (buf == NULL) {
	  fprintf(stderr, "ERROR: out of memory!\n");
	  exit(1);
	}
	*n += GROW;
	*lineptr = buf;
      }
      buf[n_read] = (char) c;
      n_read++;
      if (c == '\n')
	break;
    }

    if (c == EOF && n_read == 0)
      return -1;

    buf[n_read] = '\0';

    return n_read;
  }
  else {
    return getline(lineptr, n, fgz.fp);
  }
}

int fileOrGZ_printf(fileOrGZ &fgz, const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);
  if (fgz.isGZ) {
    char *buf;
    ret = vasprintf(&buf, format, args);
    if (ret < 0)
      return ret;
    size_t len = strlen(buf);
    ret = gzwrite(fgz.gfp, buf, len);
    free(buf);
    return ret;
  }
  else {
    ret = vfprintf(fgz.fp, format, args);
  }
  va_end(args);
  return ret;
}

int fileOrGZ_close(fileOrGZ &fgz) {
  if (fgz.isGZ) {
    return gzclose(fgz.gfp);
  }
  else {
    return fclose(fgz.fp);
  }
}

// Given the simulated break points for individuals in each pedigree/family
// stored in <theSamples> and other necessary information, reads input VCF
// format data from the file named <inVCFfile> and prints the simulated
// haplotypes for each sample to <outVCFfile> in VCF format.
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, char *inVCFfile, char *outFileBuf,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     FILE *outs[2]) {
  int inVCFlen = strlen(inVCFfile);
  bool isGZ = false;
  if (strcmp(&CmdLineOpts::inVCFfile[ inVCFlen - 3 ], ".gz") == 0) {
    isGZ = true;
    sprintf(outFileBuf, "%s.vcf.gz", CmdLineOpts::outPrefix);
  }
  else {
    sprintf(outFileBuf, "%s.vcf", CmdLineOpts::outPrefix);
  }

  // open input VCF file:
  fileOrGZ in;
  bool success = fileOrGZ_open(in, inVCFfile, "r", isGZ);
  if (!success) {
    printf("\nERROR: could not open input VCF file %s!\n", inVCFfile);
    perror("open");
    exit(1);
  }

  // open output VCF file:
  fileOrGZ out;
  success = fileOrGZ_open(out, outFileBuf, "w", isGZ);
  if (!success) {
    printf("\nERROR: could not open output VCF file %s!\n", outFileBuf);
    perror("open");
    exit(1);
  }

  bernoulli_distribution genoErr( CmdLineOpts::genoErrRate );
  bernoulli_distribution homErr( CmdLineOpts::homErrRate );
  bernoulli_distribution setMissing( CmdLineOpts::missRate );

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  // technically tab and newline; we want the latter so that the last sample id
  // on the header line doesn't include the newline character in it
  const char *tab = "\t\n";
  const char *bar = "|";
  // Below when we print the VCF output, we alternate printing tab
  // and / between successive alleles. Make this simpler with:
  char betweenAlleles[2] = { '\t', '/' };
  if (CmdLineOpts::keepPhase)
    betweenAlleles[1] = '|';

  // Have we encountered / printed a warning for markers with more than 2
  // alleles (we don't currently introduce errors at such markers)
  bool alleleCountWarnPrinted = false;

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
  vector<int> extraSamples; // Sample indexes to print for --retain_extra
  // map from haplotype index / 2 to sample_index
  int *founderSamples = new int[totalFounderHaps / 2];
  // number of elements of <extraSamples> to print (see below)
  unsigned int numToRetain = 0;

  while (fileOrGZ_getline(&buffer, &bytesRead, in) >= 0) { // lines of input VCF
    if (buffer[0] == '#' && buffer[1] == '#') {
      // header line: print to output
      fileOrGZ_printf(out, "%s", buffer);
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

      if (numInputSamples < totalFounderHaps / 2) {
	fprintf(stderr, "\nERROR: need %d founders, but input only contains %d\n",
		totalFounderHaps / 2, numInputSamples);
	exit(5);
      }

      // Do math and store sample ids for --retain_extra:
      unsigned int numExtraSamples = numInputSamples - totalFounderHaps / 2;
      bool cantRetainEnough = false;
      if (CmdLineOpts::retainExtra < 0) {
	numToRetain = numExtraSamples;
      }
      else {
	numToRetain = CmdLineOpts::retainExtra;
	if (numToRetain > numExtraSamples) {
	  cantRetainEnough = true;
	  numToRetain = numExtraSamples;
	}
      }

      // Store ids for all extra samples -- those whose shuffled haplotype
      // assignment is after all that will be used:
      for(int i = 0; i < numInputSamples; i++) {
	if (shuffHaps[i] < totalFounderHaps)
	  founderSamples[ shuffHaps[i] / 2 ] = i;
	else
	  extraSamples.push_back(i);
      }
      assert(extraSamples.size() == numExtraSamples);

      // want to randomize which samples get included, though this is only
      // relevant if we have more samples than are requested to be retained:
      if (numToRetain < numExtraSamples) {
	// will print the first <numToRetain> from this list (random subset)
	shuffle(extraSamples.begin(), extraSamples.end(), randomGen);
      }

      // open output ids file (if needed):
      FILE *idOut = NULL;
      if (CmdLineOpts::printFounderIds) {
	sprintf(outFileBuf, "%s.ids", CmdLineOpts::outPrefix);
	idOut = fopen(outFileBuf, "w");
	if (!idOut) {
	  printf("ERROR: could not open found ids file %s!\n", outFileBuf);
	  perror("open");
	  exit(1);
	}
      }

      for(int o = 0; o < 2; o++) {
	fprintf(outs[o], "done.\n"); // initial scan of VCF file (see main())
	fprintf(outs[o], "  Input contains %d samples, using %d as founders, and retaining %d\n",
		numInputSamples, totalFounderHaps / 2, numToRetain);
	if (cantRetainEnough) {
	  fprintf(outs[o], "  Note: cannot retain all requested %d samples\n",
		  CmdLineOpts::retainExtra);
	}
	if (idOut) {
	  fprintf(outs[o], "Generating founder ids file... ");
	  fflush(outs[o]);
	}
      }

      // Now print the header line indicating fields and sample ids for the
      // output VCF
      fileOrGZ_printf(out,
		       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

      // print sample ids:
      for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
	int numFam = simDetails[ped].numFam;
	int numGen = simDetails[ped].numGen;
	int *numSampsToRetain = simDetails[ped].numSampsToRetain;
	int *numBranches = simDetails[ped].numBranches;
	int **branchParents = simDetails[ped].branchParents;
	int **branchNumSpouses = simDetails[ped].branchNumSpouses;
	char *pedName = simDetails[ped].name;

	for(int fam = 0; fam < numFam; fam++)
	  for(int gen = 0; gen < numGen; gen++)
	    for(int branch = 0; branch < numBranches[gen]; branch++)
	      if (numSampsToRetain[gen] > 0 || idOut) { // need to print?
		int numNonFounders, numFounders;
		getPersonCounts(gen, numGen, branch, numSampsToRetain,
				branchParents, branchNumSpouses, numFounders,
				numNonFounders);
		int numPersons = numNonFounders + numFounders;
		for(int ind = 0; ind < numPersons; ind++) {
		  int thisBranchNumSpouses;
		  if (branchNumSpouses[gen])
		    thisBranchNumSpouses = -branchNumSpouses[gen][branch];
		  else if (gen == numGen - 1) // no spouses in last generation
		    thisBranchNumSpouses = 0;
		  else                        // one spouse by default
		    thisBranchNumSpouses = 1;
		  bool curIsFounder = false;
		  if (ind < thisBranchNumSpouses) {
		    curIsFounder = true;
		    if (numSampsToRetain[gen] > 0)
		      fileOrGZ_printf(out, "\t%s%d_g%d-b%d-s%d", pedName, fam+1,
					   gen+1, branch+1, ind+1);
		    if (idOut) // print Ped-sim id to founder id file:
		      fprintf(idOut, "%s%d_g%d-b%d-s%d", pedName, fam+1,
				     gen+1, branch+1, ind+1);
		  }
		  else {
		    int indNum = ind - thisBranchNumSpouses;
		    if (numSampsToRetain[gen] > 0)
		      fileOrGZ_printf(out, "\t%s%d_g%d-b%d-i%d", pedName, fam+1,
					   gen+1, branch+1, indNum + 1);
		    if (idOut &&
			(gen == 0 || branchParents[gen][branch*2] < 0)) {
		      assert(indNum == 0);
		      curIsFounder = true;
		      fprintf(idOut, "%s%d_g%d-b%d-i%d", pedName, fam+1,
				     gen+1, branch+1, indNum + 1);
		    }
		  }

		  if (idOut && curIsFounder) {
		    int hapNum = theSamples[ped][fam][gen][branch][ind].
				      haps[0][/*chrIdx=*/0].front().foundHapNum;
		    assert(hapNum % 2 == 0);
		    int founderIdx = founderSamples[ hapNum / 2 ];
		    fprintf(idOut, "\t%s\n", sampleIds[ founderIdx ]);
		  }

		}
	      }
      }

      // print the ids for the --retain_extra samples:
      for(unsigned int i = 0; i < numToRetain; i++) {
	int sampIdx = extraSamples[i];
	fileOrGZ_printf(out, "\t%s", sampleIds[ sampIdx ]);
      }

      fileOrGZ_printf(out, "\n");

      for(int o = 0; o < 2; o++) {
	if (idOut)
	  fprintf(outs[o], "done.\n");
	fprintf(outs[o], "Generating VCF file... ");
	fflush(outs[o]);
      }

      if (idOut)
	fclose(idOut);

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
    // TODO! ought to deal with non-trivial formats -- what if GT isn't first
    //       field
    char *otherFields[7];
    for(int i = 0; i < 7; i++)
      otherFields[i] = strtok_r(NULL, tab, &saveptr);

    // count the number of alleles present at this variant; generally this is
    // 2, but the number of (comma separated) values in the ALT field gives the
    // exact number
    char *altField = otherFields[2];
    int numAlleles = 2;
    for(int i = 0; altField[i] != '\0'; i++) {
      if (altField[i] == ',')
	numAlleles++;
    }
    if (numAlleles > 2 && CmdLineOpts::genoErrRate > 0.0 &&
						      !alleleCountWarnPrinted) {
      alleleCountWarnPrinted = true;
      fprintf(stderr, "\nWARNING: genotyping error only implemented for markers with 2 alleles\n");
      fprintf(stderr, "         will not introduce errors at any markers with >2 alleles\n");
    }

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

      for(int h = 0; h < 2; h++) {
	if (alleles[h][0] == '.') {
	  fprintf(stderr, "\nERROR: simulator currently requires all positions to be non-missing\n");
	  exit(5);
	}
	hapAlleles[numStored++] = alleles[h];
      }

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
	printf("ERROR: multiple '|' characters in data field\n");
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
    fileOrGZ_printf(out, "%s\t%s", chrom, posStr);
    for(int i = 0; i < 7; i++)
      fileOrGZ_printf(out, "\t%s", otherFields[i]);

    for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
      int numFam = simDetails[ped].numFam;
      int numGen = simDetails[ped].numGen;
      int *numSampsToRetain = simDetails[ped].numSampsToRetain;
      int *numBranches = simDetails[ped].numBranches;
      int **branchParents = simDetails[ped].branchParents;
      int **branchNumSpouses = simDetails[ped].branchNumSpouses;

      for(int fam = 0; fam < numFam; fam++)
	for(int gen = 0; gen < numGen; gen++)
	  for(int branch = 0; branch < numBranches[gen]; branch++)
	    if (numSampsToRetain[gen] > 0) {
	      int numNonFounders, numFounders;
	      getPersonCounts(gen, numGen, branch, numSampsToRetain,
			      branchParents, branchNumSpouses, numFounders,
			      numNonFounders);
	      int numPersons = numNonFounders + numFounders;

	      for(int ind = 0; ind < numPersons; ind++) {

		// set to missing (according to the rate set by the user)?
		if (setMissing( randomGen )) {
		  for(int h = 0; h < 2; h++)
		    fileOrGZ_printf(out, "%c.", betweenAlleles[h]);
		  continue; // done printing genotype data for this sample
		}

		// non-missing genotype: print, possibly with a genotyping error

		// get founder haps for the current sample
		uint32_t curFounderHaps[2];
		for(int h = 0; h < 2; h++) {
		  Haplotype &curHap = theSamples[ped][fam][gen][branch][ind].
								haps[h][chrIdx];
		  while (curHap.front().endPos < pos) {
		    pop_front(curHap);
		  }
		  assert(curHap.front().endPos >= pos);
		  curFounderHaps[h] = curHap.front().foundHapNum;
		}

		// genotyping error?
		if (genoErr( randomGen ) && numAlleles == 2) {
		  int alleles[2]; // integer allele values
		  for(int h = 0; h < 2; h++)
		    // can get character 0 from founderHaps strings: with only
		    // two alleles possible, these strings must have length 1.
		    // converting to an integer is simple: subtract '0'
		    alleles[h] = founderHaps[ curFounderHaps[h] ][0] - '0';

		  if (alleles[0] != alleles[1]) {
		    // heterozygous: choose an allele to alter
		    int alleleToFlip = coinFlip(randomGen);
		    alleles[ alleleToFlip ] ^= 1;
		  }
		  else {
		    // homozygous: determine whether to change to the opposite
		    // homozygote or to a heterozygote
		    if (homErr(randomGen)) {
		      alleles[0] ^= 1;
		      alleles[1] ^= 1;
		    }
		    else {
		      // will flip only one allele so that the sample becomes
		      // heterozygous; randomly choose which
		      int alleleToFlip = coinFlip(randomGen);
		      alleles[ alleleToFlip ] ^= 1;
		    }
		  }

		  for(int h = 0; h < 2; h++)
		    fileOrGZ_printf(out, "%c%d", betweenAlleles[h],alleles[h]);
		}
		else { // no error: print alleles from original haplotypes
		  for(int h = 0; h < 2; h++)
		    fileOrGZ_printf(out, "%c%s", betweenAlleles[h],
				     founderHaps[ curFounderHaps[h] ]);
		}
	      }
	    }
    }
    // print data for the --retain_extra samples:
    for(unsigned int i = 0; i < numToRetain; i++) {
      int sampIdx = extraSamples[i];
      for(int h = 0; h < 2; h++)
	fileOrGZ_printf(out, "%c%s", betweenAlleles[h],
			 hapAlleles[ 2*sampIdx + h ]);
    }

    fileOrGZ_printf(out, "\n");
  }

  free(buffer);
  fileOrGZ_close(out);
  fileOrGZ_close(in);
}

// print fam format file with the pedigree structure of all individuals included
// in the simulation
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
	      char *famFile) {
  // open output fam file:
  FILE *out = fopen(famFile, "w");
  if (!out) {
    printf("ERROR: could not open output fam file %s!\n", famFile);
    perror("open");
    exit(1);
  }

  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numSampsToRetain = simDetails[ped].numSampsToRetain;
    int *numBranches = simDetails[ped].numBranches;
    int **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;
    char *pedName = simDetails[ped].name;

    for(int fam = 0; fam < numFam; fam++) {
      for(int gen = 0; gen < numGen; gen++) {
	Person **prevGenSamps = NULL;
	if (gen > 0)
	  prevGenSamps = theSamples[ped][fam][gen-1];

	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  int numNonFounders, numFounders;
	  getPersonCounts(gen, numGen, branch, numSampsToRetain, branchParents,
			  branchNumSpouses, numFounders, numNonFounders);
	  int numPersons = numNonFounders + numFounders;

	  for(int ind = 0; ind < numPersons; ind++) {

	    // print family id (PLINK-specific) and sample id
	    int thisBranchNumSpouses;
	    if (branchNumSpouses[gen])
	      thisBranchNumSpouses = -branchNumSpouses[gen][branch];
	    else if (gen == numGen - 1) // no spouses in last generation
	      thisBranchNumSpouses = 0;
	    else                        // one spouse by default
	      thisBranchNumSpouses = 1;
	    if (ind < thisBranchNumSpouses)
	      fprintf(out, "%s%d %s%d_g%d-b%d-s%d ",
		      pedName, fam+1, pedName, fam+1, gen+1, branch+1,
		      ind+1);
	    else
	      fprintf(out, "%s%d %s%d_g%d-b%d-i%d ",
		      pedName, fam+1, pedName, fam+1, gen+1, branch+1,
		      ind - thisBranchNumSpouses + 1);

	    // print parents
	    if (gen == 0 || ind < numFounders) {
	      // first generation or ind >= numNonFounders are founders, so
	      // they have no parents.
	      fprintf(out, "0 0 ");
	    }
	    else {
	      int parBrch[2]; // which branch are the two parents in
	      int parIdx[2];  // index of the Person in the branch?
	      // Is the corresponding parent a spouse of a "primary" individual?
	      bool isSpouse[2] = { false, false };
	      for(int p = 0; p < 2; p++) {
		parBrch[p] = branchParents[gen][branch*2 + p];
		if (parBrch[p] < 0) {
		  assert(p == 1);
		  // founders have negative indexes that start from -1, so we
		  // add 1 to get it to be 0 based and negate to get the index
		  parIdx[p] = -(parBrch[p] + 1);
		  parBrch[1] = parBrch[0];
		  isSpouse[p] = true;
		}
		else {
		  // use "primary" person: immediately after all the spouses
		  int thisBranchNumSpouses;
		  if (branchNumSpouses[gen-1])
		    thisBranchNumSpouses = -branchNumSpouses[gen-1][parBrch[p]];
		  else if (gen == numGen - 1) // no spouse in last generation
		    thisBranchNumSpouses = 0;
		  else                        // one spouse by default
		    thisBranchNumSpouses = 1;
		  parIdx[p] = thisBranchNumSpouses;
		}
	      }

	      int par0sex = prevGenSamps[ parBrch[0] ][ parIdx[0] ].sex;
	      for(int p = 0; p < 2; p++) {
		// print parent 0 first by default, but if parent 0 is female,
		// the following will switch and print parent 1 first
		int printPar = p ^ par0sex;
		if (!isSpouse[ printPar ])
		  // must be the primary person, so i1:
		  fprintf(out, "%s%d_g%d-b%d-i1 ", pedName, fam+1,
			  gen, parBrch[ printPar ]+1);
		else
		  fprintf(out, "%s%d_g%d-b%d-s%d ", pedName, fam+1,
			  gen, parBrch[ printPar ]+1, parIdx[ printPar ]+1);
	      }
	    }

	    // print sex and phenotype:
	    int sex = theSamples[ped][fam][gen][branch][ind].sex;
	    fprintf(out, "%d -9\n", sex+1);
	  }
	}
      }
    }
  }

  fclose(out);
}

// removes the first element from <vec>. This has time that is linear in the
// number of elements. We could use a list<> instead of a vector<> to avoid this
// linear time, but having random access ability on the vectors is very
// convenient, so we'll live with this.
template<typename T>
void pop_front(vector<T> &vec) {
  vec.erase(vec.begin());
}
