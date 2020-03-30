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
#include <set>
#include <random>
#include <sys/time.h>
#include <algorithm>
#include <assert.h>
#include <zlib.h>
#include "cmdlineopts.h"
#include "cointerfere.h"

// TODO! only use sexConstraints array when there are sex-specific maps?
// TODO! make branchNumSpouses positive

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Used to store details about each simulation
struct Parent {
  // The generation and branch number of the given parent
  int gen;
  int branch;
};

struct ParentComp {
  bool operator() (const Parent &lhs, const Parent &rhs) const {
    return (lhs.gen < rhs.gen) ||
				(lhs.gen == rhs.gen && lhs.branch < rhs.branch);
  }
};

////////////////////////////////////////////////////////////////////////////////
// Used to store details about each simulation
struct SimDetails {
  SimDetails(int nFam, int nGen, int **print, int *branches, Parent **parents,
	     int **sexes, int i1FixedSex, int **spouses, char *theName) {
    numFam = nFam;
    numGen = nGen;
    numSampsToPrint = print;
    numBranches = branches;
    branchParents = parents;
    sexConstraints = sexes;
    i1Sex = i1FixedSex;
    branchNumSpouses = spouses;
    name = new char[ strlen(theName) + 1 ];
    if (name == NULL) {
      printf("ERROR: out of memory");
      exit(5);
    }
    strcpy(name, theName);
  }
  int numFam;
  int numGen;
  int **numSampsToPrint;
  int *numBranches;
  Parent **branchParents;
  int **sexConstraints;
  int i1Sex;
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
struct Segment {
  Segment() { }
  Segment(int fhn, int ep) {
    foundHapNum = fhn;
    endPos = ep;
  }
  int foundHapNum, endPos;
};

typedef vector<Segment> Haplotype;
struct Person {
  Person() {
    // by default assume using sex-averaged map: all 0
    sex = 0;
    // by default assume not using fixed COs
    fixedCOidxs[0] = fixedCOidxs[1] = UINT_MAX;
  }
  int sex;
  vector<Haplotype> haps[2]; // haplotype pair for <this>
  unsigned int fixedCOidxs[2];
};

// For the <hapCarriers> structure -- stores the sample id that inherited (or
// is the founder of) a given haplotype
struct InheritRecord {
  InheritRecord() { }
  InheritRecord(unsigned int p, int f, int g, int b, int i, int s, int e) {
    ped = p;
    fam = f;
    gen = g;
    branch = b;
    ind = i;
    startPos = s;
    endPos = e;
  }
  unsigned int ped;
  int fam;
  int gen;
  int branch;
  int ind;
  int startPos;
  int endPos;
};

struct IBDRecord {
  IBDRecord() { assert(false); }
  IBDRecord(int og, int ob, int oi, int ci, int start, int end) {
    otherGen = og;
    otherBranch = ob;
    otherInd = oi;
    chrIdx = ci;
    startPos = start;
    endPos = end;
  }
  int otherGen;
  int otherBranch;
  int otherInd;
  int chrIdx;
  int startPos;
  int endPos;
};

template<typename IO_TYPE>
class FileOrGZ {
  public:
    bool open(const char *filename, const char *mode);
    int getline();
    int printf(const char *format, ...);
    int close();

    static const int INIT_SIZE = 1024 * 50;

    // IO_TYPE is either FILE* or gzFile;
    IO_TYPE fp;

    // for I/O:
    char *buf;
    size_t buf_size;
    size_t buf_len;

  private:
    void alloc_buf();
};

struct FixedCOs {
  static void read(const char *fixedCOfile,
		   vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap);

  static vector<int> & getCOs(int sex, unsigned int idx, unsigned int chrIdx) {
    return theCOs[sex][idx][chrIdx];
  }

  // Array of size 2: one set for males, one for females
  // Outer vector is for individuals, next level for chromosomes, and the last
  // stores a vector of positions of COs
  static vector< vector< vector<int> > > theCOs[2];
};

////////////////////////////////////////////////////////////////////////////////
// Function decls
void readDef(vector<SimDetails> &simDetails, char *defFile);
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				Parent *&thisGenBranchParents, int prevGen,
				int *prevGenSpouseNum = NULL,
				vector<bool> *branchParentsAssigned = NULL);
bool readBranchSpec(int *numBranches, Parent *&thisGenBranchParents,
		    int *thisGenNumSampsToPrint, int curGen,
		    int **sexConstraints, int *&prevGenSpouseNum,
		    vector<bool> &branchParentsAssigned,
		    vector< set<Parent,ParentComp>* > &spouseDependencies,
		    const int i1Sex, const char *delim, char *&saveptr,
		    char *&endptr, int line);
void readParents(int *numBranches, int prevGen, int **sexConstraints,
		 int *&prevGenSpouseNum,
		 vector< set<Parent,ParentComp>* > &spouseDependencies,
		 char *assignBranches, char *assignPar[2], Parent pars[2],
		 char *&fullAssignPar, const int i1Sex, char *&endptr,
		 int line);
void updateSexConstraints(int **sexConstraints, Parent pars[2],
			  int *numBranches,
			  vector< set<Parent,ParentComp>*> &spouseDependencies,
			  int line);
bool intersectNonEmpty(set<Parent,ParentComp> &a, set<Parent,ParentComp> &b);
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps);
void readInterfere(vector<COInterfere> &coIntf, char *interfereFile,
		   vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		   bool &sexSpecificMaps);
int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     bool sexSpecificMaps, vector<COInterfere> &coIntf,
	     vector< vector< vector<InheritRecord> > > &hapCarriers);
void getPersonCounts(int curGen, int numGen, int branch, int **numSampsToPrint,
		     Parent **branchParents, int **branchNumSpouses,
		     int &numFounders, int &numNonFounders);
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       vector<PhysGeneticPos> *curMap,
		       vector<COInterfere> &coIntf, unsigned int chrIdx,
		       vector< vector< vector<InheritRecord> > > &hapCarriers,
		       int ped, int fam, int curGen, int branch, int ind,
		       unsigned int fixedCOidxs[2]);
void copySegs(Haplotype &toGenerate, Person &parent, int &nextSegStart,
	      int switchPos, unsigned int curSegIdx[2], int &curHap,
	      unsigned int chrIdx,
	      vector< vector< vector<InheritRecord> > > &hapCarriers,
	      int ped, int fam, int curGen, int branch, int ind);
int getBranchNumSpouses(SimDetails &pedDetails, int gen, int branch);
template<typename IO_TYPE = FILE *>
bool printSampleId(FILE *out, SimDetails &pedDetails, int fam, int gen,
		   int branch, int ind, bool printAllGens = false,
		   FileOrGZ<IO_TYPE> *gzOut = NULL);
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      char *bpFile);
bool compInheritRec(const InheritRecord &a, const InheritRecord &b);
bool compIBDRecord(const IBDRecord &a, const IBDRecord &b);
void locatePrintIBD(vector<SimDetails> &simDetails,
		    vector< vector< vector<InheritRecord> > > &hapCarriers,
		    vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		    bool sexSpecificMaps, char *ibdFile);
void printIBD(FILE *out, SimDetails &pedDetails, int fam,
	      vector< vector< vector<IBDRecord> > > *theSegs,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps);
void mergeSegments(vector<IBDRecord> &segs);
void printOneIBDSegment(FILE *out, SimDetails &pedDetails, int fam,
		  int gen, int branch, int ind, IBDRecord &seg,
		  int realStart, int realEnd, const char *type,
		  vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		  bool sexSpecificMaps);
void clearTheSegs(SimDetails &pedDetails, 
		  vector< vector< vector<IBDRecord> > > *theSegs);
template<typename IO_TYPE>
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

vector< vector< vector<int> > > FixedCOs::theCOs[2];



////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  int outPrefixLen = strlen(CmdLineOpts::outPrefix);
  char *outFile = new char[ outPrefixLen + 7 + 1 ]; //+7 for .vcf.gz, + 1 for \0
  if (outFile == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }

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
    if (CmdLineOpts::randSeed == random_device().entropy()) {
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
    fprintf(outs[o], "  Input VCF:\t\t%s\n",
	    CmdLineOpts::inVCFfile == NULL ? "[none: no genetic data]" :
					     CmdLineOpts::inVCFfile);
    fprintf(outs[o], "  Output prefix:\t%s\n\n", CmdLineOpts::outPrefix);

    fprintf(outs[o], "  Random seed:\t\t%u\n\n", CmdLineOpts::randSeed);

    if (CmdLineOpts::fixedCOfile)
      fprintf(outs[o], "  Fixed CO file:\t%s\n\n",
	      CmdLineOpts::fixedCOfile);
    else
      fprintf(outs[o], "  Interference file:\t%s\n\n",
	      CmdLineOpts::interfereFile == NULL ? "[none: Poisson model]" :
						    CmdLineOpts::interfereFile);

    if (CmdLineOpts::inVCFfile) {
      // options only relevant when generating data (so when we have input data)
      fprintf(outs[o], "  Genotype error rate:\t%.1le\n",
	      CmdLineOpts::genoErrRate);
      fprintf(outs[o], "  Opposite homozygous error rate:\t%.2lf\n",
	      CmdLineOpts::homErrRate);
      fprintf(outs[o], "  Missingness rate:\t%.1le\n",
	      CmdLineOpts::missRate);
      fprintf(outs[o], "  Pseudo-haploid rate:\t%.1lg\n\n",
	      CmdLineOpts::pseudoHapRate);

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

  if (CmdLineOpts::fixedCOfile) {
    FixedCOs::read(CmdLineOpts::fixedCOfile, geneticMap);
  }

  // The first index is the pedigree number corresponding to the description of
  // the pedigree to be simulated in the def file
  // The second index is the family: we replicate the same pedigree structure
  // some number of times as specified in the def file
  // The third index is the generation number (0-based)
  // The fourth index is the branch of the pedigree
  // The fifth index is the individual number
  Person *****theSamples;

  // Record of which founders and descendants inherited a given haplotype along
  // with the start and end positions.
  // The first index is the haplotype number
  // The second index is the chromosome
  // The third index is the record index
  vector< vector< vector<InheritRecord> > > hapCarriers;

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Simulating haplotype transmissions... ");
    fflush(outs[o]);
  }
  int totalFounderHaps = simulate(simDetails, theSamples, geneticMap,
				  sexSpecificMaps, coIntf, hapCarriers);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  if (CmdLineOpts::printBP) {
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "Printing break points... ");
      fflush(outs[o]);
    }
    sprintf(outFile, "%s.bp", CmdLineOpts::outPrefix);
    printBPs(simDetails, theSamples, geneticMap, /*bpFile=*/ outFile);
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "done.\n");
    }
  }

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Printing IBD segments... ");
    fflush(outs[o]);
  }
  sprintf(outFile, "%s.seg", CmdLineOpts::outPrefix);
  locatePrintIBD(simDetails, hapCarriers, geneticMap, sexSpecificMaps,
		 /*ibdFile=*/ outFile);
  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "done.\n");
  }

  if (CmdLineOpts::inVCFfile) {
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "Reading input VCF meta data... ");
      fflush(outs[o]);
    }
    // note: makeVCF() prints the status for generating the VCF file

    // decide whether to use gz I/O or standard, and call makeVCF() accordingly
    int inVCFlen = strlen(CmdLineOpts::inVCFfile);
    if (strcmp(&CmdLineOpts::inVCFfile[ inVCFlen - 3 ], ".gz") == 0) {
      sprintf(outFile, "%s.vcf.gz", CmdLineOpts::outPrefix);
      makeVCF<gzFile>(simDetails, theSamples, totalFounderHaps,
		      CmdLineOpts::inVCFfile, /*outVCFfile=*/ outFile,
		      geneticMap, outs);
    }
    else {
      sprintf(outFile, "%s.vcf", CmdLineOpts::outPrefix);
      makeVCF<FILE *>(simDetails, theSamples, totalFounderHaps,
		      CmdLineOpts::inVCFfile, /*outVCFfile=*/ outFile,
		      geneticMap, outs);
    }
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "done.\n");
  }

  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "Printing fam file... ");
  fflush(stdout);
  sprintf(outFile, "%s.fam", CmdLineOpts::outPrefix);
  printFam(simDetails, theSamples, /*famFile=*/ outFile);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  fclose(log);

  // NOTE: the memory isn't freed because the OS reclaims it when Ped-sim
  // finishes

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

  // def file gives the number of samples to print; we store this in a 2d array
  // with the first index being generation number and the second index the
  // branch number
  int **curNumSampsToPrint = NULL;
  // Have variable number of branches in each generation
  int *curNumBranches = NULL;
  // Who are the parents of each branch in each generation?
  // Contains <curNumGen> rows, and <2*curNumBranches[gen]> columns on each row.
  // Stores the generation and branch numbers of the two parents.
  // Negative branch values correspond to founders that are stored in the same
  // branch number as the other parent.
  Parent **curBranchParents = NULL;
  // Gives numerical values indicating dependencies of sex assignments for each
  // branch. For example, if the person in branch 1 has children with the
  // individual in branch 2 and 3, branch 2 and 3 must have the same sex.
  int **curSexConstraints = NULL;
  // If all i1 individuals are to have the same sex, the following gives its
  // value. A value of -1 corresponds to random assignment.
  int curI1Sex = -1;
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
  vector< set<Parent,ParentComp>* > spouseDependencies;

  bool warningGiven = false;

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
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
      char *i1SexStr = strtok_r(NULL, delim, &saveptr);
      // Note: leaving i1SexStr out of the next conditional because it can be
      //       NULL or non-NULL
      if (name == NULL || numFamStr == NULL || numGenStr == NULL ||
				      strtok_r(NULL, delim, &saveptr) != NULL) {
	fprintf(stderr, "ERROR: line %d in def: expect four or five fields for pedigree definition:\n",
		line);
	fprintf(stderr, "       def [name] [numFam] [numGen] <sex of i1>\n");
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

      if (i1SexStr == NULL)
	curI1Sex = -1;
      else {
	if (strcmp(i1SexStr, "F") == 0) {
	  curI1Sex = 1;
	}
	else if (strcmp(i1SexStr, "M") == 0) {
	  curI1Sex = 0;
	}
	else {
	  fprintf(stderr, "ERROR: line %d in def: allowed values for sex of i1 field are 'M' and 'F'\n",
		  line);
	  fprintf(stderr, "       got %s\n", i1SexStr);
	  exit(7);
	}
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

      curNumSampsToPrint = new int*[curNumGen];
      curNumBranches = new int[curNumGen];
      curBranchParents = new Parent*[curNumGen];
      curSexConstraints = new int*[curNumGen];
      curBranchNumSpouses = new int*[curNumGen];
      if (curNumSampsToPrint == NULL || curNumBranches == NULL ||
	  curBranchParents == NULL || curSexConstraints == NULL ||
	  curBranchNumSpouses == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
      if (lastReadGen >= 0)
	lastReadGen = -1; // reset

      for(int gen = 0; gen < curNumGen; gen++) {
	// initially
	curNumSampsToPrint[gen] = NULL;
	// set to -1 initially so we know these are unassigned; will update
	// later
	curNumBranches[gen] = -1;
	curBranchParents[gen] = NULL;
	curSexConstraints[gen] = NULL;
	curBranchNumSpouses[gen] = NULL;
      }
      simDetails.emplace_back(curNumFam, curNumGen, curNumSampsToPrint,
			      curNumBranches, curBranchParents,
			      curSexConstraints, curI1Sex,
			      curBranchNumSpouses, name);
      continue;
    }

    ///////////////////////////////////////////////////////////////////////////
    // parse line with information about a generation in the current pedigree

    // is there a current pedigree?
    if (curNumSampsToPrint == NULL) {
      fprintf(stderr, "ERROR: line %d in def: expect four or five fields for pedigree definition:\n",
	      line);
      fprintf(stderr, "       def [name] [numFam] [numGen] <sex of i1>\n");
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
      exit(5);
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
    // Will assign <numSamps> to each branch of this generation below -- first
    // need to know how many branches are in this generation

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
				   curBranchParents[i], /*prevGen=*/i-1);
      // assign default of 0 samples to print
      curNumSampsToPrint[i] = new int[ curNumBranches[i] ];
      for (int b = 0; b < curNumBranches[i]; b++)
	curNumSampsToPrint[i][b] = 0;
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

    curNumSampsToPrint[generation - 1] = new int[thisGenNumBranches];
    for(int b = 0; b < thisGenNumBranches; b++)
      curNumSampsToPrint[generation - 1][b] = numSamps;

    lastReadGen = generation - 1;

    // now read in and assign (if only using the defaults) the branch parents
    // for this generation. Note that in the first generation, all individuals
    // are necessarily founders so there should not be any specification.
    if (generation - 1 > 0) {
      bool warning = readBranchSpec(curNumBranches,
			curBranchParents[generation - 1],
			curNumSampsToPrint[generation - 1],
			/*curGen=*/generation - 1, curSexConstraints,
			/*prevSpouseNum=*/curBranchNumSpouses[generation - 2],
			branchParentsAssigned, spouseDependencies, curI1Sex,
			delim, saveptr, endptr, line);
      warningGiven = warningGiven || warning;
    }
    else if (strtok_r(NULL, delim, &saveptr) != NULL) {
      fprintf(stderr, "ERROR: line %d in def: first generation cannot have parent specifications\n",
	      line);
      exit(8);
    }
  }


  for(auto it = simDetails.begin(); it != simDetails.end(); it++) {
    bool someBranchToPrint = false;
    bool anyNoPrint = false;
    int lastGenNumBranches = it->numBranches[ it->numGen - 1 ];
    for(int b = 0; b < lastGenNumBranches; b++) {
      if (it->numSampsToPrint[ it->numGen - 1 ][b] == 0)
	anyNoPrint = true;
      else
	someBranchToPrint = true;
    }

    if (!someBranchToPrint) {
      fprintf(stderr, "ERROR: request to simulate pedigree \"%s\" with %d generations\n",
	      it->name, it->numGen);
      fprintf(stderr, "       but no request to print any samples from last generation (number %d)\n",
	      it->numGen);
      exit(4);
    }
    else if (anyNoPrint) {
      fprintf(stderr, "Warning: no-print branches in last generation of pedigree %s:\n",
	      it->name);
      fprintf(stderr, "         can omit these branches and possibly reduce number of founders needed\n");
    }
}

  if (simDetails.size() == 0) {
    fprintf(stderr, "ERROR: def file does not contain pedigree definitions;\n");
    fprintf(stderr, "       nothing to simulate\n");
    exit(3);
  }

  fclose(in);

  // don't need the spouse dependency sets anymore:
  for(unsigned int i = 0; i < spouseDependencies.size(); i++)
    if (spouseDependencies[i] != NULL)
      delete spouseDependencies[i];


  if (warningGiven)
    fprintf(stderr, "\n");
}


// Read in genetic map from <mapFile> into <geneticMap>. Also determines whether
// there are male and female maps present and sets <sexSpecificMaps> to true if
// so. If only one map is present, this is assumed to be the sex-averaged map
// and in that case, sets <sexSpecificMap> to false.
void readMap(vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     char *mapFile, bool &sexSpecificMaps) {
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
      exit(2);
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
  if (buffer == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
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

// Given a file <fixedCOfile> containing sets of crossovers labeled paternal and
// maternal, reads in the crossovers for use in simulating
void FixedCOs::read(const char *fixedCOfile,
		  vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap) {
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
	while (chrIdx < (int) geneticMap.size() - 1) {
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

    while (chrIdx == -1 || strcmp(geneticMap[chrIdx].first, chrom) != 0) {
      // new chr:
      theCOs[pm_idx].back().emplace_back();
      chrIdx++;
      if (chrIdx > (int) geneticMap.size()) {
	fprintf(stderr, "ERROR: chromosome %s in %s does not match any chromosome in genetic map\n",
		chrom, fixedCOfile);
      }
    }

    // get position
    posStr = strtok_r(NULL, delim, &saveptr);
    int pos = strtol(posStr, &endptr, 10);
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: column three of %s contains %s, which cannot be converted to an integer\n",
	  fixedCOfile, posStr);
      exit(2);
    }

    if (pos >= geneticMap[chrIdx].second->front().physPos &&
	pos <= geneticMap[chrIdx].second->back().physPos)
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

// Gives the default parent assignment for any branches that have not had
// their parents explicitly specified.
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				Parent *&thisGenBranchParents, int prevGen,
				int *prevGenSpouseNum,
				vector<bool> *branchParentsAssigned) {
  // how many new branches is each previous branch the parent of?
  int multFactor = thisGenNumBranches / prevGenNumBranches;
  if (multFactor == 0)
    multFactor = 1; // for branches that survive, map prev branch i to cur i

  // allocate space to store the parents of each branch
  if (thisGenBranchParents == NULL) {
    thisGenBranchParents = new Parent[2 * thisGenNumBranches];
    if (thisGenBranchParents == NULL) {
      printf("ERROR: out of memory");
      exit(5);
    }
  }

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
      thisGenBranchParents[curBranch*2].gen = prevGen;
      thisGenBranchParents[curBranch*2].branch = prevB;
      if (!spouseNumDefined) {
	if (prevGenSpouseNum) {
	  prevGenSpouseNum[ prevB ]--;
	  spouseNum = prevGenSpouseNum[ prevB ];
	}
	else
	  spouseNum = -1;
	spouseNumDefined = true;
      }
      thisGenBranchParents[curBranch*2 + 1].gen = prevGen;
      thisGenBranchParents[curBranch*2 + 1].branch = spouseNum;
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
    thisGenBranchParents[newB*2].gen = prevGen;
    thisGenBranchParents[newB*2].branch =
				  thisGenBranchParents[newB*2 + 1].branch = -1;
  }
}

// Reads in and performs state changes for branch specifications including
// both parent assignments and no-print directives
// Returns true iff a warning has been printed
bool readBranchSpec(int *numBranches, Parent *&thisGenBranchParents,
		    int *thisGenNumSampsToPrint, int curGen,
		    int **sexConstraints, int *&prevGenSpouseNum,
		    vector<bool> &branchParentsAssigned,
		    vector< set<Parent,ParentComp>* > &spouseDependencies,
		    const int i1Sex, const char *delim, char *&saveptr,
		    char *&endptr, int line) {
  bool warningGiven = false;
  int prevGen = curGen - 1;

  assert(prevGenSpouseNum == NULL);
  prevGenSpouseNum = new int[numBranches[prevGen]];
  if (prevGenSpouseNum == NULL) {
    printf("ERROR: out of memory\n");
    exit(5);
  }
  for(int b = 0; b < numBranches[prevGen]; b++)
    // What number have we assigned through for founder spouses of individuals
    // in the previous generation? Note that founders have an id (in the code
    // for the purposes of <curBranchParents>) that are always negative and
    // that by default we assign one founder spouse to marry one person in
    // each branch in the previous generation (see assignDefaultBranchParents())
    prevGenSpouseNum[b] = 0;

  thisGenBranchParents = new Parent[ 2 * numBranches[curGen] ];
  if (thisGenBranchParents == NULL) {
    printf("ERROR: out of memory\n");
    exit(5);
  }

  sexConstraints[prevGen] = new int[numBranches[prevGen]];
  if (sexConstraints[prevGen] == NULL) {
    printf("ERROR: out of memory\n");
    exit(5);
  }
  for(int i = 0; i < numBranches[prevGen]; i++)
    sexConstraints[prevGen][i] = -1;

  // so far, all the branches in the current generation are assigned default
  // parents; track which branches get explicitly assigned and throw an error
  // if the same branch is assigned more than once
  branchParentsAssigned.clear();
  for(int i = 0; i < numBranches[curGen]; i++)
    branchParentsAssigned.push_back(false);

  while (char *assignToken = strtok_r(NULL, delim, &saveptr)) {
    char *assignBranches = assignToken; // will add '\0' at ':'
    char *fullAssignPar = NULL;
    int i;
    Parent pars[2];

    // split on ':' to get the parent assignments on the right and the
    // branches on the left
    for(i = 0; assignToken[i] != ':' && assignToken[i] != 'n' &&
	       assignToken[i] != '\0'; i++);
    // ':' gives parents (after the ':') and n means the preceeding branches
    // should not have their members printed
    bool parentAssign = false;
    bool noPrint = false;
    if (assignToken[i] == ':')
      parentAssign = true;
    else if (assignToken[i] == 'n')
      noPrint = true;
    else {
      fprintf(stderr, "ERROR: line %d in def: improperly formatted parent assignment or no-print\n",
	      line);
      fprintf(stderr, "       field %s\n", assignToken);
      exit(8);
    }
    assignToken[i] = '\0';

    if (noPrint) {
      // expect a space after the 'n': check this
      if (assignToken[i+1] != '\0') {
	assignToken[i] = 'n';
	fprintf(stderr, "ERROR: line %d in def: improperly formatted no-print field \"%s\":\n",
		line, assignToken);
	fprintf(stderr, "       no-print character 'n' should be followed by white space\n");
	exit(8);
      }
    }
    if (parentAssign) {
      // Get the one or two parents
      char *assignPar[2];
      assignPar[0] = &(assignToken[i+1]); // will add '\0' at '_' if present
      assignPar[1] = NULL; // initially; updated just below

      readParents(numBranches, prevGen, sexConstraints, prevGenSpouseNum,
		  spouseDependencies, assignBranches, assignPar, pars,
		  fullAssignPar, i1Sex, endptr, line);
    }


    // one should be true, not both:
    assert((parentAssign || noPrint) && !(parentAssign && noPrint));

    // process the branches
    // if <parentAssign>, these will be assigned <pars> as parent OR
    // if <noPrint>, these will not be printed
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
      else if (assignBranches[i] == ',' || assignBranches[i] == '\0') {
	if (assignBranches[i] == '\0') {
	  done = true;
	  if (i == 0)
	    break; // early termination, will get error below
	}
	else
	  assignBranches[i] = '\0';

	int curBranch = strtol(assignBranches, &endptr, 10) - 1; // 0 indexed
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: line %d in def: unable to parse branch %s to ",
		  line, assignBranches);
	  if (parentAssign)
	    fprintf(stderr, "assign parent %s to\n", fullAssignPar);
	  else // no print
	    fprintf(stderr, "set as no-print\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}

	if (startBranch) {
	  int rangeEnd = curBranch;
	  int rangeStart = strtol(startBranch, &endptr, 10) - 1; // 0 indexed
	  if (errno != 0 || *endptr != '\0') {
	    fprintf(stderr, "ERROR: line %d in def: unable to parse branch %s to ",
		    line, startBranch);
	    if (parentAssign)
	      fprintf(stderr, "assign parent %s to\n", fullAssignPar);
	    else // no print
	      fprintf(stderr, "set as no-print\n");
	    if (errno != 0)
	      perror("strtol");
	    exit(2);
	  }
	  startBranch = NULL; // parsed: reset this variable

	  if (rangeStart >= rangeEnd) {
	    fprintf(stderr, "ERROR: line %d in def: non-increasing branch range %d-%d to\n",
		    line, rangeStart, rangeEnd);
	    if (parentAssign)
	      fprintf(stderr, "       assign parent %s to\n", fullAssignPar);
	    else // no print
	      fprintf(stderr, "       set as no-print\n");
	    exit(8);
	  }

	  for(int branch = rangeStart; branch <= rangeEnd; branch++) {
	    if (parentAssign) {
	      if (branchParentsAssigned[branch]) {
		fprintf(stderr, "ERROR: line %d in def: parents of branch number %d assigned multiple times\n",
			line, branch+1);
		exit(8);
	      }
	      branchParentsAssigned[branch] = true;
	      for(int p = 0; p < 2; p++)
		thisGenBranchParents[branch*2 + p] = pars[p];
	    }
	    else { // no print
	      // print 0 samples for <branch>
	      if (thisGenNumSampsToPrint[branch] > 1) {
		fprintf(stderr, "Warning: line %d in def: generation %d would print %d individuals, now set to 0\n",
			line, curGen + 1, thisGenNumSampsToPrint[branch]);
		warningGiven = true;
	      }
	      else if (thisGenNumSampsToPrint[branch] == 0) {
		fprintf(stderr, "Warning: line %d in def: generation %d branch %d, no-print is redundant\n",
			line, curGen + 1, branch + 1);
		warningGiven = true;
	      }
	      thisGenNumSampsToPrint[branch] = 0;
	    }
	  }
	}
	else {
	  if (parentAssign) {
	    if (branchParentsAssigned[curBranch]) {
	      fprintf(stderr, "ERROR: line %d in def: parents of branch number %d assigned multiple times\n",
		      line, curBranch + 1);
	      exit(8);
	    }
	    branchParentsAssigned[curBranch] = true;
	    for(int p = 0; p < 2; p++)
	      thisGenBranchParents[curBranch*2 + p] = pars[p];
	  }
	  else { // no print
	    // print 0 samples for <curBranch>
	    if (thisGenNumSampsToPrint[curBranch] > 1) {
	      fprintf(stderr, "Warning: line %d in def: generation %d would print %d individuals, now set to 0\n",
		      line, curGen + 1, thisGenNumSampsToPrint[curBranch]);
	      warningGiven = true;
	    }
	    else if (thisGenNumSampsToPrint[curBranch] == 0) {
	      fprintf(stderr, "Warning: line %d in def: generation %d branch %d, no-print is redundant\n",
		      line, curGen + 1, curBranch + 1);
	      warningGiven = true;
	    }
	    thisGenNumSampsToPrint[curBranch] = 0;
	  }
	}

	assignBranches = &(assignBranches[i+1]); // go through next in loop
      }
    }

    if (startBranch != NULL) {
      fprintf(stderr, "ERROR: line %d in def: range of branches ", line);
      if (parentAssign)
	fprintf(stderr, "to assign parents ");
      else // no print
	fprintf(stderr, "set as no-print ");
      fprintf(stderr, "does not terminate\n");
      exit(8);
    }
  }

  // TODO: As a space optimization, could compress the spouseDependencies
  //       vector more. When combining two sets, the indexes of one of the set
  //       are no longer used, but the size of the spouseDependencies list
  //       still accounts for them.
  assert(spouseDependencies.size() % 2 == 0);

  assignDefaultBranchParents(numBranches[prevGen], numBranches[curGen],
			     thisGenBranchParents, prevGen, prevGenSpouseNum,
			     &branchParentsAssigned);

  return warningGiven;
}

// In the branch specifications, read the parent assignments
void readParents(int *numBranches, int prevGen, int **sexConstraints,
		 int *&prevGenSpouseNum,
		 vector< set<Parent,ParentComp>* > &spouseDependencies,
		 char *assignBranches, char *assignPar[2], Parent pars[2],
		 char *&fullAssignPar, const int i1Sex, char *&endptr,
		 int line) {
  int i;

  // Find the second parent if present
  for(i = 0; assignPar[0][i] != '_' && assignPar[0][i] != '\0'; i++);
  if (assignPar[0][i] == '_') {
    assignPar[0][i] = '\0';
    assignPar[1] = &(assignPar[0][i+1]);
  }

  for(int p = 0; p < 2; p++) {
    pars[p].branch = -1;
    pars[p].gen = prevGen;
  }
  for(int p = 0; p < 2 && assignPar[p] != NULL && assignPar[p][0] != '\0'; p++){
    // Check for generation number
    char *genNumStr = NULL;
    for(i = 0; assignPar[p][i] != '^' && assignPar[p][i] != '\0'; i++);
    if (assignPar[p][i] == '^') {
      // Have a generation number
      if (p == 0) {
	fprintf(stderr, "ERROR: line %d in def: parent assignment for branches %s gives generation\n",
		line, assignBranches);
	fprintf(stderr, "       number for the first parent, but this is only allowed for the second\n");
	fprintf(stderr, "       parent; for example, 2:1_3^1 has branch 1 from previous generation\n");
	fprintf(stderr, "       married to branch 3 from generation 1\n");
	exit(3);
      }
      assignPar[p][i] = '\0';
      genNumStr = &(assignPar[p][i+1]);
      pars[p].gen = strtol(genNumStr, &endptr, 10) - 1; // 0 indexed => -1
      if (errno != 0 || *endptr != '\0') {
	fprintf(stderr, "ERROR: line %d in def: unable to parse parent assignment for branches %s\n",
		line, assignBranches);
	fprintf(stderr, "       malformed generation number string for second parent: %s\n",
		genNumStr);
	if (errno != 0)
	  perror("strtol");
	exit(5);
      }
      if (pars[p].gen > prevGen) {
	fprintf(stderr, "ERROR: line %d in def: unable to parse parent assignment for branches %s\n",
		line, assignBranches);
	fprintf(stderr, "       generation number %s for second parent is after previous generation\n",
		genNumStr);
	exit(-7);
      }
      else if (pars[p].gen < 0) {
	fprintf(stderr, "ERROR: line %d in def: unable to parse parent assignment for branches %s\n",
		line, assignBranches);
	fprintf(stderr, "       generation number %s for second parent is before first generation\n",
		genNumStr);
	exit(5);
      }
    }

    pars[p].branch = strtol(assignPar[p], &endptr, 10) - 1; // 0 indexed => -1
    if (errno != 0 || *endptr != '\0') {
      fprintf(stderr, "ERROR: line %d in def: unable to parse parent assignment for branches %s\n",
	      line, assignBranches);
      if (errno != 0)
	perror("strtol");
      exit(2);
    }
    if (pars[p].branch < 0) {
      fprintf(stderr, "ERROR: line %d in def: parent assignments must be of positive branch numbers\n",
	      line);
      exit(8);
    }
    else if (pars[p].branch >= numBranches[ pars[p].gen ]) {
      fprintf(stderr, "ERROR: line %d in def: parent branch number %d is more than the number of\n",
	      line, pars[p].branch+1);
      fprintf(stderr, "       branches (%d) in generation %d\n",
	      numBranches[ pars[p].gen ], pars[p].gen+1);
      exit(8);
    }
    // so that we can print the parent assignment in case of errors below
    if (genNumStr != NULL)
      genNumStr[-1] = '^';
  }
  if (pars[0].branch == -1) {
    // new founder
    assert(pars[1].branch == -1);
  }
  else if (pars[1].branch == -1) {
    // Have not yet assigned the numerical id of parent 1. Because the def
    // file doesn't specify this, it is a founder, and one that hasn't been
    // assigned before. As such, we'll get a unique number associated with a
    // spouse of pars[0].branch. Negative values correspond to founders, so
    // we decrement <prevGenSpouseNum>. (It is initialized to 0 above)
    prevGenSpouseNum[ pars[0].branch ]--;
    pars[1].branch = prevGenSpouseNum[ pars[0].branch ];
  }
  else {
    if (pars[0].branch == pars[1].branch && pars[0].gen == pars[1].gen) {
      fprintf(stderr, "ERROR: line %d in def: cannot have both parents be from same branch\n",
	      line);
      exit(8);
    }
    if (i1Sex >= 0) {
      fprintf(stderr, "ERROR: line %d in def: cannot have fixed sex for i1 samples and marriages\n",
	      line);
      fprintf(stderr, "       between branches -- i1's will have the same sex and cannot reproduce\n");
      exit(9);
    }
    updateSexConstraints(sexConstraints, pars, numBranches,
			 spouseDependencies, line);
  }

  // so that we can print the parent assignment in case of errors below
  if (assignPar[1] != NULL)
    assignPar[1][-1] = '_';
  fullAssignPar = assignPar[0];
}

// Given the branch indexes of two parents, adds constraints and error checks
// to ensure that this couple does not violate the requirement that parents
// must have opposite sex.
void updateSexConstraints(int **sexConstraints, Parent pars[2],
			  int *numBranches,
			  vector< set<Parent,ParentComp>* > &spouseDependencies,
			  int line) {
  // we check these things in the caller, but just to be sure:
  for(int p = 0; p < 2; p++) {
    assert(pars[p].branch >= 0);
    assert(pars[p].branch < numBranches[ pars[p].gen ]);
  }

  if (sexConstraints[ pars[0].gen ][ pars[0].branch ] == -1 &&
      sexConstraints[ pars[1].gen ][ pars[1].branch ] == -1) {
    // neither is a member of a spouse set: create and add to
    // <spouseDependencies>
    set<Parent,ParentComp> *sets[2];
    for(int p = 0; p < 2; p++) {
      sets[p] = new set<Parent,ParentComp>();
      if (sets[p] == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
      sets[p]->insert( pars[p] );
      // which index in spouseDependencies is the set corresponding to this
      // parent stored in? As we're about to add it, just below this, the
      // current size will be the index
      sexConstraints[pars[p].gen][pars[p].branch] = spouseDependencies.size();
      spouseDependencies.push_back( sets[p] );
    }
  }
  else if (sexConstraints[ pars[0].gen ][ pars[0].branch ] == -1 ||
	   sexConstraints[ pars[1].gen ][ pars[1].branch ] == -1) {
    // one spouse is a member of a spouse set and the other is not: add the
    // other parent to the opposite spouse set, error check, and update state
    int assignedPar = -1;
    if (sexConstraints[ pars[0].gen ][ pars[0].branch ] >= 0)
      assignedPar = 0;
    else
      assignedPar = 1;

    int otherPar = assignedPar ^ 1;
    int assignedSetIdx = sexConstraints[ pars[assignedPar].gen ]
						  [ pars[assignedPar].branch ];

    // since <otherPar> isn't in a spouse set yet, it definitely shouldn't be in
    // the same spouse set as <assignedPar>
    assert((unsigned int) assignedSetIdx < spouseDependencies.size());
    assert(spouseDependencies[ assignedSetIdx ] != NULL);
    assert(spouseDependencies[ assignedSetIdx ]->find( pars[otherPar] ) ==
				    spouseDependencies[assignedSetIdx]->end());

    // NOTE: The following is a trick that relies on the fact that sets are
    // stored as sequential pairs in <spouseDependencies>. Ex: indexes 0 and 1
    // are pairs of spouses that are dependent upon each other (those in index
    // 0 must have the same sex and must be opposite those in index 1).  To get
    // an even number that is 1 minus an odd value or the odd number that is 1
    // plus an even value, it suffices to flip the lowest order bit:
    int otherSetIdx = assignedSetIdx ^ 1;
    // add <pars[otherPar]> to the set containing spouses of <pars[assignedPar]>
    spouseDependencies[ otherSetIdx ]->insert( pars[otherPar] );
    sexConstraints[ pars[otherPar].gen ][ pars[otherPar].branch ] = otherSetIdx;
  }
  else {
    assert(sexConstraints[ pars[0].gen ][ pars[0].branch ] >= 0 &&
			  sexConstraints[ pars[1].gen ][ pars[1].branch ] >= 0);
    // both are members of spouse sets

    if (sexConstraints[ pars[0].gen ][ pars[0].branch ] / 2 ==
			  sexConstraints[ pars[1].gen ][ pars[1].branch ] / 2) {
      // spouse sets from the same pair of sets
      if (sexConstraints[ pars[0].gen ][ pars[0].branch ] ==
			      sexConstraints[ pars[1].gen ][ pars[1].branch ]) {
	fprintf(stderr, "ERROR: line %d in def: assigning branch %d from generation %d and branch %d from\n",
		line, pars[0].branch+1, pars[0].gen+1, pars[1].branch+1);
	fprintf(stderr, "       generation %d as parents is impossible due to other parent assignments:\n",
		pars[1].gen+1);
	fprintf(stderr, "       they necessarily have same sex\n");
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
	setIdxes[p][0] = sexConstraints[ pars[p].gen ][ pars[p].branch ];
	// see comment denoted with NOTE a little ways above for why this works:
	setIdxes[p][1] = setIdxes[p][0] ^ 1;
      }

      for(int i = 0; i < 2; i++)
	for(int j = 0; j < 2; j++)
	  assert(spouseDependencies[ setIdxes[i][j] ] != NULL);

      // This is legacy code, but I'll keep it anyway:
      // ensure the values are different pointers (old code used to update
      // the pointers directly instead of setting one of the two merged sets
      // to NULL)
      assert(spouseDependencies[ setIdxes[0][0] ] !=
					spouseDependencies[ setIdxes[1][1] ]);
      assert(spouseDependencies[ setIdxes[0][1] ] !=
					  spouseDependencies[ setIdxes[1][0] ]);

      // take the union such that spouses of each are contained in both
      spouseDependencies[ setIdxes[0][0] ]->insert(
				  spouseDependencies[ setIdxes[1][1] ]->begin(),
				  spouseDependencies[ setIdxes[1][1] ]->end());
      spouseDependencies[ setIdxes[0][1] ]->insert(
				  spouseDependencies[ setIdxes[1][0] ]->begin(),
				  spouseDependencies[ setIdxes[1][0] ]->end());

      // intersection of these should be empty, otherwise there's
      // inconsistencies and individuals that are meant to be different sexes
      // simultaneously
      if (intersectNonEmpty(*spouseDependencies[ setIdxes[0][0] ],
			    *spouseDependencies[ setIdxes[0][1] ])) {
	fprintf(stderr, "ERROR: line %d in def: assigning branch %d from generation %d and branch %d from\n",
		line, pars[0].branch+1, pars[0].gen+1, pars[1].branch+1);
	fprintf(stderr, "       generation %d as parents is impossible due to other parent assignments:\n",
		pars[1].gen+1);
	fprintf(stderr, "       they necessarily have same sex\n");
	exit(4);
      }

      // replace set values stored in sexConstraints to the newly merged set
      // indexes
      for(int p = 0; p < 2; p++) {
	// update all set indexes for samples in
	// <spouseDependencies[ setIdxes[1][p] ]> to setIdxes[0][ p^1 ];
	set<Parent,ParentComp> *toUpdate = spouseDependencies[ setIdxes[1][p] ];
	for(auto it = toUpdate->begin(); it != toUpdate->end(); it++)
	  sexConstraints[ it->gen ][ it->branch ] = setIdxes[0][p^1];

	delete toUpdate;
	spouseDependencies[ setIdxes[1][p] ] = NULL;
      }
    }
  }
}

// Returns true if the intersection of <a> and <b> is non-empty
bool intersectNonEmpty(set<Parent,ParentComp> &a, set<Parent,ParentComp> &b) {
  // Similar to the std::set_intersection code
  auto a_it = a.begin();
  auto b_it = b.begin();
  while(a_it != a.end() && b_it != b.end()) {
    if (ParentComp()(*a_it, *b_it))
      a_it++;
    else if (ParentComp()(*b_it, *a_it))
      b_it++;
    else {
      // have *a_it == *b_it, so intersection is non-empty:
      return true;
    }
  }
  return false;
}

// Simulate data for each specified pedigree type for the number of requested
// families. Returns the number of founder haplotypes used to produce these
// simulated samples.
int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     bool sexSpecificMaps, vector<COInterfere> &coIntf,
	     vector< vector< vector<InheritRecord> > > &hapCarriers) {
  // Note: throughout we use 0-based values for generations though the input
  // def file is 1-based

  int totalFounderHaps = 0;
  // Stores the random assignments of sex for the parents in a given generation
  // This relates to the constraints on sexes of individuals. Individuals with
  // the same index -- an index into the <sexAssignments> vector -- will be
  // assigned the same sex. And any pair of indexes i and i^1 will have
  // opposite sex.
  vector<int> sexAssignments;

  // For randomly assigning a set of fixed (read in, probably from real data)
  // COs to each person. We'll shuffle a list of indexes of fixed crossovers and
  // assign one to each simulated non-founder (for the COs they inherit) in turn
  vector<unsigned int> randFixedCOs[2];
  int curFixedCOidx = -1; // initially: -1 indicates we're not using fixed COs
  if (!FixedCOs::theCOs[0].empty()) {
    curFixedCOidx = 0; // will assign fixed COs to non-founders

    for(int s = 0; s < 2; s++) {
      unsigned int numFixedCOs = FixedCOs::theCOs[s].size();
      for(unsigned int i = 0; i < numFixedCOs; i++) {
	randFixedCOs[s].push_back(i);
      }
      shuffle(randFixedCOs[s].begin(), randFixedCOs[s].end(), randomGen);
    }
  }

  theSamples = new Person****[simDetails.size()];
  if (theSamples == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  for(unsigned int ped = 0; ped < simDetails.size(); ped++) { // for each ped
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **sexConstraints = simDetails[ped].sexConstraints;
    int i1Sex = simDetails[ped].i1Sex;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    ////////////////////////////////////////////////////////////////////////////
    // Allocate space and make Person objects for all those we will simulate,
    // assigning sex if <sexSpecificMaps> is true
    theSamples[ped] = new Person***[numFam];
    if (theSamples[ped] == NULL) {
      printf("ERROR: out of memory");
      exit(5);
    }
    for (int fam = 0; fam < numFam; fam++) {

      // ready to make sex assignments for this family
      sexAssignments.clear();

      theSamples[ped][fam] = new Person**[numGen];
      if (theSamples[ped][fam] == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
      for(int curGen = 0; curGen < numGen; curGen++) {

	theSamples[ped][fam][curGen] = new Person*[ numBranches[curGen] ];
	if (theSamples[ped][fam][curGen] == NULL) {
	  printf("ERROR: out of memory");
	  exit(5);
	}

	// allocate Persons for each branch of <curGen>, assign their sex,
	// and do the simulation for these allocated individuals
	for(int branch = 0; branch < numBranches[curGen]; branch++) {
	  // Determine how many founders and non-founders we need data for in
	  // <branch>:
	  int numFounders, numNonFounders;
	  getPersonCounts(curGen, numGen, branch, numSampsToPrint,
			  branchParents, branchNumSpouses, numFounders,
			  numNonFounders);

	  // Will use convention that all founder spouses are stored as indexes
	  // 0 through <branchNumSpouses>, the "primary" samples is next and
	  // all other non-founders follow.
	  // Note that when the branch is new and has no parents, all
	  // individuals are founders.
	  int numPersons = numFounders + numNonFounders;
	  theSamples[ped][fam][curGen][branch] = new Person[numPersons];
	  if (theSamples[ped][fam][curGen][branch] == NULL) {
	    printf("ERROR: out of memory");
	    exit(5);
	  }

	  if (sexSpecificMaps) {
	    int branchAssign;
	    if (i1Sex >= 0)
	      branchAssign = i1Sex;
	    else if (sexConstraints[curGen] == NULL ||
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

	  Segment trivialSeg;

	  // for each chromosome:
	  unsigned int numChrs = geneticMap.size();
	  for(unsigned int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
	    vector<PhysGeneticPos> *curMap = geneticMap[chrIdx].second;

	    // Chromosome start/ends
	    int chrStart = curMap->front().physPos;
	    int chrEnd = curMap->back().physPos;

	    // Trivial haplotypes for founders in current generation:
	    // no crossovers in founders
	    trivialSeg.endPos = chrEnd;

	    // Simulate the founders for this chromosome:
	    if (curGen != numGen - 1) { // no founders in the last generation
	      for(int ind = 0; ind < numFounders; ind++) {
		for(int h = 0; h < 2; h++) { // 2 founder haplotypes per founder
		  int foundHapNum;
		  if (chrIdx == 0) {
		    foundHapNum = totalFounderHaps++;

		    hapCarriers.emplace_back();
		    hapCarriers[ foundHapNum ].reserve( numChrs );
		    hapCarriers[ foundHapNum ].resize( numChrs );
		  }
		  else
		    // want the same founder on all chromosomes, so access
		    // the haplotype number assigned to the previous
		    // chromosome for this person:
		    foundHapNum =
			theSamples[ped][fam][curGen][branch][ind].haps[h].
						      back().back().foundHapNum;

		  trivialSeg.foundHapNum = foundHapNum;

		  // the following copies <trivialSeg>, so we can reuse it
		  theSamples[ped][fam][curGen][branch][ind].haps[h].
								 emplace_back();
		  theSamples[ped][fam][curGen][branch][ind].haps[h].back().
							  push_back(trivialSeg);

		  // print this branch?
		  if (numSampsToPrint[curGen][branch] > 0) {
		    hapCarriers[ foundHapNum ][ chrIdx ].emplace_back(
			ped, fam, curGen, branch, ind, chrStart, chrEnd);
		  }

		}
	      }
	    }

	    if (numNonFounders == 0) {
	      assert(curGen == 0 ||
				  branchParents[curGen][branch*2].branch == -1);
	      continue; // no non-founders in first generation
	    }

	    assert(curGen > 0 && branchParents[curGen][branch*2].branch >= 0);

	    // Now simulate the non-founders in <branch> of <curGen>
	    //
	    // First figure out who the parents are:
	    Parent pars[2];
	    int parIdx[2];  // index of the Person in the branch
	    for(int p = 0; p < 2; p++) {
	      pars[p] = branchParents[curGen][branch*2 + p];
	      if (pars[p].branch < 0) {
		assert(p == 1 && pars[p].gen == curGen - 1);
		// founders have negative indexes that start from -1, so we
		// add 1 to get it to be 0 based and negate to get the index
		parIdx[p] = -(pars[p].branch + 1);
		pars[1].branch = pars[0].branch;
	      }
	      else {
		// use "primary" person: immediately after all the spouses
		int thisBranchNumSpouses;
		if (branchNumSpouses[curGen-1])
		  thisBranchNumSpouses =
		      -branchNumSpouses[ pars[p].gen ][ pars[p].branch ];
		else if (curGen == numGen - 1) // no spouses in last generation
		  thisBranchNumSpouses = 0;
		else                           // one spouse by default
		  thisBranchNumSpouses = 1;
		parIdx[p] = thisBranchNumSpouses;
	      }
	    }

	    Person ***curFamSamps = theSamples[ped][fam];
	    // the non-founders are stored just after the founders
	    for(int ind = numFounders; ind < numPersons; ind++) {
	      // If we're using sex-specific maps, the two parents' sexes should
	      // differ
	      if (sexSpecificMaps) {
		assert(curFamSamps[ pars[0].gen ][ pars[0].branch ][ parIdx[0] ].sex !=
		       curFamSamps[ pars[1].gen ][ pars[1].branch ][ parIdx[1] ].sex);
	      }

	      for(int p = 0; p < 2; p++) { // meioses from each parent index <p>
		Person &theParent = curFamSamps[ pars[p].gen ][ pars[p].branch ][ parIdx[p] ];

		int hapIdx; // haplotype index for the simulated sample
		if (sexSpecificMaps)
		  // match the sex of the parent if using sex-specific maps
		  hapIdx = theParent.sex;
		else
		  hapIdx = p;

		// Make space for this haplotype in the current sample:
		curFamSamps[curGen][branch][ind].haps[hapIdx].emplace_back();
		Person &thePerson = curFamSamps[curGen][branch][ind];
		if (p == 0 && chrIdx == 0 && curFixedCOidx >= 0) {
		  // assign fixed COs for <thePerson>:
		  for(int s = 0; s < 2; s++)
		    thePerson.fixedCOidxs[s] = randFixedCOs[s][curFixedCOidx];
		  curFixedCOidx++; // next person's indexes
		}
		Haplotype &toGen = thePerson.haps[hapIdx].back();
		generateHaplotype(toGen, theParent, curMap, coIntf, chrIdx,
				  hapCarriers,
				  (numSampsToPrint[curGen][branch] > 0) ? ped
									: -1,
				  fam, curGen, branch, ind,
				  thePerson.fixedCOidxs);
	      } // <parIdx> (simulate each transmitted haplotype for <ind>)
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
void getPersonCounts(int curGen, int numGen, int branch, int **numSampsToPrint,
		     Parent **branchParents, int **branchNumSpouses,
		     int &numFounders, int &numNonFounders) {
  if (curGen > 0 &&
      branchParents[curGen][branch*2].branch >= 0) { // have parent(s)?
    numNonFounders = numSampsToPrint[curGen][branch];
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
    assert(curGen == 0 || branchParents[curGen][branch*2 + 1].branch == -1);

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
		       vector<COInterfere> &coIntf, unsigned int chrIdx,
		       vector< vector< vector<InheritRecord> > > &hapCarriers,
		       int ped, int fam, int curGen, int branch, int ind,
		       unsigned int fixedCOidxs[2]) {
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

  // first segment starts at first valid map position
  int nextSegStart = curMap->front().physPos;

  if (fixedCOidxs[0] == UINT_MAX) { // no fixed crossovers -- simulate:
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

    // initially assume we'll recombine between first and second positions with
    // map info
    int switchIdx = 0;

    int mapNumPos = curMap->size();
    for(auto it = coLocations.begin(); it != coLocations.end(); it++) {
      // Multiply by 100 to get cM:
      double cMPosNextCO = firstcMPos + (*it * 100);

      int left = 0, right = curMap->size() - 1;
      while (true) {
	if (right - left == 1) {
	  switchIdx = left; // want <switchIdx> <= than <cMPosNextCO>
	  break;
	}
	int mid = (left + right) / 2;
	if ((*curMap)[mid].mapPos[mapIdx] < cMPosNextCO)
	  left = mid;
	else if ((*curMap)[mid].mapPos[mapIdx] > cMPosNextCO)
	  right = mid;
	else {
	  // equal: exact map position
	  switchIdx = mid;
	  break;
	}
      }
      if (switchIdx == mapNumPos - 1)
	break; // let code below this while loop insert the final segments

      // get physical position using linear interpolation:
      double frac = (cMPosNextCO - (*curMap)[switchIdx].mapPos[mapIdx]) /
	      ((*curMap)[switchIdx+1].mapPos[mapIdx] -
					   (*curMap)[switchIdx].mapPos[mapIdx]);
      assert(frac >= 0.0 && frac <= 1.0);
      int switchPos = (*curMap)[switchIdx].physPos +
	frac * ((*curMap)[switchIdx+1].physPos - (*curMap)[switchIdx].physPos);

      // copy Segments from <curHap>
      copySegs(toGenerate, parent, nextSegStart, switchPos, curSegIdx, curHap,
	       chrIdx, hapCarriers, ped, fam, curGen, branch, ind);
    }
  }
  else {
    vector<int> &theCOs = FixedCOs::getCOs(parent.sex, fixedCOidxs[parent.sex],
					   chrIdx);

    for(auto it = theCOs.begin(); it != theCOs.end(); it++) {
      // copy Segments from <curHap>
      copySegs(toGenerate, parent, nextSegStart, /*switchPos=*/ *it, curSegIdx,
	       curHap, chrIdx, hapCarriers, ped, fam, curGen, branch, ind);
    }
  }
  

  // copy through to the end of the chromosome:
  for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
    Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
    toGenerate.push_back(seg);

    if (ped >= 0) {
      hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	  ped, fam, curGen, branch, ind, nextSegStart, seg.endPos);
      nextSegStart = seg.endPos + 1;
    }
  }
}

// Copies the <Segment>s between <nextSegStart> and <switchPos> from <parent>'s
// <curHap> haplotype to <toGenerate>.
void copySegs(Haplotype &toGenerate, Person &parent, int &nextSegStart,
	      int switchPos, unsigned int curSegIdx[2], int &curHap,
	      unsigned int chrIdx,
	      vector< vector< vector<InheritRecord> > > &hapCarriers,
	      int ped, int fam, int curGen, int branch, int ind) {
  for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
    Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
    if (seg.endPos >= switchPos) {
      // last segment to copy, and we will break it at <switchPos>
      toGenerate.emplace_back(seg.foundHapNum, switchPos);
      if (seg.endPos == switchPos)
	curSegIdx[curHap]++;

      if (ped >= 0) {
	hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	    ped, fam, curGen, branch, ind, nextSegStart, switchPos);
	nextSegStart = switchPos + 1;
      }
      break; // done copying
    }
    else {
      toGenerate.push_back(seg);

      if (ped >= 0) {
	hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	    ped, fam, curGen, branch, ind, nextSegStart, seg.endPos);
	nextSegStart = seg.endPos + 1;
      }
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

// Returns the number of spouses a given generation <gen> and <branch> has
int getBranchNumSpouses(SimDetails &pedDetails, int gen, int branch) {
  if (pedDetails.branchNumSpouses[gen])
    return -pedDetails.branchNumSpouses[gen][branch];
  else if (gen == pedDetails.numGen - 1) // no spouses in last generation
    return 0;
  else // one spouse by default
    return 1;
}

// Prints the sample id of the given sample to <out>.
// Returns true if the sample is a founder, false otherwise.
template<class IO_TYPE>
bool printSampleId(FILE *out, SimDetails &pedDetails, int fam, int gen,
		   int branch, int ind, bool printAllGens,
		   FileOrGZ<IO_TYPE> *gzOut) {
  int thisBranchNumSpouses = getBranchNumSpouses(pedDetails, gen, branch);
  bool shouldPrint = pedDetails.numSampsToPrint[gen][branch] >0 || printAllGens;

  if (ind < thisBranchNumSpouses) {
    if (shouldPrint) {
      if (!gzOut)
	fprintf(out, "%s%d_g%d-b%d-s%d", pedDetails.name, fam+1, gen+1,
		branch+1, ind+1);
      else
	gzOut->printf("%s%d_g%d-b%d-s%d", pedDetails.name, fam+1, gen+1,
		      branch+1, ind+1);
    }
    return true; // is a founder
  }
  else {
    if (shouldPrint) {
      if (!gzOut)
	fprintf(out, "%s%d_g%d-b%d-i%d", pedDetails.name, fam+1, gen+1,
		branch+1, ind - thisBranchNumSpouses + 1);
      else
	gzOut->printf("%s%d_g%d-b%d-i%d", pedDetails.name, fam+1, gen+1,
		      branch+1, ind - thisBranchNumSpouses + 1);
    }
    if (gen == 0 || pedDetails.branchParents[gen][branch*2].branch < 0) {
      assert(ind - thisBranchNumSpouses == 0);
      return true; // is a founder
    }
    return false;
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
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    for(int fam = 0; fam < numFam; fam++) {
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  if (numSampsToPrint[gen][branch] > 0) {
	    int numNonFounders, numFounders;
	    getPersonCounts(gen, numGen, branch, numSampsToPrint,
			    branchParents, branchNumSpouses, numFounders,
			    numNonFounders);
	    int numPersons = numNonFounders + numFounders;

	    for(int ind = 0; ind < numPersons; ind++) {
	      for(int h = 0; h < 2; h++) {
		int sex = theSamples[ped][fam][gen][branch][ind].sex;
		printSampleId(out, simDetails[ped], fam, gen, branch, ind);
		fprintf(out, " s%d h%d", sex, h);

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

bool compInheritRec(const InheritRecord &a, const InheritRecord &b) {
  return a.startPos < b.startPos;
}

bool compIBDRecord(const IBDRecord &a, const IBDRecord &b) {
  return (a.otherGen < b.otherGen) ||
	 (a.otherGen == b.otherGen && a.otherBranch < b.otherBranch) ||
	 (a.otherGen == b.otherGen && a.otherBranch == b.otherBranch &&
	 a.otherInd < b.otherInd) ||
	 (a.otherGen == b.otherGen && a.otherBranch == b.otherBranch &&
	 a.otherInd == b.otherInd && a.chrIdx < b.chrIdx) ||
	 (a.otherGen == b.otherGen && a.otherBranch == b.otherBranch &&
	 a.otherInd == b.otherInd && a.chrIdx == b.chrIdx &&
	 a.startPos < b.startPos);
}

// Locates and prints IBD segments using <hapCarriers>
void locatePrintIBD(vector<SimDetails> &simDetails,
		    vector< vector< vector<InheritRecord> > > &hapCarriers,
		    vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		    bool sexSpecificMaps, char *ibdFile) {
  FILE *out = fopen(ibdFile, "w");
  if (!out) {
    printf("ERROR: could not open output file %s!\n", ibdFile);
    perror("open");
    exit(1);
  }

  int maxNumGens = -1; // how many generations in the largest pedigree?
  for(auto it = simDetails.begin(); it != simDetails.end(); it++) {
    if (it->numGen > maxNumGens)
      maxNumGens = it->numGen;
  }
  assert(maxNumGens > 0);

  // place to store IBD segments identified below
  vector< vector< vector<IBDRecord> > > *theSegs =
			  new vector< vector< vector<IBDRecord> > >[maxNumGens];
  if (theSegs == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  int curPed = -1;
  int curFam = -1;

  int totalFounderHaps = hapCarriers.size();
  unsigned int numChrs = geneticMap.size();
  vector<InheritRecord> overlapRecs;

  // Have a record of all individuals that inherited each founder haplotype
  // Go through this one founder haplotype and one chromosome at a time to
  // find the overlapping IBD segments
  for (int foundHapNum = 0; foundHapNum < totalFounderHaps; foundHapNum++) {
    for(unsigned int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
      sort(hapCarriers[foundHapNum][chrIdx].begin(),
	   hapCarriers[foundHapNum][chrIdx].end(), compInheritRec);
      for(auto it = hapCarriers[foundHapNum][chrIdx].begin();
	       it != hapCarriers[foundHapNum][chrIdx].end();
	       it++) {
	InheritRecord &curRec = *it;

	// <theSegs> stores segments for a given pedigree and family, and to
	// save space, it gets reused across these.
	// If we're about to start analyzing a new pedigree or family, print
	// the IBD segments found below and clear out <theSegs>.
	if ((int) curRec.ped != curPed || (int) curRec.fam != curFam) {
	  // print stored segments, locating any IBD2
	  if (curPed >= 0)
	    printIBD(out, simDetails[curPed], curFam, theSegs, geneticMap,
		     sexSpecificMaps);
	  // update:
	  curPed = curRec.ped;
	  curFam = curRec.fam;

	  clearTheSegs(simDetails[curPed], theSegs);
	}

	// Iterate over InheritRecords that came before <curRec> to search for
	// overlap (which implies an IBD segment)
	auto overIt = overlapRecs.begin();
	while (overIt != overlapRecs.end()) {
	  if (curRec.startPos < overIt->endPos) {
	    // Have IBD segment! would print, but we have to find IBD2 regions
	    // before we can do so.
	    assert(curRec.ped == overIt->ped && curRec.fam == overIt->fam);
	    assert(curRec.startPos >= overIt->startPos);

	    int endPos = min(curRec.endPos, overIt->endPos);
	    // Store away in the entry associated with the numerically lower id
	    if (curRec.gen < overIt->gen ||
		(curRec.gen == overIt->gen && curRec.branch < overIt->branch) ||
		(curRec.gen == overIt->gen && curRec.branch == overIt->branch &&
		 curRec.ind < overIt->ind)) {
	      theSegs[ curRec.gen ][ curRec.branch ][ curRec.ind ].
		emplace_back(overIt->gen, overIt->branch, overIt->ind,
		    chrIdx, curRec.startPos, endPos);
	    }
	    else {
	      theSegs[ overIt->gen ][ overIt->branch ][ overIt->ind ].
		emplace_back(curRec.gen, curRec.branch, curRec.ind,
		    chrIdx, curRec.startPos, endPos);
	    }
	    overIt++;
	  }
	  else {
	    // Segment ends before the current start: cannot match this or any
	    // upcoming segments
	    // TODO: this is a bit slow
	    overIt = overlapRecs.erase(overIt);
	  }
	}
	overlapRecs.push_back(curRec);
      }

      overlapRecs.clear();
    }
  }

  if (curPed >= 0)
    printIBD(out, simDetails[curPed], curFam, theSegs, geneticMap,
	     sexSpecificMaps);

  fclose(out);

  delete [] theSegs;
}

// print stored segments, locating any IBD2 regions
void printIBD(FILE *out, SimDetails &pedDetails, int fam,
	      vector< vector< vector<IBDRecord> > > *theSegs,
	      vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	      bool sexSpecificMaps) {
  // Go through <theSegs> and print segments for samples that were listed as
  // printed in the def file
  for(int gen = 0; gen < pedDetails.numGen; gen++) {
    for(int branch = 0; branch < pedDetails.numBranches[gen]; branch++) {
      if (pedDetails.numSampsToPrint[gen][branch] <= 0)
	continue; // no need to print IBD segment (generation not printed)

      int numNonFounders, numFounders;
      getPersonCounts(gen, pedDetails.numGen, branch,
		      pedDetails.numSampsToPrint, pedDetails.branchParents,
		      pedDetails.branchNumSpouses, numFounders, numNonFounders);
      int numPersons = numNonFounders + numFounders;
      for(int ind = 0; ind < numPersons; ind++) {
	// segments for current individual:
	vector<IBDRecord> &segs = theSegs[gen][branch][ind];
	if (segs.size() == 0)
	  continue;

	// reasoning below used to locate IBD2 requires the segments be sorted.
	// also nice to have them sorted in the output:
	sort(segs.begin(), segs.end(), compIBDRecord);

	// merge segments that are adjacent to each other
	mergeSegments(segs);
	int numSegs = segs.size();

	// Now print segments, identifying any IBD2 and printing it as such
	for(int i = 0; i < numSegs; i++) {

	  int nextI = 1; // shift for next segment
	  bool done = false;
	  while (!done) {
	    // by default should only execute the below code once
	    // certain conditions below set this to false
	    done = true;

	    // is the <i + nextI>th segment part of an IBD2 segment with the
	    // current one?
	    if (i + nextI < numSegs && // valid next segment?
		segs[i].otherGen == segs[i + nextI].otherGen &&
		segs[i].otherBranch == segs[i + nextI].otherBranch &&
		segs[i].otherInd == segs[i + nextI].otherInd &&
		segs[i].chrIdx == segs[i + nextI].chrIdx &&
		segs[i].endPos >= segs[i + nextI].startPos) {
	      // IBD2 region

	      if (pedDetails.numSampsToPrint[ segs[i].otherGen ]
					    [ segs[i].otherBranch ] <= 0)
		continue; // don't to print IBD segment (branch not printed)

	      // ensure this doesn't look like a HBD segment: shouldn't happen
	      assert(gen != segs[i].otherGen || branch != segs[i].otherBranch ||
		     ind != segs[i].otherInd);

	      // <seg[i]> should start before <seg[i + nextI]>: they're sorted
	      assert(segs[i].startPos <= segs[i + nextI].startPos);

	      // any preceding IBD1 segment?
	      if (segs[i].startPos < segs[i + nextI].startPos) {
		printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				   segs[i], /*realStart=*/ segs[i].startPos,
				   /*realEnd=*/ segs[i + nextI].startPos - 1,
				   /*type=*/ "IBD1", geneticMap,
				   sexSpecificMaps);
	      }

	      // now the IBD2 segment:
	      int ibd2End = min(segs[i].endPos, segs[i + nextI].endPos);
	      printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				 segs[i],
				 /*realStart=*/ segs[i + nextI].startPos,
				 /*realEnd=*/ ibd2End,
				 /*type=*/ "IBD2", geneticMap,
				 sexSpecificMaps);

	      // likely another segment just after the IBD2 end, but that region
	      // must be checked against later <segs>. To ensure this is done
	      // properly, we update the startPos of the continuing segment:
	      if (segs[i].endPos == segs[i + nextI].endPos) {
		// corner case: both finished, increment i so that
		// <seg[i+1]> through <seg[i + nextI]> gets skipped
		// (they're already printed)
		// NOTE: EFFECT WITH CODE JUST AFTER WHILE LOOP IS i += nextI
		i++;
		//done = true; // no need to loop (commented b/c done is true)
	      }
	      else if (segs[i].endPos > ibd2End) {
		// <segs[i]> continues; should increment nextI and loop:
		nextI++;
		// also update start position of <segs[i]> to exclude the
		// printed portion
		segs[i].startPos = ibd2End + 1;
		done = false; // loop
	      }
	      else {
		// <segs[i + nextI]> continues. This is the intuitive case: <i>
		// will increment so <segs[i + (nextI - 1) + 1]> will be
		// considered next (see code after the while loop)
		// Only require that the start of that segment doesn't include
		// the already-printed regions
		segs[i + nextI].startPos = ibd2End + 1;
		assert(segs[i + nextI].startPos <= segs[i + nextI].endPos);

		// Annoying thing is that now segs[i + nextI], with its new
		// start position, may be out of order.
		// This loop moves it to the right position:
		int initI = i + nextI;
		for(int j = 0; initI + j + 1 < numSegs &&
			       !compIBDRecord(segs[initI + j],
					       segs[initI + j + 1]) &&
			       compIBDRecord(segs[initI + j + 1],
					     segs[initI + j]); j++) {
		  // do the swap:
		  IBDRecord tmp = segs[initI + j];
		  segs[initI + j] = segs[initI + j + 1];
		  segs[initI + j + 1] = tmp;
		}
		//done = true; // no need to loop (commented b/c done is true)
	      }
	    }
	    else {
	      // IBD1 or HBD region:

	      if (pedDetails.numSampsToPrint[ segs[i].otherGen ]
					    [ segs[i].otherBranch ] <= 0)
		continue; // don't to print IBD segment (branch not printed)

	      if (gen == segs[i].otherGen && branch == segs[i].otherBranch &&
		  ind == segs[i].otherInd)
		// HBD
		printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				   segs[i],
				   /*realStart=standard=*/ segs[i].startPos,
				   /*realEnd=standard=*/ segs[i].endPos,
				   /*type=*/ "HBD", geneticMap,
				   sexSpecificMaps);
	      else
		// IBD1
		printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				   segs[i],
				   /*realStart=standard=*/ segs[i].startPos,
				   /*realEnd=standard=*/ segs[i].endPos,
				   /*type=*/ "IBD1", geneticMap,
				   sexSpecificMaps);
	    }
	  } // looping over <nextI> values

	  i += nextI - 1; // skip the already-processed IBD records
	}
      }
    }
  }
}

// Helper for printIBD(): merges adjacent IBD segments
void mergeSegments(vector<IBDRecord> &segs) {
  int numSegs = segs.size();

  int numRemoved = 0; // how many segments removed?
  for(int i = 0; i < numSegs - numRemoved; i++) {
    // Note: next valid segment is not really index <i> but
    // index <i + numRemoved>. Values from <i> to <i + numRemoved - 1> are
    // stale and present in multiple copies until the loop ends and
    // <segs> gets resized
    if (numRemoved)
      segs[i] = segs[i + numRemoved]; // shift 

    while (segs[i].otherGen == -1) { // segment that has been merged?
      // remove it!
      numRemoved++;

      if (i + numRemoved < numSegs)
	// can/should copy future segment to this position
	segs[i] = segs[i + numRemoved];
      else
	// no need to do any copying for last segment that's been merged
	// with a prior one
	break; // avoid infinite loop on last element if it was merged
    }

    // from <i + numRemoved + 1> search for a segment to merge
    // with current index <i>
    for(int j = numRemoved + 1; (i + j) < numSegs; j++) {
      if (segs[i + j].otherGen == -1)
	continue; // skip already-merged segments

      if (segs[i].otherGen != segs[i + j].otherGen ||
	  segs[i].otherBranch != segs[i + j].otherBranch ||
	  segs[i].otherInd != segs[i + j].otherInd ||
	  segs[i].chrIdx != segs[i + j].chrIdx)
	// must be same other person/chromosome in order to merge
	// due to sorting, if we encounter a different person, we know that
	// no future segments will match
	break;
      if (segs[i].endPos + 1 < segs[i + j].startPos)
	// next segment beyond end of current one, and all
	// subsequent segments necessarily start at least as far away
	// because of sorting: can't merge
	break;
      if (segs[i].endPos + 1 == segs[i + j].startPos) {
	// merge!
	segs[i].endPos = segs[i+j].endPos;
	segs[i + j].otherGen = -1;
	// continue searching for additional segments to merge wth <i>
      }
    }
  }
  // Resize <segs> post-merging
  numSegs -= numRemoved;
  segs.resize(numSegs);

  // Check for bugs
  for(int i = 0; i < numSegs; i++) {
    assert(segs[i].otherGen != -1);
    if (i < numSegs - 1)
      assert(compIBDRecord(segs[i], segs[i+1]) ||
	     !compIBDRecord(segs[i+1], segs[i]));
  }
}

// Prints the IBD segment described by the parameters to <out>
void printOneIBDSegment(FILE *out, SimDetails &pedDetails, int fam,
		  int gen, int branch, int ind, IBDRecord &seg,
		  int realStart, int realEnd, const char *type,
		  vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
		  bool sexSpecificMaps) {
  printSampleId(out, pedDetails, fam, gen, branch, ind);
  fprintf(out, "\t");
  printSampleId(out, pedDetails, fam, seg.otherGen, seg.otherBranch,
		seg.otherInd);

  fprintf(out, "\t%s\t%d\t%d\t%s", /*chrName=*/geneticMap[seg.chrIdx].first,
	  realStart, realEnd, type);

  // Find the genetic positions of the start and ends
  int ibdPhys[2] = { realStart, realEnd };
  double ibdGenet[2];
  vector<PhysGeneticPos> &thisChrMap = *geneticMap[ seg.chrIdx ].second;

  for(int pos = 0; pos < 2; pos++) {
    int left = 0, right = thisChrMap.size() - 1;

    while (true) {
      if (right - left == 1) {
	// have the left and right side: interpolate
	double interpFrac =
	  (double) (ibdPhys[pos] - thisChrMap[left].physPos) /
			(thisChrMap[right].physPos - thisChrMap[left].physPos);
	// start from the left position
	if (sexSpecificMaps)
	  ibdGenet[pos] =
		  (thisChrMap[left].mapPos[0] + thisChrMap[left].mapPos[1]) / 2;
	else
	  ibdGenet[pos] = thisChrMap[left].mapPos[0];

	// and add the factor for the distance between <left> and <right> that
	// the position is:
	if (sexSpecificMaps)
	  ibdGenet[pos] += interpFrac *
	    ((thisChrMap[right].mapPos[0] + thisChrMap[right].mapPos[1]) / 2 -
	     ibdGenet[pos]);
	else
	  ibdGenet[pos] += interpFrac * (thisChrMap[right].mapPos[0] -
								 ibdGenet[pos]);
	break;
      }
      int mid = (left + right) / 2;
      if (ibdPhys[pos] < thisChrMap[mid].physPos)
	right = mid;
      else if (ibdPhys[pos] > thisChrMap[mid].physPos)
	left = mid;
      else {
	// equal: have exact position in map:
	if (sexSpecificMaps)
	  ibdGenet[pos] =
		    (thisChrMap[mid].mapPos[0] + thisChrMap[mid].mapPos[1]) / 2;
	else
	  ibdGenet[pos] = thisChrMap[mid].mapPos[0];
	break;
      }
    }
  }

  fprintf(out, "\t%lf\t%lf\t%lf\n", ibdGenet[0], ibdGenet[1],
	  ibdGenet[1] - ibdGenet[0]);
}

// Helper for locatePrintIBD() to clear information in <theSegs>
void clearTheSegs(SimDetails &pedDetails, 
		  vector< vector< vector<IBDRecord> > > *theSegs) {
  int numGen = pedDetails.numGen;
  for(int gen = 0; gen < numGen; gen++) {
    int numBranches = pedDetails.numBranches[gen];
    if ((int) theSegs[gen].size() < numBranches)
      theSegs[gen].resize(numBranches);
    for(int branch = 0; branch < numBranches; branch++) {
      int numNonFounders, numFounders;
      getPersonCounts(gen, numGen, branch, pedDetails.numSampsToPrint,
		      pedDetails.branchParents,
		      pedDetails.branchNumSpouses, numFounders,
		      numNonFounders);
      int numPersons = numNonFounders + numFounders;
      if ((int) theSegs[gen][branch].size() < numPersons)
	theSegs[gen][branch].resize(numPersons);
      for(int ind = 0; ind < numPersons; ind++) {
	theSegs[gen][branch][ind].clear();
      }
    }
  }
}

template<typename IO_TYPE>
void FileOrGZ<IO_TYPE>::alloc_buf() {
  buf = (char *) malloc(INIT_SIZE);
  if (buf == NULL) {
    fprintf(stderr, "ERROR: out of memory\n");
    exit(1);
  }
  buf_size = INIT_SIZE;
  buf_len = 0;
}

// open <filename> using standard FILE *
template<>
bool FileOrGZ<FILE *>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = fopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

// open <filename> as a gzipped file
template<>
bool FileOrGZ<gzFile>::open(const char *filename, const char *mode) {
  // First allocate a buffer for I/O:
  alloc_buf();

  fp = gzopen(filename, mode);
  if (!fp)
    return false;
  else
    return true;
}

template<>
int FileOrGZ<FILE *>::getline() {
  return ::getline(&buf, &buf_size, fp);
}

template<>
int FileOrGZ<gzFile>::getline() {
  int n_read = 0;
  int c;

  while ((c = gzgetc(fp)) != EOF) {
    // About to have read one more, so n_read + 1 needs to be less than *n.
    // Note that we use >= not > since we need one more space for '\0'
    if (n_read + 1 >= (int) buf_size) {
      const size_t GROW = 1024;
      char *tmp_buf = (char *) realloc(buf, buf_size + GROW);
      if (tmp_buf == NULL) {
	fprintf(stderr, "ERROR: out of memory!\n");
	exit(1);
      }
      buf_size += GROW;
      buf = tmp_buf;
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

template<>
int FileOrGZ<FILE *>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  ret = vfprintf(fp, format, args);

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<gzFile>::printf(const char *format, ...) {
  va_list args;
  int ret;
  va_start(args, format);

  // NOTE: one can get automatic parallelization (in a second thread) for
  // gzipped output by opening a pipe to gzip (or bgzip). For example:
  //FILE *pipe = popen("gzip > output.vcf.gz", "w");
  // Can then fprintf(pipe, ...) as if it were a normal file.

  // gzvprintf() is slower than the code below that buffers the output.
  // Saw 13.4% speedup for processing a truncated VCF with ~50k lines and
  // 8955 samples.
//  ret = gzvprintf(fp, format, args);
  ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  if (ret < 0) {
    printf("ERROR: could not print\n");
    perror("printf");
    exit(10);
  }

  if (buf_len + ret > buf_size - 1) {
    // didn't fit the text in buf
    // first print what was in buf before the vsnprintf() call:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
    // now ensure that redoing vsnprintf() will fit in buf:
    if ((size_t) ret > buf_size - 1) {
      do { // find the buffer size that fits the last vsnprintf() call
	buf_size += INIT_SIZE;
      } while ((size_t) ret > buf_size - 1);
      free(buf);
      buf = (char *) malloc(buf_size);
      if (buf == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }
    }
    // redo:
    ret = vsnprintf(buf + buf_len, buf_size - buf_len, format, args);
  }

  buf_len += ret;
  if (buf_len >= buf_size - 1024) { // within a tolerance of MAX_BUF?
    // flush:
    gzwrite(fp, buf, buf_len);
    buf_len = 0;
  }

  va_end(args);
  return ret;
}

template<>
int FileOrGZ<FILE *>::close() {
  assert(buf_len == 0);
  // should free buf, but I know the program is about to end, so won't
  return fclose(fp);
}

template<>
int FileOrGZ<gzFile>::close() {
  if (buf_len > 0)
    gzwrite(fp, buf, buf_len);
  // should free buf, but I know the program is about to end, so won't
  return gzclose(fp);
}

// Given the simulated break points for individuals in each pedigree/family
// stored in <theSamples> and other necessary information, reads input VCF
// format data from the file named <inVCFfile> and prints the simulated
// haplotypes for each sample to <outVCFfile> in VCF format.
template<typename IO_TYPE>
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, char *inVCFfile, char *outFileBuf,
	     vector< pair<char*, vector<PhysGeneticPos>* > > &geneticMap,
	     FILE *outs[2]) {
  // open input VCF file:
  FileOrGZ<IO_TYPE> in;
  bool success = in.open(inVCFfile, "r");
  if (!success) {
    printf("\nERROR: could not open input VCF file %s!\n", inVCFfile);
    perror("open");
    exit(1);
  }

  // open output VCF file:
  FileOrGZ<IO_TYPE> out;
  success = out.open(outFileBuf, "w");
  if (!success) {
    printf("\nERROR: could not open output VCF file %s!\n", outFileBuf);
    perror("open");
    exit(1);
  }

  bernoulli_distribution genoErr( CmdLineOpts::genoErrRate );
  bernoulli_distribution homErr( CmdLineOpts::homErrRate );
  bernoulli_distribution setMissing( CmdLineOpts::missRate );
  bernoulli_distribution isPseudoHap( CmdLineOpts::pseudoHapRate );

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
  if (founderHaps == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }

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
  if (founderSamples == NULL) {
    printf("ERROR: out of memory");
    exit(5);
  }
  // number of elements of <extraSamples> to print (see below)
  unsigned int numToRetain = 0;

  while (in.getline() >= 0) { // lines of input VCF
    if (in.buf[0] == '#' && in.buf[1] == '#') {
      // header line: print to output
      out.printf("%s", in.buf);
      continue;
    }

    if (in.buf[0] == '#') {
      // header line with sample ids

      // skip all the header fields relating to meta-data:
      char *saveptr;
      char *token = strtok_r(in.buf, tab, &saveptr);
      for(int i = 1; i < 9; i++)
	token = strtok_r(NULL, tab, &saveptr);

      // now parse / store the sample ids:
      vector<char*> sampleIds;
      while ((token = strtok_r(NULL, tab, &saveptr))) {
	sampleIds.push_back(token);
      }
      numInputSamples = sampleIds.size();
      hapAlleles = new char*[numInputSamples * 2]; // 2 for diploid samples
      if (hapAlleles == NULL) {
	printf("ERROR: out of memory");
	exit(5);
      }

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
	fprintf(stderr, "\nERROR: need %d founders, but input only contains %d samples\n",
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
      out.printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

      // print sample ids:
      for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
	int numFam = simDetails[ped].numFam;
	int numGen = simDetails[ped].numGen;
	int **numSampsToPrint = simDetails[ped].numSampsToPrint;
	int *numBranches = simDetails[ped].numBranches;
	Parent **branchParents = simDetails[ped].branchParents;
	int **branchNumSpouses = simDetails[ped].branchNumSpouses;

	for(int fam = 0; fam < numFam; fam++)
	  for(int gen = 0; gen < numGen; gen++)
	    for(int branch = 0; branch < numBranches[gen]; branch++)
	      if (numSampsToPrint[gen][branch] > 0 || idOut) { // need to print?
		int numNonFounders, numFounders;
		getPersonCounts(gen, numGen, branch, numSampsToPrint,
				branchParents, branchNumSpouses, numFounders,
				numNonFounders);
		int numPersons = numNonFounders + numFounders;
		for(int ind = 0; ind < numPersons; ind++) {
		  out.printf("\t");
		  bool curIsFounder = printSampleId(NULL, simDetails[ped],
						    fam, gen, branch, ind,
						    /*printAllGens=*/ false,
						    /*gzOut=*/ &out);
		  if (idOut && curIsFounder) {
		    // print Ped-sim id to founder id file:
		    printSampleId(idOut, simDetails[ped], fam, gen, branch,ind);

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
	out.printf("\t%s", sampleIds[ sampIdx ]);
      }

      out.printf("\n");

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
    char *chrom = strtok_r(in.buf, tab, &saveptr);

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
	  fprintf(stderr, "         see variant on chromosome/contig %s, position %d\n",
		  chrom, pos);
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
    out.printf("%s\t%s", chrom, posStr);
    for(int i = 0; i < 7; i++)
      out.printf("\t%s", otherFields[i]);

    for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
      int numFam = simDetails[ped].numFam;
      int numGen = simDetails[ped].numGen;
      int **numSampsToPrint = simDetails[ped].numSampsToPrint;
      int *numBranches = simDetails[ped].numBranches;
      Parent **branchParents = simDetails[ped].branchParents;
      int **branchNumSpouses = simDetails[ped].branchNumSpouses;

      for(int fam = 0; fam < numFam; fam++)
	for(int gen = 0; gen < numGen; gen++)
	  for(int branch = 0; branch < numBranches[gen]; branch++)
	    if (numSampsToPrint[gen][branch] > 0) {
	      int numNonFounders, numFounders;
	      getPersonCounts(gen, numGen, branch, numSampsToPrint,
			      branchParents, branchNumSpouses, numFounders,
			      numNonFounders);
	      int numPersons = numNonFounders + numFounders;

	      for(int ind = 0; ind < numPersons; ind++) {

		// set to missing (according to the rate set by the user)?
		if (setMissing( randomGen )) {
		  for(int h = 0; h < 2; h++)
		    out.printf("%c.", betweenAlleles[h]);
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

		// make this a pseudo haploid genotype?
		if (CmdLineOpts::pseudoHapRate > 0) {
		  if (isPseudoHap( randomGen )) {
		    // pseudo-haploid; pick one haplotype to print
		    int printHap = coinFlip(randomGen);
		    for(int h = 0; h < 2; h++)
		      out.printf("%c%s", betweenAlleles[h],
				 founderHaps[ curFounderHaps[printHap] ]);
		  }
		  else { // not pseudo-haploid => both alleles missing:
		    for(int h = 0; h < 2; h++)
		      out.printf("%c.", betweenAlleles[h]);
		  }
		  continue;
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
		    out.printf("%c%d", betweenAlleles[h],alleles[h]);
		}
		else { // no error: print alleles from original haplotypes
		  for(int h = 0; h < 2; h++)
		    out.printf("%c%s", betweenAlleles[h],
			       founderHaps[ curFounderHaps[h] ]);
		}
	      }
	    }
    }
    // print data for the --retain_extra samples:
    for(unsigned int i = 0; i < numToRetain; i++) {
      int sampIdx = extraSamples[i];
      for(int h = 0; h < 2; h++)
	out.printf("%c%s", betweenAlleles[h], hapAlleles[ 2*sampIdx + h ]);
    }

    out.printf("\n");
  }

  out.close();
  in.close();
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
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;
    char *pedName = simDetails[ped].name;

    for(int fam = 0; fam < numFam; fam++) {
      Person ***curFamSamps = theSamples[ped][fam];
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  int numNonFounders, numFounders;
	  getPersonCounts(gen, numGen, branch, numSampsToPrint, branchParents,
			  branchNumSpouses, numFounders, numNonFounders);
	  int numPersons = numNonFounders + numFounders;

	  for(int ind = 0; ind < numPersons; ind++) {

	    // print family id (PLINK-specific) and sample id
	    fprintf(out, "%s%d ", pedName, fam+1); // family id first
	    printSampleId(out, simDetails[ped], fam, gen, branch, ind,
			  /*printAllGens=*/ true);
	    fprintf(out, " ");

	    // print parents
	    if (gen == 0 || ind < numFounders) {
	      // first generation or ind >= numNonFounders are founders, so
	      // they have no parents.
	      fprintf(out, "0 0 ");
	    }
	    else {
	      Parent pars[2]; // which branch are the two parents in
	      int parIdx[2];  // index of the Person in the branch?
	      // Is the corresponding parent a spouse of a "primary" individual?
	      bool isSpouse[2] = { false, false };
	      for(int p = 0; p < 2; p++) {
		pars[p] = branchParents[gen][branch*2 + p];
		if (pars[p].branch < 0) {
		  assert(p == 1);
		  // founders have negative indexes that start from -1, so we
		  // add 1 to get it to be 0 based and negate to get the index
		  parIdx[p] = -(pars[p].branch + 1);
		  pars[1].branch = pars[0].branch;
		  isSpouse[p] = true;
		}
		else {
		  // use "primary" person: immediately after all the spouses
		  int thisBranchNumSpouses= getBranchNumSpouses(simDetails[ped],
								pars[p].gen,
								pars[p].branch);
		  parIdx[p] = thisBranchNumSpouses;
		}
	      }

	      int par0sex = curFamSamps[ pars[0].gen ][ pars[0].branch ][ parIdx[0] ].sex;
	      for(int p = 0; p < 2; p++) {
		// print parent 0 first by default, but if parent 0 is female,
		// the following will switch and print parent 1 first
		int printPar = p ^ par0sex;
		// TODO: use printSampleId()
		if (!isSpouse[ printPar ])
		  // must be the primary person, so i1:
		  fprintf(out, "%s%d_g%d-b%d-i1 ", pedName, fam+1,
			  pars[ printPar ].gen+1, pars[ printPar ].branch+1);
		else
		  fprintf(out, "%s%d_g%d-b%d-s%d ", pedName, fam+1,
			  pars[ printPar ].gen+1, pars[ printPar ].branch+1,
			  parIdx[ printPar ]+1);
	      }
	    }

	    // print sex and phenotype; phenotype depends on whether the same
	    // gets printed:
	    int sex = theSamples[ped][fam][gen][branch][ind].sex;
	    int pheno = (numSampsToPrint[gen][branch] > 0) ? 1 : -9;
	    fprintf(out, "%d %d\n", sex+1, pheno);
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
