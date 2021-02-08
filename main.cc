// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <random>
#include <sys/time.h>
#include "cmdlineopts.h"
#include "readdef.h"
#include "geneticmap.h"
#include "cointerfere.h"
#include "simulate.h"
#include "bpvcffam.h"
#include "ibdseg.h"
#include "fixedcos.h"

using namespace std;

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  // +13 for -everyone.fam, + 1 for \0
  int outFileLen = strlen(CmdLineOpts::outPrefix)+ 13 + 1;
  char *outFile = new char[outFileLen];
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
    if (strcmp(CmdLineOpts::chrX, "X") != 0)
      fprintf(outs[o], "  X chromosome:\t\t%s\n\n", CmdLineOpts::chrX);
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

  bool sexSpecificMaps;
  GeneticMap map(CmdLineOpts::mapFile, sexSpecificMaps); // read the genetic map

  vector<COInterfere> coIntf;
  if (CmdLineOpts::interfereFile) {
    COInterfere::read(coIntf, CmdLineOpts::interfereFile, map, sexSpecificMaps);
  }

#ifndef NOFIXEDCO
  if (CmdLineOpts::fixedCOfile) {
    FixedCOs::read(CmdLineOpts::fixedCOfile, map);
  }
#endif // NOFIXEDCO

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
  int totalFounderHaps = simulate(simDetails, theSamples, map, sexSpecificMaps,
				  coIntf, hapCarriers);
  for(int o = 0; o < 2; o++)
    fprintf(outs[o], "done.\n");

  if (CmdLineOpts::printBP) {
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "Printing break points... ");
      fflush(outs[o]);
    }
    sprintf(outFile, "%s.bp", CmdLineOpts::outPrefix);
    printBPs(simDetails, theSamples, map, /*bpFile=*/ outFile);
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "done.\n");
    }
  }

  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "Printing IBD segments... ");
    fflush(outs[o]);
  }
  sprintf(outFile, "%s.seg", CmdLineOpts::outPrefix);
  locatePrintIBD(simDetails, hapCarriers, map, sexSpecificMaps,
		 /*ibdFile=*/ outFile, /*onlyGenetLen=*/ false);
  for(int o = 0; o < 2; o++) {
    fprintf(outs[o], "done.\n");
  }

  if (CmdLineOpts::printFam) {
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "Printing fam file... ");
    fflush(stdout);
    sprintf(outFile, "%s-everyone.fam", CmdLineOpts::outPrefix);
    printFam(simDetails, theSamples, /*famFile=*/ outFile);
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "done.  (Do not use with PLINK data: see README.md)\n");
  }

  if (CmdLineOpts::inVCFfile) {
    for(int o = 0; o < 2; o++) {
      fprintf(outs[o], "Reading input VCF meta data... ");
      fflush(outs[o]);
    }
    // note: printVCF()'s callee makeVCF() print the status for generating the
    // VCF file

    int ret = printVCF(simDetails, theSamples, totalFounderHaps,
		       CmdLineOpts::inVCFfile, outFile, map, outs);

    if (ret == 0)
      for(int o = 0; o < 2; o++)
	fprintf(outs[o], "done.\n");
  }
  else {
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "\nTo simulate genetic data, must use an input VCF with %d founders.\n",
	      totalFounderHaps / 2);
  }

  fclose(log);

  // NOTE: the memory isn't freed because the OS reclaims it when Ped-sim
  // finishes

  return 0;
}
