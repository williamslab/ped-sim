// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>
#include "cmdlineopts.h"

////////////////////////////////////////////////////////////////////////////////
// define/initialize static members
char  *CmdLineOpts::defFile = NULL;
char  *CmdLineOpts::mapFile = NULL;
char  *CmdLineOpts::interfereFile = NULL;
char  *CmdLineOpts::inVCFfile = NULL;
char  *CmdLineOpts::outPrefix = NULL;
bool   CmdLineOpts::autoSeed = true;
unsigned int CmdLineOpts::randSeed;
int    CmdLineOpts::printFam = 0;
int    CmdLineOpts::printBP = 0;
int    CmdLineOpts::printMRCA = 0;
int    CmdLineOpts::nogz = 0;
double CmdLineOpts::genoErrRate = 1e-3;
double CmdLineOpts::homErrRate = 0;
double CmdLineOpts::missRate = 1e-3;
double CmdLineOpts::pseudoHapRate = 0.0;
int    CmdLineOpts::keepPhase = 0;
int    CmdLineOpts::retainExtra = 0;
int    CmdLineOpts::printFounderIds = 0;
char  *CmdLineOpts::fixedCOfile = NULL;
char  *CmdLineOpts::chrX = NULL;
char  *CmdLineOpts::vcfSexesFile = NULL;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    RAND_SEED = CHAR_MAX + 1,
    INTERFERENCE,
    RETAIN_EXTRA,
    ERR_RATE,
    ERR_HOM_RATE,
    MISS_RATE,
    PSEUDO_HAP_RATE,
    FIXED_CO,
    SEXES,
  };

  // This is a local variable because whenever <interfereFile> is NULL, the
  // program uses a Poisson model. Variable needed to determine which model
  // is selected (and ensure that simulation is not done using both models).
  int poisson = 0;

  static struct option const longopts[] =
  {
    {"intf", required_argument, NULL, INTERFERENCE},
    {"pois", no_argument, &poisson, 1},
    {"seed", required_argument, NULL, RAND_SEED},
    {"sexes", required_argument, NULL, SEXES},
    {"fam", no_argument, &CmdLineOpts::printFam, 1},
    {"bp", no_argument, &CmdLineOpts::printBP, 1},
    {"mrca", no_argument, &CmdLineOpts::printMRCA, 1},
    {"nogz", no_argument, &CmdLineOpts::nogz, 1},
    {"keep_phase", no_argument, &CmdLineOpts::keepPhase, 1},
    {"founder_ids", no_argument, &CmdLineOpts::printFounderIds, 1},
    {"retain_extra", required_argument, NULL, RETAIN_EXTRA},
    {"err_rate", required_argument, NULL, ERR_RATE},
    {"err_hom_rate", required_argument, NULL, ERR_HOM_RATE},
    {"miss_rate", required_argument, NULL, MISS_RATE},
    {"pseudo_hap", required_argument, NULL, PSEUDO_HAP_RATE},
#ifndef NOFIXEDCO
    {"fixed_co", required_argument, NULL, FIXED_CO},
#endif // NOFIXEDCO
    {0, 0, 0, 0}
  };

  // option index for getopt_long()
  int optionIndex = 0;
  int c;

  bool haveGoodArgs = true;
  bool setMissRate = false;

  char optstring[80] = "d:m:i:o:X:";
  while ((c = getopt_long(argc, argv, optstring, longopts, &optionIndex))
									!= -1) {
    errno = 0;    // initially: may get errors from strtol
    char *endptr; // for strtol
    switch (c) {
      case 0:
	// flag set by getopt_long()
	break;

      case 'd':
	if (defFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of def filename\n");
	  haveGoodArgs = false;
	}
	defFile = optarg;
	break;
      case 'm':
	if (mapFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of map filename\n");
	  haveGoodArgs = false;
	}
	mapFile = optarg;
	break;
      case 'i':
	if (inVCFfile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input VCF filename\n");
	  haveGoodArgs = false;
	}
	inVCFfile = optarg;
	break;
      case 'o':
	outPrefix = optarg;
	break;
      case 'X':
	chrX = optarg;
	break;
      case RAND_SEED:
	autoSeed = false;
	randSeed = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse random seed as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case INTERFERENCE:
	interfereFile = optarg;
	break;
      case ERR_RATE:
	genoErrRate = strtod(optarg, &endptr);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --err_rate argument as floating point value\n");
	  if (errno != 0)
	    perror("strtod");
	  exit(2);
	}
	if (genoErrRate < 0 || genoErrRate > 1) {
	  fprintf(stderr, "ERROR: --err_rate value must be between 0 and 1\n");
	  exit(5);
	}
	break;
      case ERR_HOM_RATE:
	homErrRate = strtod(optarg, &endptr);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --err_hom_rate argument as floating point value\n");
	  if (errno != 0)
	    perror("strtod");
	  exit(2);
	}
	if (homErrRate < 0 || homErrRate > 1) {
	  fprintf(stderr, "ERROR: --err_hom_rate value must be between 0 and 1\n");
	  exit(5);
	}
	break;
      case MISS_RATE:
	missRate = strtod(optarg, &endptr);
	setMissRate = true;
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --miss_rate argument as floating point value\n");
	  if (errno != 0)
	    perror("strtod");
	  exit(2);
	}
	if (missRate < 0 || missRate > 1) {
	  fprintf(stderr, "ERROR: --miss_rate value must be between 0 and 1\n");
	  exit(5);
	}
	break;
      case PSEUDO_HAP_RATE:
	pseudoHapRate = strtod(optarg, &endptr);
	if (!setMissRate)
	  missRate = 0.0;
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --pseudo_hap argument as floating point value\n");
	  if (errno != 0)
	    perror("strtod");
	  exit(2);
	}
	if (pseudoHapRate < 0 || pseudoHapRate > 1) {
	  fprintf(stderr, "ERROR: --miss_rate value must be between 0 and 1\n");
	  exit(5);
	}
	break;
      case RETAIN_EXTRA:
	retainExtra = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --retain_extra argument as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
	break;
      case FIXED_CO:
	if (fixedCOfile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of fixed CO file\n");
	  haveGoodArgs = false;
	}
	fixedCOfile = optarg;
	break;
      case SEXES:
	if (vcfSexesFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of the VCF sexes file\n");
	  haveGoodArgs = false;
	}
	vcfSexesFile = optarg;
	break;

      case '?':
	// bad option; getopt_long already printed error message
        printUsage(stderr, argv[0]);
	exit(1);
	break;

      default:
	exit(1);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Check for errors in command line options

  if (defFile == NULL || mapFile == NULL || outPrefix == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: def, map, and output prefix names required\n");
    haveGoodArgs = false;
  }
  if (!poisson && !interfereFile && !fixedCOfile) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
#ifndef NOFIXEDCO
    fprintf(stderr, "ERROR: must specify crossover model, --pois or --intf, or use --fixedCOfile\n");
#else
    fprintf(stderr, "ERROR: must specify crossover model, --pois or --intf\n");
#endif // NOFIXEDCO
    haveGoodArgs = false;
  }
  else if ((poisson && interfereFile) || (poisson && fixedCOfile) ||
	   (interfereFile && fixedCOfile)) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
#ifndef NOFIXEDCO
    fprintf(stderr, "ERROR: can only use one crossover model, --pois or --intf, or --fixedCOfile\n");
#else
    fprintf(stderr, "ERROR: can only use one crossover model, --pois or --intf\n");
#endif // NOFIXEDCO
    haveGoodArgs = false;
  }

  if (missRate > 0 && pseudoHapRate > 0) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: can only use --miss_rate or --pseudo_hap for missingness, not both\n");
    haveGoodArgs = false;
  }

  if (chrX == NULL && haveGoodArgs) {
    chrX = new char[2];
    chrX[0] = 'X';
    chrX[1] = '\0';
  }

  if (!haveGoodArgs) {
    printUsage(stderr, argv[0]);
  }

  return haveGoodArgs;
}

// Prints usage message to <out>.  <programName> should be argv[0]
void CmdLineOpts::printUsage(FILE *out, char *programName) {
  fprintf(out, "\n");
  fprintf(out, "Pedigree simulator!  v%s    (Released %s)\n\n",
	  VERSION_NUMBER, RELEASE_DATE);
  fprintf(out, "Usage:\n");
  fprintf(out, "%s [ARGUMENTS]\n", programName);
  fprintf(out, "\n");
  fprintf(out, "REQUIRED ARGUMENTS:\n");
  fprintf(out, "  -d <filename>\t\tdef file describing pedigree structures to simulate\n");
  fprintf(out, "  -m <filename>\t\tgenetic map file containing either a sex averaged map\n");
  fprintf(out, "\t\t\t  or both male and female maps (format in README.md)\n");
  fprintf(out, "  -o <prefix>\t\toutput prefix (creates <prefix>.vcf, <prefix>.bp, etc.)\n");
  fprintf(out, "\t\t\t  if input VCF is gzipped, output is too\n");
  fprintf(out, " AND EITHER:\n");
  fprintf(out, "  --intf <filename>\tshape, escape values for interference model RECOMMENDED\n");
  fprintf(out, " OR:\n");
  fprintf(out, "  --pois\t\tPoisson crossover model (no interference)\n");
#ifndef NOFIXEDCO
  fprintf(out, " OR:\n");
  fprintf(out, "  --fixed_co <filename>\tfixed crossovers to use for simulating\n");
#endif // NOFIXEDCO
  fprintf(out, "\n\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  -i <filename>\t\tinput VCF containing phased samples to use as founders\n");
  fprintf(out, "\t\t\t  can be gzipped (with .gz extension) or not\n");
  fprintf(out, "\t\t\t  required for genetic data output\n");
  fprintf(out, "  -X <string>\t\tassign the label of the X chromosome (default: X)\n");
  fprintf(out, "  --sexes <filename>\tsexes (M/F) for all samples in the input VCF\n");
  fprintf(out, "\t\t\t  required if the input VCF contains X chromosome data\n");
  fprintf(out, "\t\t\t  (otherwise data for the X chromosome is not output and\n");
  fprintf(out, "\t\t\t   founder genotypes are assigned irrespective of sex)\n");
  fprintf(out, "  --fam\t\t\tprint PLINK fam file (see README.md before use)\n");
  fprintf(out, "  --bp\t\t\tprint BP file (complete haplotype transmission info)\n");
  fprintf(out, "  --mcra\t\tprint MRCA file (founder each IBD segment coalesces in)\n");
  fprintf(out, "  --nogz\t\talways print uncompressed VCF files\n");
  fprintf(out, "\n");
  fprintf(out, "  --seed <#>\t\tspecify random seed\n");
  fprintf(out, "\n");
  fprintf(out, " USED WITH -i:\n");
  fprintf(out, "  --err_rate <#>\tgenotyping error rate (default 1e-3; 0 disables)\n");
  fprintf(out, "  --err_hom_rate <#>\trate of opposite homozygote errors conditional on a\n");
  fprintf(out, "\t\t\t  genotyping error at the marker (default 0)\n");
  fprintf(out, "  --miss_rate <#>\tmissingness rate (default 1e-3; 0 disables)\n");
  fprintf(out, "  --pseudo_hap <#>\trate of pseudo-haploid sites; all other sites missing\n");
  fprintf(out, "\n");
  fprintf(out, "  --keep_phase\t\toutput VCF with phase information (defaults to unphased)\n");
  fprintf(out, "\n");
  fprintf(out, "  --founder_ids\t\tprint ids of founders to output file <prefix>.ids\n");
  fprintf(out, "\n");
  fprintf(out, "  --retain_extra <#>\toutput samples not used as founders to VCF file\n");
  fprintf(out, "\t\t\t  numeric argument indicates number to retain\n");
  fprintf(out, "\t\t\t  a negative argument will retain all unused samples\n");
  fprintf(out, "\n");
}
