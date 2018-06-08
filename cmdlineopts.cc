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
char  *CmdLineOpts::inVCFfile = 0;
char  *CmdLineOpts::outPrefix = NULL;
bool   CmdLineOpts::autoSeed = true;
unsigned int CmdLineOpts::randSeed;
double CmdLineOpts::genoErrRate = 1e-3;
double CmdLineOpts::homErrRate = .1;
double CmdLineOpts::missRate = 1e-3;
int    CmdLineOpts::keepPhase = 0;
int    CmdLineOpts::retainExtra = 0;
int    CmdLineOpts::printFounderIds = 0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    RAND_SEED = CHAR_MAX + 1,
    INTERFERENCE,
    RETAIN_EXTRA,
    ERR_RATE,
    ERR_HOM_RATE,
    MISS_RATE,
  };

  static struct option const longopts[] =
  {
    {"seed", required_argument, NULL, RAND_SEED},
    {"intf", required_argument, NULL, INTERFERENCE},
    {"keep_phase", no_argument, &CmdLineOpts::keepPhase, 1},
    {"founder_ids", no_argument, &CmdLineOpts::printFounderIds, 1},
    {"retain_extra", required_argument, NULL, RETAIN_EXTRA},
    {"err_rate", required_argument, NULL, ERR_RATE},
    {"err_hom_rate", required_argument, NULL, ERR_HOM_RATE},
    {"miss_rate", required_argument, NULL, MISS_RATE},
    {0, 0, 0, 0}
  };

  // option index for getopt_long()
  int optionIndex = 0;
  int c;

  bool haveGoodArgs = true;

  char optstring[80] = "d:m:i:o:";
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
      case RETAIN_EXTRA:
	retainExtra = strtol(optarg, &endptr, 10);
	if (errno != 0 || *endptr != '\0') {
	  fprintf(stderr, "ERROR: unable to parse --retain_extra argument as integer\n");
	  if (errno != 0)
	    perror("strtol");
	  exit(2);
	}
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

  if (defFile == NULL || mapFile == NULL || inVCFfile == NULL ||
							    outPrefix == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: def, map, input VCF, and output prefix names required\n");
    haveGoodArgs = false;
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
  fprintf(out, "\t\t\tor both male and female maps (format given in README.md)\n");
  fprintf(out, "  -i <filename>\t\tinput VCF containing phased samples to use as founders\n");
  fprintf(out, "\t\t\tcan be gzipped (with .gz extension) or not\n");
  fprintf(out, "  -o <prefix>\t\toutput prefix (creates <prefix>.vcf, <prefix>.bp, etc.)\n");
  fprintf(out, "\t\t\tif input VCF is gzipped, output is too\n");
  fprintf(out, "\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  --seed <#>\t\tspecify random seed\n");
  fprintf(out, "\n");
  fprintf(out, "  --intf <filename>\tshape and escape fraction values for interference model\n");
  fprintf(out, "\t\t\t(recommended)\n");
  fprintf(out, "\n");
  fprintf(out, "  --err_rate <#>\tgenotyping error rate (default 1e-3; 0 disables)\n");
  fprintf(out, "  --err_hom_rate <#>\trate of opposite homozygote errors conditional on a\n");
  fprintf(out, "\t\t\tgenotyping error at the marker (default .1)\n");
  fprintf(out, "  --miss_rate <#>\tmissingness rate (default 1e-3; 0 disables)\n");
  fprintf(out, "\n");
  fprintf(out, "  --keep_phase\t\toutput VCF with phase information (defaults to unphased)\n");
  fprintf(out, "\n");
  fprintf(out, "  --founder_ids\t\tprint ids of founders to output file <prefix>.ids\n");
  fprintf(out, "\n");
  fprintf(out, "  --retain_extra <#>\toutput samples not used as founders to VCF file\n");
  fprintf(out, "\t\t\tnumeric argument indicates number to retain\n");
  fprintf(out, "\t\t\ta negative argument will retain all unused samples\n");
  fprintf(out, "\n");
}
