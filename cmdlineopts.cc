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
char *CmdLineOpts::datFile = NULL;
char *CmdLineOpts::mapFile = NULL;
char *CmdLineOpts::inVCFfile = 0;
char *CmdLineOpts::outPrefix = NULL;
bool  CmdLineOpts::autoSeed = true;
unsigned int CmdLineOpts::randSeed;
int   CmdLineOpts::keepPhase = 0;
int   CmdLineOpts::retainExtra = 0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    RAND_SEED = CHAR_MAX + 1,
    RETAIN_EXTRA,
  };

  static struct option const longopts[] =
  {
    {"seed", required_argument, NULL, RAND_SEED},
    {"keep_phase", no_argument, &CmdLineOpts::keepPhase, 1},
    {"retain_extra", required_argument, NULL, RETAIN_EXTRA},
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
	if (datFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of dat filename\n");
	  haveGoodArgs = false;
	}
	datFile = optarg;
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

  if (datFile == NULL || mapFile == NULL || inVCFfile == NULL ||
							    outPrefix == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: dat, map, input VCF, and output prefix names required\n");
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
  fprintf(out, "  -d <filename>\t\tdat file describing pedigree structures to simulate\n");
  fprintf(out, "  -m <filename>\t\tgenetic map file containing either a sex averaged map\n");
  fprintf(out, "\t\t\tor both male and female maps (format given in README.md)\n");
  fprintf(out, "  -i <filename>\t\tinput VCF containing phased samples to use as founders\n");
  fprintf(out, "  -o <prefix>\t\toutput prefix (creates <prefix>.vcf, <prefix>.bp, etc.)\n");
  fprintf(out, "\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  --seed <#>\t\tspecify random seed\n");
  fprintf(out, "\n");
  fprintf(out, "  --keep_phase\t\toutput VCF with phase information (defaults to unphased)\n");
  fprintf(out, "\n");
  fprintf(out, "  --retain_extra <#>\toutput samples not used as founders to VCF file\n");
  fprintf(out, "\t\t\tnumeric argument indicates number to retain\n");
  fprintf(out, "\t\t\ta negative argument will retain all unused samples\n");
  fprintf(out, "\n");
}
