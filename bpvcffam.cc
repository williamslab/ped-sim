// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <random>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "bpvcffam.h"
#include "cmdlineopts.h"
#include "datastructs.h"
#include "fileorgz.h"
#include "simulate.h"

// Reads the file input with the `--sexes` option that specifies the sex of
// individuals in the input VCF
void readSexes(unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
	       uint32_t sexCount[2], const char *sexesFile) {
  FILE *in = fopen(sexesFile, "r");
  if (!in) {
    fprintf(stderr, "ERROR: could not open sexes file %s!\n", sexesFile);
    perror("open");
    exit(1);
  }

  size_t bytesRead = 1024;
  char *buffer = (char *) malloc(bytesRead + 1);
  if (buffer == NULL) {
    fprintf(stderr, "ERROR: out of memory\n");
    exit(5);
  }
  const char *delim = " \t\n";

  int line = 0;
  while (getline(&buffer, &bytesRead, in) >= 0) {
    line++;

    char *id, *sex, *saveptr;
    id = strtok_r(buffer, delim, &saveptr);
    if (id == NULL)
      // blank line
      continue;

    sex = strtok_r(NULL, delim, &saveptr);
    if (sex == NULL || strtok_r(NULL, delim, &saveptr) != NULL) {
      fprintf(stderr, "ERROR: line %d in sexes file: expect two fields per line:\n",
	      line);
      fprintf(stderr, "       [id] [M/F]\n");
      exit(6);
    }
    if (sex[1] != '\0' || (sex[0] != 'M' && sex[0] != 'F')) {
      fprintf(stderr, "ERROR: line %d has a sex of %s but only 'M' or 'F' are valid\n",
	      line, sex);
      exit(6);
    }

    char *storeId = new char[ strlen(id) + 1 ]; // +1 for '\0'
    if (storeId == NULL) {
      fprintf(stderr, "ERROR: out of memory\n");
      exit(5);
    }
    strcpy(storeId, id);

    uint8_t sexIdx = (sex[0] == 'M') ? 0 : 1;
    sexes[ storeId ] = sexIdx;
    sexCount[ sexIdx ]++;
  }

  free(buffer);
  fclose(in);
}

// Prints the sample id of the given sample to <out>.
// Returns true if the sample is a founder, false otherwise.
// TODO: ideally refactor with getSampleId using templates for the first
//       argument and, e.g., `if (std::is_same<T, FILE *>::value) {  }`
template<class IO_TYPE>
bool printSampleId(FILE *out, SimDetails &pedDetails, int rep, int gen,
		   int branch, int ind, bool printAllGens,
		   FileOrGZ<IO_TYPE> *gzOut) {
  int thisBranchNumSpouses = getBranchNumSpouses(pedDetails, gen, branch);
  bool shouldPrint = pedDetails.numSampsToPrint[gen][branch] >0 || printAllGens;

  if (ind < thisBranchNumSpouses) {
    if (shouldPrint) {
      if (!gzOut)
	fprintf(out, "%s%d_g%d-b%d-s%d", pedDetails.name, rep+1, gen+1,
		branch+1, ind+1);
      else
	gzOut->printf("%s%d_g%d-b%d-s%d", pedDetails.name, rep+1, gen+1,
		      branch+1, ind+1);
    }
    return true; // is a founder
  }
  else {
    if (shouldPrint) {
      if (!gzOut)
	fprintf(out, "%s%d_g%d-b%d-i%d", pedDetails.name, rep+1, gen+1,
		branch+1, ind - thisBranchNumSpouses + 1);
      else
	gzOut->printf("%s%d_g%d-b%d-i%d", pedDetails.name, rep+1, gen+1,
		      branch+1, ind - thisBranchNumSpouses + 1);
    }
    if (gen == 0 || pedDetails.branchParents[gen][branch*2].branch < 0) {
      assert(ind - thisBranchNumSpouses == 0);
      return true; // is a founder
    }
    return false;
  }
}

// Prints the sample id of the given sample to the string <id>.
// Returns true if the sample is a founder, false otherwise.
bool getSampleId(char *id, const int MAX_LEN, SimDetails &pedDetails, int rep,
		 int gen, int branch, int ind) {
  int thisBranchNumSpouses = getBranchNumSpouses(pedDetails, gen, branch);

  if (ind < thisBranchNumSpouses) {
    int n = snprintf(id, MAX_LEN, "%s%d_g%d-b%d-s%d", pedDetails.name, rep+1,
		     gen+1, branch+1, ind+1);
    if (n >= MAX_LEN) {
      fprintf(stderr, "ERROR: a Ped-sim id is longer than %d bytes and exceeds allotted space\n",
	      MAX_LEN);
      exit(10);
    }
    return true; // is a founder
  }
  else {
    int n = snprintf(id, MAX_LEN, "%s%d_g%d-b%d-i%d", pedDetails.name, rep+1,
		     gen+1, branch+1, ind - thisBranchNumSpouses + 1);
    if (n >= MAX_LEN) {
      fprintf(stderr, "ERROR: a Ped-sim id is longer than %d bytes and exceeds allotted space\n",
	      MAX_LEN);
      exit(10);
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
	      GeneticMap &map, char *bpFile) {
  // shouldn't be possible to get here with --dry_run
  assert(!CmdLineOpts::dryRun);

  FILE *out = fopen(bpFile, "w");
  if (!out) {
    fprintf(stderr, "ERROR: could not open output file %s!\n", bpFile);
    perror("open");
    exit(1);
  }

  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    int numReps = simDetails[ped].numReps;
    int numGen = simDetails[ped].numGen;
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    for(int rep = 0; rep < numReps; rep++) {
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
		int sex = theSamples[ped][rep][gen][branch][ind].sex;
		printSampleId(out, simDetails[ped], rep, gen, branch, ind);
		fprintf(out, " s%d h%d", sex, h);

		for(unsigned int chr = 0; chr < map.size(); chr++) {
		  if (map.isX(chr) &&
		      theSamples[ped][rep][gen][branch][ind].sex == 0 &&
		      h == 0)
		    continue; // no paternal X chromosome in males

		  // print chrom name and starting position
		  fprintf(out, " %s|%d", map.chromName(chr),
			  map.chromStartPhys(chr));
		  Haplotype &curHap = theSamples[ped][rep][gen][branch][ind].
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

// Given an input VCF filename, determines whether the output should be gzipped
// or not and calls makeVCF()
int printVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, const char *inVCFfile, char *outFile,
	     GeneticMap &map, FILE *outs[2], vector<int> hapNumsBySex[2],
	     unordered_map<const char*,uint8_t,HashString,EqString> &sexes) {
  // shouldn't be possible to get here with --dry_run
  assert(!CmdLineOpts::dryRun);

  // decide whether to use gz I/O or standard, and call makeVCF() accordingly
  int inVCFlen = strlen(inVCFfile);
  if (strcmp(&CmdLineOpts::inVCFfile[ inVCFlen - 3 ], ".gz") == 0) {
    if (CmdLineOpts::nogz) {
      sprintf(outFile, "%s.vcf", CmdLineOpts::outPrefix);
      return makeVCF<gzFile, FILE *>(simDetails, theSamples, totalFounderHaps,
				     CmdLineOpts::inVCFfile,
				     /*outVCFfile=*/ outFile, map, outs,
				     hapNumsBySex, sexes);
    }
    else {
      sprintf(outFile, "%s.vcf.gz", CmdLineOpts::outPrefix);
      return makeVCF<gzFile, gzFile>(simDetails, theSamples, totalFounderHaps,
				     CmdLineOpts::inVCFfile,
				     /*outVCFfile=*/ outFile, map, outs,
				     hapNumsBySex, sexes);
    }
  }
  else {
    sprintf(outFile, "%s.vcf", CmdLineOpts::outPrefix);
    return makeVCF<FILE *, FILE *>(simDetails, theSamples, totalFounderHaps,
				   CmdLineOpts::inVCFfile,
				   /*outVCFfile=*/ outFile, map, outs,
				   hapNumsBySex, sexes);
  }
}

// Given the simulated break points for individuals in each pedigree/replicate
// stored in <theSamples> and other necessary information, reads input VCF
// format data from the file named <inVCFfile> and prints the simulated
// haplotypes for each sample to <outVCFfile> in VCF format.
template<typename I_TYPE, typename O_TYPE>
int makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	    int totalFounderHaps, const char *inVCFfile, char *outFileBuf,
	    GeneticMap &map, FILE *outs[2], vector<int> hapNumsBySex[2],
	    unordered_map<const char*,uint8_t,HashString,EqString> &sexes) {
  // shouldn't be possible to get here with --dry_run
  assert(!CmdLineOpts::dryRun);

  // open input VCF file:
  FileOrGZ<I_TYPE> in;
  bool success = in.open(inVCFfile, "r");
  if (!success) {
    fprintf(stderr, "\nERROR: could not open input VCF file %s!\n", inVCFfile);
    perror("open");
    exit(1);
  }

  // open output VCF file:
  FileOrGZ<O_TYPE> out;
  success = out.open(outFileBuf, "w");
  if (!success) {
    fprintf(stderr, "\nERROR: could not open output VCF file %s!\n",
	    outFileBuf);
    perror("open");
    exit(1);
  }

  bernoulli_distribution genoErr( CmdLineOpts::genoErrRate );
  bernoulli_distribution homErr( CmdLineOpts::homErrRate );
  bernoulli_distribution setMissing( CmdLineOpts::missRate );
  bernoulli_distribution isPseudoHap( CmdLineOpts::pseudoHapRate );
  // below, if a male is heterozygous on the X chromosome, we choose a random
  // allele. If this code is run with different values for the distributions
  // above, that can affect the outcome of the heterozygous male allele choice
  // (e.g., if the error or missing rates are 0) since it modifies the
  // state of the random number generator at the next decision point for the
  // male heterozygote. To avoid these issues, we use a different random
  // generator for the het male choice, and we define the generator now before
  // any further random numbers are generated
  mt19937 hetMaleXRandGen( randomGen );

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
    fprintf(stderr, "ERROR: out of memory\n");
    exit(5);
  }

  // iterate over chromosomes in the genetic map
  unsigned int chrIdx = 0; // index of current chromosome number;
  const char *chrName = map.chromName(chrIdx);
  int chrBegin = map.chromStartPhys(chrIdx);
  int chrEnd = map.chromEndPhys(chrIdx);

  bool gotSomeData = false;

  int numInputSamples = 0;
  // maps VCF sample index to sex (0 => M, 1 => F, 2 => unknown)
  vector<uint8_t> vcfSampleSexes; // to check X genotypes in males
  bool warnedHetMaleX = false;
  // For randomizing the assigned haplotypes or assigning VCF ids to founders
  // based on --set_founders
  vector<vector<int>> vcfIdx2FounderHap;
  vector<int> extraSamples; // Sample indexes to print for --retain_extra
  // map from founder haplotype index / 2 to VCF sample index
  // (inverse of vcfIdx2FounderHap [though stores founder haplotype index/2)
  int *founderHap2VcfIdx = new int[totalFounderHaps / 2];
  if (founderHap2VcfIdx == NULL) {
    fprintf(stderr, "ERROR: out of memory\n");
    exit(5);
  }
  else {
    for (int i = 0 ; i < totalFounderHaps / 2; i++)
      founderHap2VcfIdx[i] = -1;
  }
  // number of elements of <extraSamples> to print (see below)
  unsigned int numToRetain = 0;
  bool readMeta = false;

  while (in.getline() >= 0) { // lines of input VCF
    if (in.buf[0] == '#' && in.buf[1] == '#') {
      // header line: print to output
      out.printf("%s", in.buf);
      continue;
    }

    if (in.buf[0] == '#') {
      // header line with sample ids
      
      if (readMeta) {
	fprintf(stderr, "\n");
	fprintf(stderr, "ERROR: multiple copies of line giving sample ids: please remove all headers\n");
	fprintf(stderr, "       not at the beginning of the input VCF\n");
	exit(2);
      }
      readMeta = true;

      // skip all the header fields relating to meta-data:
      char *saveptr;
      // these fields aren't important: aren't stored
      strtok_r(in.buf, tab, &saveptr);
      for(int i = 1; i < 9; i++)
	strtok_r(NULL, tab, &saveptr);

      // map of Ped-sim sample ids to their Person* entry (only for founders)
      // this is only needed and non-empty if the below condition is true
      unordered_map<const char*,Person*,HashString,EqString> simId2Person;
      if (CmdLineOpts::setFoundersFile != NULL)
	assignSimId2Person(simId2Person, simDetails, theSamples);

      // now parse / store the VCF sample ids and either assign founders to
      // the specific VCF ids (based on the --set_founders argument) or
      // randomize haplotypes in a way that respects the sex of the input
      // samples when necessary (or do both if some founders are not listed in
      // the --set_founders file)
      vector<char*> vcfIds;
      uint32_t numDupVcfIds = getSampleIdsShuffHaps(vcfIds, vcfSampleSexes,
                                simId2Person, vcfIdx2FounderHap, outs,
                                hapNumsBySex, sexes, saveptr, totalFounderHaps,
                                tab);

      numInputSamples = vcfIds.size();
      hapAlleles = new char*[numInputSamples * 2]; // 2 for diploid samples
      if (hapAlleles == NULL) {
	fprintf(stderr, "ERROR: out of memory\n");
	exit(5);
      }

      // Store ids for all samples not assigned to a simulated founder
      // and fill in founderHap2VcfIdx array
      for(int i = 0; i < numInputSamples; i++) {
	bool assignedToFounder = false;
	for (int founderHap : vcfIdx2FounderHap[i]) {
	  if (founderHap < totalFounderHaps) {
	    assignedToFounder = true;
	    founderHap2VcfIdx[founderHap / 2] = i;
	  }
	}

	if (!assignedToFounder)
	  extraSamples.push_back(i);
      }

      uint32_t numExtraSamples = extraSamples.size();
      assert(numExtraSamples == numInputSamples -
                                          (totalFounderHaps/2 - numDupVcfIds));

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
	  fprintf(stderr, "ERROR: could not open found ids file %s!\n",
		  outFileBuf);
	  perror("open");
	  exit(1);
	}
      }

      for(int o = 0; o < 2; o++) {
	fprintf(outs[o], "done.\n"); // initial scan of VCF file (see main())
	fprintf(outs[o], "  Input contains %d samples, using %d as founders, and retaining %d\n",
		numInputSamples, totalFounderHaps/2-numDupVcfIds, numToRetain);
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
	int numReps = simDetails[ped].numReps;
	int numGen = simDetails[ped].numGen;
	int **numSampsToPrint = simDetails[ped].numSampsToPrint;
	int *numBranches = simDetails[ped].numBranches;
	Parent **branchParents = simDetails[ped].branchParents;
	int **branchNumSpouses = simDetails[ped].branchNumSpouses;

	for(int rep = 0; rep < numReps; rep++)
	  for(int gen = 0; gen < numGen; gen++)
	    for(int branch = 0; branch < numBranches[gen]; branch++)
	      if (numSampsToPrint[gen][branch] > 0 || idOut) { // need to print?
		int numNonFounders, numFounders;
		getPersonCounts(gen, numGen, branch, numSampsToPrint,
				branchParents, branchNumSpouses, numFounders,
				numNonFounders);
		int numPersons = numNonFounders + numFounders;
		for(int ind = 0; ind < numPersons; ind++) {
		  if (numSampsToPrint[gen][branch] > 0)
		    out.printf("\t");
		  bool curIsFounder = printSampleId(NULL, simDetails[ped],
						    rep, gen, branch, ind,
						    /*printAllGens=*/ false,
						    /*gzOut=*/ &out);
		  if (idOut && curIsFounder) {
		    // print Ped-sim id to founder id file:
		    printSampleId(idOut, simDetails[ped], rep, gen, branch, ind,
				  /*(always print)=*/true);

		    // since males on the X chromosome only have a defined
		    // haplotype for haps index 1, we use that index
		    int hapNum = theSamples[ped][rep][gen][branch][ind].
				      haps[1][/*chrIdx=*/0].front().foundHapNum;
		    hapNum--; // hap index 1 is an odd number, so we decrement
		    assert(hapNum % 2 == 0);
		    int founderIdx = founderHap2VcfIdx[ hapNum / 2 ];
		    assert(founderIdx >= 0);
		    fprintf(idOut, "\t%s\n", vcfIds[ founderIdx ]);
		  }

		}
	      }
      }

      // print the ids for the --retain_extra samples:
      for(unsigned int i = 0; i < numToRetain; i++) {
	int sampIdx = extraSamples[i];
	out.printf("\t%s", vcfIds[ sampIdx ]);
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
	if (chrIdx == map.size())
	  // no more chromosomes to process; will ignore remainder of VCF
	  break;
	chrName = map.chromName(chrIdx);
      }
      if (!gotSomeData || strcmp(chrom, chrName) != 0) {
	fprintf(stderr, "\nERROR: chromosome %s in VCF file either out of order or not present\n",
		chrom);
	fprintf(stderr, "       in genetic map\n");
	exit(5);
      }

      // update beginning / end positions for this chromosome
      chrBegin = map.chromStartPhys(chrIdx);
      chrEnd = map.chromEndPhys(chrIdx);
    }

    if (sexes.size() == 0 && map.isX(chrIdx))
      continue;

    gotSomeData = true;

    char *posStr = strtok_r(NULL, tab, &saveptr);
    int pos = atoi(posStr);
    if (pos < chrBegin || pos > chrEnd)
      continue; // no genetic map information for this position: skip

    // read/save the ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT fields:
    char *otherFields[6];
    for(int i = 0; i < 6; i++)
      otherFields[i] = strtok_r(NULL, tab, &saveptr);

    // Find index of GT field in format string
    char *formatStr = strtok_r(NULL, tab, &saveptr);
    char *saveptr2;
    char *formatCurField = strtok_r(formatStr, ":", &saveptr2);
    int gtField = 0;
    for ( ; formatCurField != NULL && strcmp(formatCurField, "GT") != 0;
         gtField++) {
      formatCurField = strtok_r(NULL, ":", &saveptr2);
    }
    if (formatCurField == NULL) {
      fprintf(stderr, "ERROR: in VCF: no GT field at CHROM %s, POS %s\n",
	      chrom, posStr);
      exit(6);
    }

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
      for(int o = 0; o < 2; o++) {
	fprintf(outs[o], "\nWARNING: genotyping error only implemented for markers with 2 alleles\n");
	fprintf(outs[o], "         will not introduce errors at any markers with >2 alleles\n");
      }
    }

    // read in/store the haplotypes
    int inputIndex = 0;
    int numStored = 0;
    char *sampStr;
    while((sampStr = strtok_r(NULL, tab, &saveptr)) &&
					      numStored < numInputSamples * 2) {
      // tokenize <sampStr> on ":" until we reach the <gtField>th entry
      char *saveptrGT;
      char *theGT = strtok_r(sampStr, ":", &saveptrGT);
      for(int tokenIndex = 0; tokenIndex < gtField; tokenIndex++)
	theGT = strtok_r(NULL, ":", &saveptrGT);

      for(int c = 0; theGT[c] != '\0'; c++) {
	if (theGT[c] == '/') {
	  fprintf(stderr, "\n\nERROR: detected unphased genotype; input VCF must be phased.\n");
	  fprintf(stderr, "       Prematurely truncating output VCF.\n");
	  fprintf(stderr, "       See variant on chromosome/contig %s, position %d\n",
		  chrom, pos);
	  out.close();
	  in.close();
	  return 1;
	}
      }

      // Now break apart the genotype into the alleles of the two haplotypes
      char *alleles[2];
      char *saveptrAlleles;
      alleles[0] = strtok_r(theGT, bar, &saveptrAlleles);
      alleles[1] = strtok_r(NULL, bar, &saveptrAlleles);

      if (map.isX(chrIdx) && vcfSampleSexes[inputIndex] == 0) {
	if (alleles[1] == NULL) {
	  alleles[1] = alleles[0];
	}
	else if (strcmp(alleles[0], alleles[1]) != 0) {
	  if (!warnedHetMaleX) {
	    for(int o = 0; o < 2; o++) {
	      fprintf(outs[o], "\nWARNING: heterozygous male X genotype found\n");
	      fprintf(outs[o], "         will pick allele from first haplotype\n");
	    }
	    warnedHetMaleX = true;
	  }
	  // Picking a random allele creates switch errors; although males
	  // shouldn't be heterozygous on the X, the haplotypes they have are
	  // valid and Ped-sim should not switch between them willy nilly.
	  //int keepHap = coinFlip(hetMaleXRandGen);
	  int keepHap = 0;
	  alleles[ 1 ^ keepHap ] = alleles[keepHap];
	}
      }

      // error check:
      if (alleles[1] == NULL) {
	fprintf(stderr, "ERROR: VCF contains genotype %s, which is not phased or contains one haplotype\n",
		theGT);
	fprintf(stderr, "       this is only allowed for males (input with --sexes) on the X chromosome\n");
	exit(5);
      }
      if (strtok_r(NULL, bar, &saveptrAlleles) != NULL) {
	fprintf(stderr, "ERROR: multiple '|' characters in data field\n");
	exit(5);
      }

      for(int h = 0; h < 2; h++) {
	if (alleles[h][0] == '.') {
	  fprintf(stderr, "\nERROR: simulator currently requires all positions to be non-missing\n");
	  fprintf(stderr, "         see variant on chromosome/contig %s, position %d\n",
		  chrom, pos);
	  exit(5);
	}
	hapAlleles[numStored++] = alleles[h];
      }

      // if this VCF sample has been assigned to >= 1 Ped-sim founders, store
      // the alleles we read in that/those founder haplotypes
      vector<int> &founderIndices = vcfIdx2FounderHap[ inputIndex ];
      inputIndex++;

      for(int founderIndex : founderIndices) {
	if (founderIndex < totalFounderHaps) {
	  for(int h = 0; h < 2; h++)
	    founderHaps[founderIndex + h] = alleles[h];
	}
      }

    }

    bool fewer = numStored < numInputSamples * 2;
    bool more = sampStr != NULL;
    if (fewer || more) {
      fprintf(stderr, "ERROR: line in VCF file has data for %s than the indicated %d samples\n",
	      (more) ? "more" : "fewer", numInputSamples);
      exit(6);
    }

    // Print this line to the output file
    out.printf("%s\t%s", chrom, posStr);
    for(int i = 0; i < 6; i++)
      out.printf("\t%s", otherFields[i]);
    out.printf("\tGT");

    for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
      int numReps = simDetails[ped].numReps;
      int numGen = simDetails[ped].numGen;
      int **numSampsToPrint = simDetails[ped].numSampsToPrint;
      int *numBranches = simDetails[ped].numBranches;
      Parent **branchParents = simDetails[ped].branchParents;
      int **branchNumSpouses = simDetails[ped].branchNumSpouses;

      for(int rep = 0; rep < numReps; rep++)
	for(int gen = 0; gen < numGen; gen++)
	  for(int branch = 0; branch < numBranches[gen]; branch++)
	    if (numSampsToPrint[gen][branch] > 0) {
	      int numNonFounders, numFounders;
	      getPersonCounts(gen, numGen, branch, numSampsToPrint,
			      branchParents, branchNumSpouses, numFounders,
			      numNonFounders);
	      int numPersons = numNonFounders + numFounders;

	      for(int ind = 0; ind < numPersons; ind++) {
		int numHaps = 2;
		bool maleX = false;
		if (map.isX(chrIdx) &&
		    theSamples[ped][rep][gen][branch][ind].sex == 0) {
		  maleX = true;
		  // male X: haploid output per VCF spec
		  numHaps = 1;
		}

		// set to missing (according to the rate set by the user)?
		if (setMissing( randomGen )) {
		  for(int h = 0; h < numHaps; h++)
		    out.printf("%c.", betweenAlleles[h]);
		  continue; // done printing genotype data for this sample
		}

		// non-missing genotype: print, possibly with a genotyping error

		// get founder haps for the current sample
		uint32_t curFounderHaps[2];
		for(int h = 0; h < 2; h++) {
		  if (maleX && h == 0) {
		    // only h == 1 valid for males on chrX
		    curFounderHaps[h] = (uint32_t) -1; // placeholder: fix below
		    continue;
		  }

		  Haplotype &curHap = theSamples[ped][rep][gen][branch][ind].
								haps[h][chrIdx];
		  while (curHap.front().endPos < pos) {
		    pop_front(curHap);
		  }
		  assert(curHap.front().endPos >= pos);
		  curFounderHaps[h] = curHap.front().foundHapNum;
		}
		if (maleX)
		  curFounderHaps[0] = curFounderHaps[1];

		// make this a pseudo haploid genotype?
		if (CmdLineOpts::pseudoHapRate > 0) {
		  if (isPseudoHap( randomGen )) {
		    // pseudo-haploid; pick one haplotype to print
		    int printHap = coinFlip(randomGen);
		    for(int h = 0; h < numHaps; h++)
		      out.printf("%c%s", betweenAlleles[h],
				 founderHaps[ curFounderHaps[printHap] ]);
		  }
		  else { // not pseudo-haploid => both alleles missing:
		    for(int h = 0; h < numHaps; h++)
		      out.printf("%c.", betweenAlleles[h]);
		  }
		  continue;
		}

		// genotyping error?
		if (genoErr( randomGen ) && numAlleles == 2) {
		  int alleles[2]; // integer allele values
		  for(int h = 0; h < numHaps; h++)
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
		      if (maleX) // ensure male homozygous on X
			alleles[ 1^alleleToFlip ] ^= 1;
		    }
		  }

		  for(int h = 0; h < numHaps; h++)
		    out.printf("%c%d", betweenAlleles[h],alleles[h]);
		}
		else { // no error: print alleles from original haplotypes
		  for(int h = 0; h < numHaps; h++)
		    out.printf("%c%s", betweenAlleles[h],
			       founderHaps[ curFounderHaps[h] ]);
		}
	      }
	    }
    }
    // print data for the --retain_extra samples:
    for(unsigned int i = 0; i < numToRetain; i++) {
      int sampIdx = extraSamples[i];
      int numHaps = 2;
      if (map.isX(chrIdx) && vcfSampleSexes[sampIdx] == 0)
	// male X: haploid output per VCF spec
	numHaps = 1;
      for(int h = 0; h < numHaps; h++)
	out.printf("%c%s", betweenAlleles[h], hapAlleles[ 2*sampIdx + h ]);
    }

    out.printf("\n");
  }

  out.close();
  in.close();

  return 0;
}

// Fills <simId2Person>, a map from Ped-sim founder ids (strings) to their
// Person* entry in <theSamples>.
void assignSimId2Person(
    unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
    vector<SimDetails> &simDetails, Person *****theSamples) {
  // all Ped-sim ids should fit in this; if not, we'll throw an error
  static const int BUF_SIZE = 1024 * 50;
  char idBuf[BUF_SIZE];

  // map of Ped-sim sample ids to their founder haplotype number (only for
  // founders)
  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    int numReps = simDetails[ped].numReps;
    int numGen = simDetails[ped].numGen;
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    for(int rep = 0; rep < numReps; rep++) {
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  int numNonFounders, numFounders;
	  getPersonCounts(gen, numGen, branch, numSampsToPrint,
			  branchParents, branchNumSpouses, numFounders,
			  numNonFounders);
	  int numPersons = numNonFounders + numFounders;
	  // TODO: optimization: should we just go to <numFounders>?
	  for(int ind = 0; ind < numPersons; ind++) {
	    // print Ped-sim id to <idBuf>
	    bool curIsFounder = getSampleId(idBuf, BUF_SIZE, simDetails[ped],
					    rep, gen, branch, ind);
            if (!curIsFounder)
              continue;

            char *storeId = new char[ strlen(idBuf) + 1 ]; // +1 for '\0'
            if (storeId == NULL) {
              fprintf(stderr, "ERROR: out of memory\n");
              exit(5);
            }
            strcpy(storeId, idBuf);

            Person *thisPerson = &theSamples[ped][rep][gen][branch][ind];
            simId2Person[storeId] = thisPerson;
	  }
	}
      }
    }
  }

}

// Reads and does associated bookkeeping for the founder sample assignments
// given in <CmdLineOpts::setFoundersFile> (if present).
// Checks that the VCF sample counts (including sexes where supplied) are
// sufficient for the simulation.
// Places founder haplotypes not assigned via <CmdLineOpts::setFoundersFile> in
// <hapNumsToShuf> for shuffling and assigning to VCF ids by the caller.
// Returns the number of founders that were assigned to a VCF id already
// assigned to another founder
uint32_t setFoundersHapsToShuf(
          unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
          unordered_map<const char*,uint32_t,HashString,EqString> &vcfId2index,
          unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
          vector<vector<int>> &vcfIdx2FounderHap,
          vector<int> hapNumsToShuf[3], vector<int> hapNumsBySex[2],
          int totalFounderHaps, uint32_t vcfSexCounts[3]) {
  // Should we assign haplotypes according to the sexes of the samples?
  bool respectSexes = sexes.size() > 0;

  unordered_set<int> assignedFounderHaps;
  int vcfIdsAssigned[3] = { 0, 0, 0 };
  // number of multiple assignments of VCF id; if a sample is assigned twice,
  // that leads to only 1 count below
  uint32_t dupVcfIds[3] = { 0, 0, 0 };

  if (CmdLineOpts::setFoundersFile) {
    FILE *in = fopen(CmdLineOpts::setFoundersFile, "r");
    if (!in) {
      fprintf(stderr, "\nERROR: could not open set_founders file %s!\n",
              CmdLineOpts::setFoundersFile);
      perror("open");
      exit(1);
    }

    size_t bytesRead = 1024;
    char *buffer = (char *) malloc(bytesRead + 1);
    if (buffer == NULL) {
      fprintf(stderr, "ERROR: out of memory\n");
      exit(5);
    }
    const char *delim = " \t\n";

    int line = 0;
    while (getline(&buffer, &bytesRead, in) >= 0) {
      line++;

      char *simId, *vcfId, *saveptr;
      // read Ped-sim id
      simId = strtok_r(buffer, delim, &saveptr);
      if (simId == NULL || simId[0] == '#')
        // blank line or comment -- skip
        continue;

      // ... and VCF id
      vcfId = strtok_r(NULL, delim, &saveptr);
      if (vcfId == NULL || strtok_r(NULL, delim, &saveptr) != NULL) {
        fprintf(stderr, "\nERROR: line %d in set_founders file: expect two fields per line:\n",
                line);
        fprintf(stderr, "       [Ped-sim id] [VCF id]\n");
        exit(6);
      }

      // look up founder haplotype for <simId>
      auto simIdEntry = simId2Person.find(simId);
      if (simIdEntry == simId2Person.end()) {
        fprintf(stderr, "\nERROR: line %d in set_founders file: %s is either not a\n",
                line, simId);
        fprintf(stderr, "       valid Ped-sim id or is not a founder\n");
        exit(5);
      }
      Person *thisPerson = simIdEntry->second;
      // since males on the X chromosome only have a defined
      // haplotype for haps index 1, we use that index
      int foundHapNum = thisPerson->haps[1][/*chrIdx=*/0].front().foundHapNum;
      foundHapNum--; // hap index 1 is an odd number, so we decrement
      assert(foundHapNum % 2 == 0);

      // ensure this founder haplotype (i.e., simId) wasn't assigned already
      auto assignEntry = assignedFounderHaps.find(foundHapNum); 
      if (assignEntry != assignedFounderHaps.end()) {
        fprintf(stderr, "\nERROR: line %d in set_founders file: Ped-sim founder id %s\n",
                line, simId);
        fprintf(stderr, "       assigned multiple times\n");
        exit(5);
      }
      assignedFounderHaps.insert(foundHapNum);

      // look up the index for <vcfId>
      auto vcfIdEntry = vcfId2index.find(vcfId);
      if (vcfIdEntry == vcfId2index.end()) {
        fprintf(stderr, "\nERROR: line %d in set_founders file: %s is not an id in the input VCF\n",
                line, vcfId);
        exit(5);
      }
      uint32_t vcfIndex = vcfIdEntry->second;

      uint8_t vcfSex = 2; // updated below if we're respecting sexes
      if (respectSexes) {
        // check that the sex of the Ped-sim id is the same as the person in the
        // VCF file
        auto sexEntry = sexes.find(vcfId);
        if (sexEntry != sexes.end()) {
          vcfSex = sexEntry->second;
          if (vcfSex != thisPerson->sex) {
            fprintf(stderr, "\nERROR: line %d in set_founders file: sexes of Ped-sim sample\n",
                    line);
            fprintf(stderr, "       %s and VCF sample %s do not match\n",
                    simId, vcfId);
            exit(5);
          }
        }
        else {
          fprintf(stderr, "\nERROR: line %d in set_founders file: --sexes option given, but VCF id\n",
                  line);
          fprintf(stderr, "       %s does not have a sex assigned in the sexes file so\n",
                  vcfId);
          fprintf(stderr, "       cannot be assigned to a founder\n");
          exit(5);
        }
      }

      vcfIdx2FounderHap[vcfIndex].push_back(foundHapNum);
      if (vcfIdx2FounderHap[vcfIndex].size() == 1)
        // only count this VCF id if <founderHapNum> was the first assignment to
        // this id
        vcfIdsAssigned[vcfSex]++;
      else
        dupVcfIds[vcfSex]++;
    }

    free(buffer);
    fclose(in);
  }

  // check for errors in sex counts
  if (respectSexes) {
    for (int sex = 0; sex <= 1; sex++) {
      // we subtract dupVcfIds[sex] since the haplotypes so assigned do not need
      // to be counted
      if (vcfSexCounts[sex] < hapNumsBySex[sex].size() - dupVcfIds[sex]) {
        fprintf(stderr, "\nERROR: need at least %lu females and %lu males to simulate, but the VCF\n",
            hapNumsBySex[1].size() - dupVcfIds[1],
            hapNumsBySex[0].size() - dupVcfIds[0]);
        fprintf(stderr, "       only contains %d females and %d males as specified by the sexes file\n",
            vcfSexCounts[1], vcfSexCounts[0]);
        fprintf(stderr, "       Note: it is always possible to run without an input VCF or to get\n");
        fprintf(stderr, "       autosomal genotypes by running without the --sexes option\n");
        exit(5);
      }
    }
  }
  else {
    // Simulating hap count / 2 individuals
    // subtract dupVcfIds[2] -- haplotypes so assigned do not need to be counted
    if (vcfIdx2FounderHap.size() < totalFounderHaps / 2 - dupVcfIds[2]) {
      fprintf(stderr, "\nERROR: need %d founders, but input only contains %lu samples\n",
	  totalFounderHaps / 2 - dupVcfIds[2], vcfIdx2FounderHap.size());
      exit(5);
    }
  }


  // assumed below
  assert((uint32_t) totalFounderHaps/2 == hapNumsBySex[0].size() + hapNumsBySex[1].size());

  // For all Ped-sim founders that weren't assigned above, add them to
  // <hapNumsToShuf>:
  for (int sex = 0; sex <= 1; sex++) {
    // If <respectSexes>, we'll put the haplotype numbers in the corresponding
    // sex index in <hapNumsToShuf>; otherwise, we're not worried about the
    // sex of the VCF smaples, and all individuals go in <hapNumsToShuf[2]>:
    uint8_t sexToAssign = (respectSexes) ? sex : 2;
    uint32_t numHapsThisSex = hapNumsBySex[sex].size();
    for (uint32_t i = 0; i < numHapsThisSex; i++) {
      // check whether this haplotype has been assigned
      auto assignEntry = assignedFounderHaps.find(hapNumsBySex[sex][i]); 
      if (assignEntry != assignedFounderHaps.end())
        // already assigned: skip
        continue;

      hapNumsToShuf[sexToAssign].push_back(hapNumsBySex[sex][i]);
    }
    assert(hapNumsToShuf[sexToAssign].size() <= vcfSexCounts[sexToAssign] - vcfIdsAssigned[sexToAssign]);
  }

  // Need the same number of entries in <hapNumsToShuf> as there are individuals
  // of each sex, though we must subtract the number of VCF individuals already
  // assigned above
  int numVCFids = vcfSexCounts[0] + vcfSexCounts[1] + vcfSexCounts[2];
  int numVCFidsAssigned = vcfIdsAssigned[0] + vcfIdsAssigned[1] + vcfIdsAssigned[2];
  int numShufHaps = (int) hapNumsToShuf[0].size() + hapNumsToShuf[1].size() + hapNumsToShuf[2].size();
  for(uint8_t sexToAssign = 0; sexToAssign <= 2; sexToAssign++) {
    while (true) {
      if (numShufHaps >= numVCFids - numVCFidsAssigned) {
        // done: should have sum(hapNumsToShuf[*].size()) == numVCFids (minus
        // VCF ids already assigned)
        assert(hapNumsToShuf[0].size() + hapNumsToShuf[1].size() + hapNumsToShuf[2].size() ==
               (uint32_t) numVCFids - numVCFidsAssigned);
        break;
      }
      if (sexToAssign < 2 &&
          hapNumsToShuf[sexToAssign].size() >=
                    vcfSexCounts[sexToAssign] - vcfIdsAssigned[sexToAssign])
        // done with <sexToAssign>
        break;

      // haplotype numbers >= totalFounderHaps are not assigned to founders, so
      // VCF ids that are assigned this number won't be assigned to founders
      hapNumsToShuf[sexToAssign].push_back(totalFounderHaps);
      numShufHaps++;
    }

    assert(hapNumsToShuf[sexToAssign].size() == vcfSexCounts[sexToAssign] - vcfIdsAssigned[sexToAssign]);
  }

  return dupVcfIds[0] + dupVcfIds[1] + dupVcfIds[2];
}

// Reads the sample ids from the VCF and stores them for printing the .ids file
// and any --retain_extra samples. Also maps them to sexes if --sexes was
// supplied. Next assigns VCF samples to founders by:
// (a) performing the assignments given in the --set_founders file, and/or
// (b) randomly assigning these samples to founders (i.e., founders not assigned
//     in any --set_founders file).
// For both the above, the code ensures that the sexes of the VCF sample and the
// founder are same when --sexes is supplied.
uint32_t getSampleIdsShuffHaps(vector<char*> &vcfIds,
        vector<uint8_t> &vcfSampleSexes,
        unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
        vector<vector<int>> &vcfIdx2FounderHap, FILE *outs[2],
        vector<int> hapNumsBySex[2],
        unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
        char *&saveptr, int totalFounderHaps, const char *tab) {
  // Need to:
  // (a) perform any assignments given in the --set_founders file and
  // (b) randomly assign the VCF samples to any unspecified founder haplotypes
  //
  // <hapNumsBySex> contains all the founder haplotype numbers distinguished
  // by founder sex. If the sexes of the VCF samples are known (in <sexes>), we
  // assign consistent with sex. Otherwise, we ignore sexes.

  // Should we assign haplotypes according to the sexes of the samples? Only
  // if we have sexes (via --sexes)
  bool respectSexes = sexes.size() > 0;

  // Is there a file specifying which VCF ids should be assigned to founders?
  bool settingFounders = CmdLineOpts::setFoundersFile != NULL;

  // How many VCF samples of each sex do we have?
  uint32_t vcfSexCounts[3] = { 0, 0, 0 };
  char *id;
  // map of VCF sample id to the sample's index in the VCF (i.e., the column);
  // only defined if <settingFounders>
  unordered_map<const char*,uint32_t,HashString,EqString> vcfId2index;
  while ((id = strtok_r(NULL, tab, &saveptr))) { // read next VCF id
    int vcfIdIndex = vcfIds.size();
    vcfIds.push_back(id);

    if (settingFounders) {
      char *storeId = new char[ strlen(id) + 1 ]; // +1 for '\0'
      if (storeId == NULL) {
	fprintf(stderr, "ERROR: out of memory\n");
	exit(5);
      }
      strcpy(storeId, id);

      vcfId2index[ storeId ] = vcfIdIndex;
    }

    uint8_t sexIdx = 2;
    if (respectSexes) {
      auto sexEntry = sexes.find(id);
      if (sexEntry != sexes.end())
	sexIdx = sexEntry->second;
    }
    vcfSampleSexes.push_back(sexIdx);
    vcfSexCounts[sexIdx]++;
  }
  int numVCFids = (int) vcfIds.size();

  if (respectSexes && vcfSexCounts[2] > 0) {
    for(int o = 0; o < 2; o++)
      fprintf(outs[o], "\nWARNING: %u samples in the VCF did not have sexes assigned\n",
              vcfSexCounts[2]);
  }

  // Make space in <vcfIdx2FounderHap> to store the mapping for each VCF sample
  vcfIdx2FounderHap.reserve(numVCFids);
  for (int i = 0; i < numVCFids; i++)
    vcfIdx2FounderHap.emplace_back();

  // For use with shuffling the founder haplotype assignments to individuals in
  // the VCF: store the sex-specific founder haplotype numbers. Once full:
  // - Elements [0] and [1] will have the same number of elements as the
  //   corresponding numbers of males and females in the VCF and include all
  //   founders (currently in <hapNumsBySex>).
  // - Element [2] contains either all individuals (if !<respectSexes>) or
  //   has the same number of entries as those in the VCF without a sex.
  // - NOTE: use of --set_founders changes the above slightly: we only include
  //   founders in <hapNumsBySex> that aren't assigned and the number of
  //   entries is reduced by the number of VCF ids assigned.
  vector<int> hapNumsToShuf[3];

  uint32_t numDupVcfIds = setFoundersHapsToShuf(simId2Person, vcfId2index,
                              sexes, vcfIdx2FounderHap, hapNumsToShuf,
                              hapNumsBySex, totalFounderHaps, vcfSexCounts);

  // Now shuffle the haplotypes. We begin by shuffling within the various sex
  // groups. We then insert these into <vcfIdx2FounderHap> below, respecting the
  // order of the sexes of the input VCF samples (as now stored in
  // <vcfSampleSexes>).
  // The corresponding sample in the VCF will be assigned the simulated founder
  // haplotype number, so we are logically shuffling the simulated founders, not
  // the input samples. (The input samples must be read in in a fixed order --
  // and we keep that order for the indexes of <vcfIdx2FounderHap>.)
  // Ultimately only the haplotype indexes that are < totalFounderHaps are for
  // simulated samples. (We may print some others if --retain_extra is in
  // place.)
  for(int sexIdx = 0; sexIdx < 3; sexIdx++)
    shuffle(hapNumsToShuf[sexIdx].begin(), hapNumsToShuf[sexIdx].end(),
	    randomGen);

  int curSSHapIdx[3] = { 0, 0, 0 };
  for(uint32_t i = 0; i < vcfSampleSexes.size(); i++) {
    if (vcfIdx2FounderHap[i].size() > 0)
      // already assigned founder(s) to this VCF id using --set_founders
      continue;
    uint8_t curSex = vcfSampleSexes[i];
    // map VCF sample index <i> to the given founder haplotype
    vcfIdx2FounderHap[i].push_back(hapNumsToShuf[curSex][curSSHapIdx[curSex]]);
    curSSHapIdx[curSex]++;
  }

  assert(vcfIdx2FounderHap.size() == vcfIds.size());

  return numDupVcfIds;
}

// print fam format file with the pedigree structure of all individuals included
// in the simulation
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
	      const char *famFile) {
  // open output fam file:
  FILE *out = fopen(famFile, "w");
  if (!out) {
    fprintf(stderr, "ERROR: could not open output fam file %s!\n", famFile);
    perror("open");
    exit(1);
  }

  for(unsigned int ped = 0; ped < simDetails.size(); ped++) {
    int numReps = simDetails[ped].numReps;
    int numGen = simDetails[ped].numGen;
    int **numSampsToPrint = simDetails[ped].numSampsToPrint;
    int *numBranches = simDetails[ped].numBranches;
    Parent **branchParents = simDetails[ped].branchParents;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;
    char *pedName = simDetails[ped].name;

    if (CmdLineOpts::dryRun)
      // for --dry_run, only generated one replicate per pedigree
      numReps = 1;

    for(int rep = 0; rep < numReps; rep++) {
      Person ***curRepSamps = theSamples[ped][rep];
      for(int gen = 0; gen < numGen; gen++) {
	for(int branch = 0; branch < numBranches[gen]; branch++) {
	  int numNonFounders, numFounders;
	  getPersonCounts(gen, numGen, branch, numSampsToPrint, branchParents,
			  branchNumSpouses, numFounders, numNonFounders);
	  int numPersons = numNonFounders + numFounders;

	  for(int ind = 0; ind < numPersons; ind++) {

	    // print family id (PLINK-specific) and sample id
	    fprintf(out, "%s%d ", pedName, rep+1); // family id first
	    printSampleId(out, simDetails[ped], rep, gen, branch, ind,
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

	      int par0sex = curRepSamps[ pars[0].gen ][ pars[0].branch ][ parIdx[0] ].sex;
	      for(int p = 0; p < 2; p++) {
		// print parent 0 first by default, but if parent 0 is female,
		// the following will switch and print parent 1 first
		int printPar = p ^ par0sex;
		// TODO: use printSampleId()
		if (!isSpouse[ printPar ])
		  // must be the primary person, so i1:
		  fprintf(out, "%s%d_g%d-b%d-i1 ", pedName, rep+1,
			  pars[ printPar ].gen+1, pars[ printPar ].branch+1);
		else
		  fprintf(out, "%s%d_g%d-b%d-s%d ", pedName, rep+1,
			  pars[ printPar ].gen+1, pars[ printPar ].branch+1,
			  parIdx[ printPar ]+1);
	      }
	    }

	    // print sex and phenotype; phenotype depends on whether the same
	    // gets printed:
	    int sex = theSamples[ped][rep][gen][branch][ind].sex;
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
