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
	      GeneticMap &map, char *bpFile) {
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

		for(unsigned int chr = 0; chr < map.size(); chr++) {
		  // print chrom name and starting position
		  fprintf(out, " %s|%d", map.chromName(chr),
			  map.chromStartPhys(chr));
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

// Given an input VCF filename, determines whether the output should be gzipped
// or not and calls makeVCF()
void printVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	      int totalFounderHaps, const char *inVCFfile, char *outFile,
	      GeneticMap &map, FILE *outs[2]) {
  // decide whether to use gz I/O or standard, and call makeVCF() accordingly
  int inVCFlen = strlen(inVCFfile);
  if (strcmp(&CmdLineOpts::inVCFfile[ inVCFlen - 3 ], ".gz") == 0) {
    sprintf(outFile, "%s.vcf.gz", CmdLineOpts::outPrefix);
    makeVCF<gzFile>(simDetails, theSamples, totalFounderHaps,
		    CmdLineOpts::inVCFfile, /*outVCFfile=*/ outFile, map,
		    outs);
  }
  else {
    sprintf(outFile, "%s.vcf", CmdLineOpts::outPrefix);
    makeVCF<FILE *>(simDetails, theSamples, totalFounderHaps,
		    CmdLineOpts::inVCFfile, /*outVCFfile=*/ outFile, map,
		    outs);
  }
}

// Given the simulated break points for individuals in each pedigree/family
// stored in <theSamples> and other necessary information, reads input VCF
// format data from the file named <inVCFfile> and prints the simulated
// haplotypes for each sample to <outVCFfile> in VCF format.
template<typename IO_TYPE>
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, const char *inVCFfile, char *outFileBuf,
	     GeneticMap &map, FILE *outs[2]) {
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
  const char *chrName = map.chromName(chrIdx);
  int chrBegin = map.chromStartPhys(chrIdx);
  int chrEnd = map.chromEndPhys(chrIdx);

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
		  if (numSampsToPrint[gen][branch] > 0)
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
	if (chrIdx == map.size())
	  // no more chromosomes to process; will ignore remainder of VCF
	  break;
	chrName = map.chromName(chrIdx);
      }
      if (!gotSomeData || strcmp(chrom, chrName) != 0) {
	printf("ERROR: chromosome %s in VCF file either out of order or not present\n",
	       chrom);
	printf("       in genetic map\n");
	exit(5);
      }

      // update beginning / end positions for this chromosome
      chrBegin = map.chromStartPhys(chrIdx);
      chrEnd = map.chromEndPhys(chrIdx);
    }

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
      fprintf(stderr, "\nWARNING: genotyping error only implemented for markers with 2 alleles\n");
      fprintf(stderr, "         will not introduce errors at any markers with >2 alleles\n");
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

      // Now break apart the genotype into the alleles of the two haplotypes
      char *alleles[2];
      char *saveptrAlleles;
      alleles[0] = strtok_r(theGT, bar, &saveptrAlleles);
      alleles[1] = strtok_r(NULL, bar, &saveptrAlleles);

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
		theGT);
	exit(5);
      }
      if (strtok_r(NULL, bar, &saveptrAlleles) != NULL) {
	printf("ERROR: multiple '|' characters in data field\n");
	exit(5);
      }
    }

    bool fewer = numStored < numInputSamples * 2;
    bool more = sampStr != NULL;
    if (fewer || more) {
      printf("ERROR: line in VCF file has data for %s than the indicated %d samples\n",
	     (more) ? "more" : "fewer", numInputSamples);
      exit(6);
    }

    // Print this line to the output file
    out.printf("%s\t%s", chrom, posStr);
    for(int i = 0; i < 6; i++)
      out.printf("\t%s", otherFields[i]);
    out.printf("\tGT");

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
