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
#include "simulate.h"
#include "cmdlineopts.h"
#include "fixedcos.h"

mt19937 randomGen;
uniform_int_distribution<int> coinFlip(0,1);
exponential_distribution<double> crossoverDist(1.0);


// Simulate data for each specified pedigree type for the number of requested
// families. Returns the number of founder haplotypes used to produce these
// simulated samples.
int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     GeneticMap &map, bool sexSpecificMaps, vector<COInterfere> &coIntf,
	     vector< vector< vector<InheritRecord> > > &hapCarriers,
	     vector<int> hapNumsBySex[2]) {
  // Note: throughout we use 0-based values for generations though the input
  // def file is 1-based

  int totalFounderHaps = 0;
  // Stores the random assignments of sex for the parents in a given generation
  // This relates to the constraints on sexes of individuals. Individuals with
  // the same index -- an index into the <sexAssignments> vector -- will be
  // assigned the same sex. And any pair of indexes i and i^1 will have
  // opposite sex.
  vector<int> sexAssignments;

#ifndef NOFIXEDCO
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
#endif // NOFIXEDCO

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
    SexConstraint **sexConstraints = simDetails[ped].sexConstraints;
    int i1Sex = simDetails[ped].i1Sex;
    int **branchNumSpouses = simDetails[ped].branchNumSpouses;

    simDetails[ped].founderOffset = totalFounderHaps;

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
	    if (sexConstraints[curGen] != NULL &&
				    sexConstraints[curGen][branch].theSex >= 0)
	      branchAssign = sexConstraints[curGen][branch].theSex;
	    else if (i1Sex >= 0)
	      branchAssign = i1Sex;
	    else if (sexConstraints[curGen] == NULL ||
				    sexConstraints[curGen][branch].set == -1) {
	      // no dependencies, just pick randomly
	      branchAssign = coinFlip(randomGen);
	    }
	    else {
	      while (sexConstraints[curGen][branch].set >=
						  (int) sexAssignments.size()) {
		int rand = coinFlip(randomGen);
		sexAssignments.push_back(rand);
		sexAssignments.push_back(rand ^ 1);
	      }
	      branchAssign =
			  sexAssignments[ sexConstraints[curGen][branch].set ];
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
	  unsigned int numChrs = map.size();
	  for(unsigned int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
	    // Chromosome start/ends
	    int chrStart = map.chromStartPhys(chrIdx);
	    int chrEnd = map.chromEndPhys(chrIdx);

	    // Trivial haplotypes for founders in current generation:
	    // no crossovers in founders
	    trivialSeg.endPos = chrEnd;

	    // Simulate the founders for this chromosome:
	    if (curGen != numGen - 1) { // no founders in the last generation
	      for(int ind = 0; ind < numFounders; ind++) {
		if (fam == 0 && CmdLineOpts::printMRCA) {
		  // number of digits for the numbers is 1 + (value+1) / 10
		  // 6 for 'g-b-i\0'
		  int branchNumSpouses =
			  getBranchNumSpouses(simDetails[ped], curGen, branch);
		  int idLength = 6 + 1 + (curGen + 1) / 10 +
				     1 + (branch + 1) / 10 +
				     1 + (ind + 1) / 10;
		  char *idSuffix = new char[ idLength ];
		  if (idSuffix == NULL) {
		    printf("ERROR: out of memory\n");
		    exit(5);
		  }
		  if (ind < branchNumSpouses)
		    sprintf(idSuffix, "g%d-b%d-s%d", curGen + 1, branch + 1,
			    ind + 1);
		  else
		    sprintf(idSuffix, "g%d-b%d-i%d", curGen + 1, branch + 1,
			    ind - branchNumSpouses + 1);
		  // same founder for two ids in a row for the two haplotypes:
		  simDetails[ped].founderIdSuffix.push_back(idSuffix);
		  simDetails[ped].founderIdSuffix.push_back(idSuffix);
		}

		for(int h = 0; h < 2; h++) { // 2 founder haplotypes per founder
		  int foundHapNum;
		  if (chrIdx == 0) {
		    foundHapNum = totalFounderHaps++;
		    int founderSex =
				  theSamples[ped][fam][curGen][branch][ind].sex;
		    if (h == 0)
		      // track which hap numbers (only the even indexed ones)
		      // are of each sex (for use when outputting VCFs):
		      hapNumsBySex[founderSex].push_back( foundHapNum );

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

		  if (map.isX(chrIdx) && h == 0 &&
		      theSamples[ped][fam][curGen][branch][ind].sex == 0)
		    // only one haplotype (maternal) for males on X
		    continue;

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

		Person &thePerson = curFamSamps[curGen][branch][ind];

		if (map.isX(chrIdx) && thePerson.sex == 0 && hapIdx == 0)
		  // only one haplotype for males on X: the one from Mom
		  // skip dad
		  continue;

		// Make space for this haplotype in the current sample:
		curFamSamps[curGen][branch][ind].haps[hapIdx].emplace_back();
#ifndef NOFIXEDCO
		if (p == 0 && chrIdx == 0 && curFixedCOidx >= 0) {
		  // assign fixed COs for <thePerson>:
		  for(int s = 0; s < 2; s++)
		    thePerson.fixedCOidxs[s] = randFixedCOs[s][curFixedCOidx];
		  curFixedCOidx++; // next person's indexes
		}
#endif // NOFIXEDCO
		Haplotype &toGen = thePerson.haps[hapIdx].back();
		generateHaplotype(toGen, theParent, map, coIntf, chrIdx,
				  hapCarriers,
				  (numSampsToPrint[curGen][branch] > 0) ? ped
									: -1,
				  fam, curGen, branch, ind,
				  thePerson.fixedCOidxs);
	      } // <parIdx> (simulate each transmitted haplotype for <ind>)
	    } // <ind>
	  } // <map> (chroms)
	} // <branch>
      } // <curGen>

      if (fam == 0)
	simDetails[ped].numFounders = totalFounderHaps -
						  simDetails[ped].founderOffset;
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
      // that person is married to the <-(*prevGenSpouseNum)[curGen][branch]>
      // other founders
      numFounders = -branchNumSpouses[curGen][branch] + 1;
    else if (curGen != numGen - 1)
      // though we haven't initialized <branchNumSpouses>,
      // will necessarily have 1 spouse (except for last generation)
      numFounders++;
  }
}

// Simulate one haplotype <toGenerate> by sampling crossovers and switching
// between the two haplotypes stored in <parent>.
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       GeneticMap &map, vector<COInterfere> &coIntf,
		       unsigned int chrIdx,
		       vector< vector< vector<InheritRecord> > > &hapCarriers,
		       int ped, int fam, int curGen, int branch, int ind,
		       unsigned int fixedCOidxs[2]) {
  // For the two haplotypes in <parent>, which segment index (in
  // parent.haps[].back()) is the current <switchMarker> position contained in?
  unsigned int curSegIdx[2] = { 0, 0 };

  // Pick haplotype for the beginning of the transmitted one:
  int curHap = coinFlip(randomGen);

  if (map.isX(chrIdx) && parent.sex == 0)
    // only one haplotype on X (the maternal) if the parent is male
    curHap = 1;


  double firstcMPos = map.chromStartGenet(chrIdx, parent.sex);
  double lastcMPos = map.chromEndGenet(chrIdx, parent.sex);
  // Get the genetic length of this chromosome (for the appropriate sex).
  // This is the last element in the genetic maps for the corresponding
  // chromosome:
  // In centiMorgans:
  double chrLength = lastcMPos - firstcMPos;
  // Locations of the crossovers
  vector<double> coLocations; // in Morgans

  // first segment starts at first valid map position
  int nextSegStart = map.chromStartPhys(chrIdx);

#ifndef NOFIXEDCO
  if (fixedCOidxs[0] == UINT_MAX) { // no fixed crossovers -- simulate:
#endif // NOFIXEDCO
    if (chrLength > 0.0 && // any genetic length? (is 0 on chrX for males)
	CmdLineOpts::interfereFile) {
      coIntf[chrIdx].simStahl(coLocations, parent.sex, randomGen);
    }
    else if (chrLength > 0.0) { // any genetic length? (is 0 on chrX for males)
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

    int mapNumPos = map.chromNumPos(chrIdx);
    for(auto it = coLocations.begin(); it != coLocations.end(); it++) {
      // Multiply by 100 to get cM:
      double cMPosNextCO = firstcMPos + (*it * 100);

      int left = 0, right = mapNumPos - 1;
      while (true) {
	if (right - left == 1) {
	  switchIdx = left; // want <switchIdx> <= than <cMPosNextCO>
	  break;
	}
	int mid = (left + right) / 2;
	double midGenet = map.chromGenetPos(chrIdx, parent.sex, mid);
	if (midGenet < cMPosNextCO)
	  left = mid;
	else if (midGenet > cMPosNextCO)
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
      double frac = (cMPosNextCO -
		     map.chromGenetPos(chrIdx, parent.sex, switchIdx)) /
	      (map.chromGenetPos(chrIdx, parent.sex, switchIdx+1) -
			      map.chromGenetPos(chrIdx, parent.sex, switchIdx));
      assert(frac >= 0.0 && frac <= 1.0);
      int switchPos = map.chromPhysPos(chrIdx, switchIdx) +
		      frac * (map.chromPhysPos(chrIdx, switchIdx+1) -
					  map.chromPhysPos(chrIdx, switchIdx));

      // copy Segments from <curHap>
      copySegs(toGenerate, parent, nextSegStart, switchPos, curSegIdx, curHap,
	       chrIdx, hapCarriers, ped, fam, curGen, branch, ind);
    }
#ifndef NOFIXEDCO
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
#endif // NOFIXEDCO


  // copy through to the end of the chromosome:
  for( ; curSegIdx[curHap] < parent.haps[curHap][chrIdx].size();
							  curSegIdx[curHap]++) {
    Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
    toGenerate.push_back(seg);

    if (ped >= 0) {
      hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	  ped, fam, curGen, branch, ind, nextSegStart, seg.endPos);
    }
    nextSegStart = seg.endPos + 1;
  }

  int prevEndPos = -1;
  for(auto it = toGenerate.begin(); it != toGenerate.end(); it++) {
    assert(it->endPos > prevEndPos);
    prevEndPos = it->endPos;
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
    if (nextSegStart > switchPos)
      // no distance to copy the segment over
      break;
    Segment &seg = parent.haps[curHap][chrIdx][ curSegIdx[curHap] ];
    if (seg.endPos >= switchPos) {
      // last segment to copy, and we will break it at <switchPos>
      toGenerate.emplace_back(seg.foundHapNum, switchPos);
      if (seg.endPos == switchPos)
	curSegIdx[curHap]++;

      if (ped >= 0) {
	hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	    ped, fam, curGen, branch, ind, nextSegStart, switchPos);
      }
      nextSegStart = switchPos + 1;
      break; // done copying
    }
    else {
      toGenerate.push_back(seg);

      if (ped >= 0) {
	hapCarriers[ seg.foundHapNum ][ chrIdx ].emplace_back(
	    ped, fam, curGen, branch, ind, nextSegStart, seg.endPos);
      }
      nextSegStart = seg.endPos + 1;
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

void deleteTheSamples(vector<SimDetails> &simDetails, Person *****theSamples) {
  for(unsigned int ped = 0; ped < simDetails.size(); ped++) { // for each ped
    int numFam = simDetails[ped].numFam;
    int numGen = simDetails[ped].numGen;
    int *numBranches = simDetails[ped].numBranches;
    for (int fam = 0; fam < numFam; fam++) {
      for(int curGen = 0; curGen < numGen; curGen++) {
	for(int branch = 0; branch < numBranches[curGen]; branch++) {
	  delete [] theSamples[ped][fam][curGen][branch];
	}
	delete [] theSamples[ped][fam][curGen];
      }
      delete [] theSamples[ped][fam];
    }
    delete [] theSamples[ped];
  }
  delete [] theSamples;
}
