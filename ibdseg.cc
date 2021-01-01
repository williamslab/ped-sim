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
#include "ibdseg.h"
#include "cmdlineopts.h"
#include "datastructs.h"
#include "simulate.h"
#include "bpvcffam.h"

bool compInheritRecSamp(const InheritRecord &a, const InheritRecord &b) {
  return (a.ped < b.ped) ||
	 (a.ped == b.ped && a.fam < b.fam) ||
	 (a.ped == b.ped && a.fam == b.fam && a.gen < b.gen) ||
	 (a.ped == b.ped && a.fam == b.fam && a.gen == b.gen &&
	  a.branch < b.branch) ||
	 (a.ped == b.ped && a.fam == b.fam && a.gen == b.gen &&
	  a.branch == b.branch && a.ind < b.ind) ||
	 (a.ped == b.ped && a.fam == b.fam && a.gen == b.gen &&
	  a.branch == b.branch && a.ind == b.ind && a.startPos < b.startPos);
}

bool compInheritRecStart(const InheritRecord &a, const InheritRecord &b) {
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
// If <onlyGenetLen>, omits the start and end genetic positions and only gives
// the segment length
void locatePrintIBD(vector<SimDetails> &simDetails,
		    vector< vector< vector<InheritRecord> > > &hapCarriers,
		    GeneticMap &map, bool sexSpecificMaps, char *ibdFile,
		    bool onlyGenetLen) {
  FILE *out;
  if (ibdFile != NULL) {
    out = fopen(ibdFile, "w");
    if (!out) {
      printf("ERROR: could not open output file %s!\n", ibdFile);
      perror("open");
      exit(1);
    }
  }
  else {
    out = stdout;
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
  unsigned int numChrs = map.size();

  // Have a record of all individuals that inherited each founder haplotype
  // Go through this one founder haplotype and one chromosome at a time to
  // find the overlapping IBD and also HBD segments
  for (int foundHapNum = 0; foundHapNum < totalFounderHaps; foundHapNum++) {
    for(unsigned int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
      // First find all HBD regions
      // Do this first because multiple InheritRecords for the same sample that
      // span the same region will result in multiple IBD segments to that
      // region. To prevent this, we merge overlapping InheritRecords and
      // store away the <HBD> record. Actually do add back some IBD segments
      // so that when both samples have HBD regions, we get an IBD2 segment,
      // but without the code below, we'd get four IBD segments at such a region

      // sort by sample id to make finding HBD regions easy
      sort(hapCarriers[foundHapNum][chrIdx].begin(),
	   hapCarriers[foundHapNum][chrIdx].end(), compInheritRecSamp);

      // find the HBD regions
      for(auto it1 = hapCarriers[foundHapNum][chrIdx].begin();
	       it1 != hapCarriers[foundHapNum][chrIdx].end();
	       it1++) {
	for(auto it2 = it1 + 1;
		 it2 != hapCarriers[foundHapNum][chrIdx].end(); ) {
	  if (it1->ped != it2->ped || it1->fam != it2->fam ||
	      it1->gen != it2->gen || it1->branch != it2->branch ||
	      it1->ind != it2->ind)
	    // sorting means there are no records after <it2> that are the same
	    // sample as <it1> -- stop looping
	    break;

	  // <it1> and <it2> are for the same sample, so the same person
	  // inherited <foundHapNum>
	  // ... do the regions overlap?
	  if (it2->startPos <= it1->endPos) {
	    // <it1> and <it2> include an HBD section
	    // store the HBD region:
	    int hbdStart = max(it1->startPos, it2->startPos);
	    int hbdEnd = min(it1->endPos, it2->endPos);
	    it1->hbd.emplace_back(hbdStart, hbdEnd);

	    // Now ensure that the retained InheritRecord spans the whole region
	    it1->endPos = max(it1->endPos, it2->endPos);

	    // lastly, remove <it2> so that other samples only include one IBD
	    // segment over any HBD region
	    // TODO: this is a bit slow
	    it2 = hapCarriers[foundHapNum][chrIdx].erase(it2);
	  }
	  else {
	    it2++;
	  }
	}
      }
    }
  }

  // Now find and store all IBD (and HBD) segments
  for (int foundHapNum = 0; foundHapNum < totalFounderHaps; foundHapNum++) {
    for(unsigned int chrIdx = 0; chrIdx < numChrs; chrIdx++) {
      // here we sort by segment start to find IBD segments
      sort(hapCarriers[foundHapNum][chrIdx].begin(),
	   hapCarriers[foundHapNum][chrIdx].end(), compInheritRecStart);

      // find all the IBD segments
      for(auto it1 = hapCarriers[foundHapNum][chrIdx].begin();
	       it1 != hapCarriers[foundHapNum][chrIdx].end();
	       it1++) {
	// <theSegs> stores segments for a given pedigree and family, and to
	// save space, it gets reused across these.
	// If we're about to start analyzing a new pedigree or family, print
	// the IBD segments found below and clear out <theSegs>.
	if ((int) it1->ped != curPed || (int) it1->fam != curFam) {
	  // print stored segments, locating any IBD2
	  if (curPed >= 0)
	    printIBD(out, simDetails[curPed], curFam, theSegs, map,
		     sexSpecificMaps, onlyGenetLen);
	  // update:
	  curPed = it1->ped;
	  curFam = it1->fam;

	  clearTheSegs(simDetails[curPed], theSegs);
	}

	// Add in the HBD regions
	for(auto hbdIt = it1->hbd.begin();
		 hbdIt != it1->hbd.end();
		 hbdIt++) {
	  theSegs[ it1->gen ][ it1->branch ][ it1->ind ].
	    emplace_back(it1->gen, it1->branch, it1->ind,
			 chrIdx, hbdIt->first, hbdIt->second);
	}

	// Iterate over InheritRecords after <it1> to search for overlap
	// (which implies an IBD segment)
	for(auto it2 = it1 + 1;
		 it2 != hapCarriers[foundHapNum][chrIdx].end();
		 it2++) {
	  if (it2->startPos > it1->endPos)
	    // <it2> and later segments don't overlap <it1>: they're sorted by
	    // start position
	    break;

	  // Have IBD segment! would print, but we have to find IBD2 regions
	  // before we can do so.
	  assert(it1->ped == it2->ped && it1->fam == it2->fam);
	  assert(it1->startPos <= it2->startPos);

	  int startPos = it2->startPos;
	  int endPos = min(it1->endPos, it2->endPos);

	  auto samp1 = it1;
	  auto samp2 = it2;

	  // Store away in the entry associated with the numerically lower id
	  if (it2->gen < it1->gen ||
	      (it2->gen == it1->gen && it2->branch < it1->branch) ||
	      (it2->gen == it1->gen && it2->branch == it1->branch &&
	       it2->ind < it1->ind)) {
	    samp1 = it2;
	    samp2 = it1;
	  }

	  // We got all HBD segments above, so the two samples should be
	  // different
	  assert(it1->gen != it2->gen || it1->branch != it2->branch ||
		 it1->ind != it2->ind);

	  theSegs[ samp1->gen ][ samp1->branch ][ samp1->ind ].
	    emplace_back(samp2->gen, samp2->branch, samp2->ind,
			 chrIdx, startPos, endPos);

	  // if HBD in both <samp1> and <samp2> spans part of this region,
	  // then there is IBD2 and we need more IBD segments
	  for(auto samp1HBD = samp1->hbd.begin();
		   samp1HBD != samp1->hbd.end();
		   samp1HBD++) {

	    for(auto samp2HBD = samp2->hbd.begin();
		     samp2HBD != samp2->hbd.end();
		     samp2HBD++) {
	      if (samp2HBD->first > samp1HBD->second)
		// this and later HBD segments in <samp2> don't overlap
		// <samp1HBD>
		break;

	      if (samp2HBD->second >= samp1HBD->first &&
		  samp2HBD->first <= samp1HBD->second) {
		// overlapping HBD: add another IBD segment for IBD2
		int thisStart = max(samp1HBD->first, samp2HBD->first);
		int thisEnd = min(samp1HBD->second, samp2HBD->second);

		theSegs[ samp1->gen ][ samp1->branch ][ samp1->ind ].
		  emplace_back(samp2->gen, samp2->branch, samp2->ind,
			       chrIdx, thisStart, thisEnd);
	      }
	    } // HBD in samp2 loop
	  } // HBD in samp1 loop
	} // it2 loop
      } // it1 loop
    } // chrIdx loop
  } // foundHapNum loop

  if (curPed >= 0)
    printIBD(out, simDetails[curPed], curFam, theSegs, map, sexSpecificMaps,
	     onlyGenetLen);

  if (out != stdout)
    fclose(out);

  delete [] theSegs;
}

// print stored segments, locating any IBD2 regions
// If <onlyGenetLen>, omits the start and end genetic positions and only gives
// the segment length
void printIBD(FILE *out, SimDetails &pedDetails, int fam,
	      vector< vector< vector<IBDRecord> > > *theSegs,
	      GeneticMap &map, bool sexSpecificMaps, bool onlyGenetLen) {
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
				   /*type=*/ "IBD1", map, sexSpecificMaps,
				   onlyGenetLen);
	      }

	      // now the IBD2 segment:
	      int ibd2End = min(segs[i].endPos, segs[i + nextI].endPos);
	      printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				 segs[i],
				 /*realStart=*/ segs[i + nextI].startPos,
				 /*realEnd=*/ ibd2End,
				 /*type=*/ "IBD2", map, sexSpecificMaps,
				 onlyGenetLen);

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
				   /*type=*/ "HBD", map, sexSpecificMaps,
				   onlyGenetLen);
	      else
		// IBD1
		printOneIBDSegment(out, pedDetails, fam, gen, branch, ind,
				   segs[i],
				   /*realStart=standard=*/ segs[i].startPos,
				   /*realEnd=standard=*/ segs[i].endPos,
				   /*type=*/ "IBD1", map, sexSpecificMaps,
				   onlyGenetLen);
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
// If <onlyGenetLen>, omits the start and end genetic positions and only gives
// the segment length
void printOneIBDSegment(FILE *out, SimDetails &pedDetails, int fam,
			int gen, int branch, int ind, IBDRecord &seg,
			int realStart, int realEnd, const char *type,
			GeneticMap &map, bool sexSpecificMaps,
			bool onlyGenetLen) {
  printSampleId(out, pedDetails, fam, gen, branch, ind);
  fprintf(out, "\t");
  printSampleId(out, pedDetails, fam, seg.otherGen, seg.otherBranch,
		seg.otherInd);

  fprintf(out, "\t%s\t%d\t%d\t%s", map.chromName(seg.chrIdx),realStart, realEnd,
	  type);

  // Find the genetic positions of the start and ends
  int ibdPhys[2] = { realStart, realEnd };
  double ibdGenet[2];

  for(int pos = 0; pos < 2; pos++) {
    int left = 0, right = map.chromNumPos(seg.chrIdx) - 1;

    while (true) {
      if (right - left == 1) {
	// have the left and right side: interpolate
	double interpFrac =
	  (double) (ibdPhys[pos] - map.chromPhysPos(seg.chrIdx, left)) /
			(map.chromPhysPos(seg.chrIdx, right) -
					    map.chromPhysPos(seg.chrIdx, left));
	// start from the left position
	if (sexSpecificMaps)
	  ibdGenet[pos] =
		  (map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, left) +
		   map.chromGenetPos(seg.chrIdx, /*sex=*/ 1, left)) / 2;
	else
	  ibdGenet[pos] = map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, left);

	// and add the factor for the distance between <left> and <right> that
	// the position is:
	if (sexSpecificMaps)
	  ibdGenet[pos] += interpFrac *
	    ((map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, right) +
		   map.chromGenetPos(seg.chrIdx, /*sex=*/ 1, right)) / 2 -
	     ibdGenet[pos]);
	else
	  ibdGenet[pos] += interpFrac *
	    (map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, right) - ibdGenet[pos]);
	break;
      }
      int mid = (left + right) / 2;
      int midPhys = map.chromPhysPos(seg.chrIdx, mid);
      if (ibdPhys[pos] < midPhys)
	right = mid;
      else if (ibdPhys[pos] > midPhys)
	left = mid;
      else {
	// equal: have exact position in map:
	if (sexSpecificMaps)
	  ibdGenet[pos] =
	    (map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, mid) +
	     map.chromGenetPos(seg.chrIdx, /*sex=*/ 1, mid)) / 2;
	else
	  ibdGenet[pos] = map.chromGenetPos(seg.chrIdx, /*sex=*/ 0, mid);
	break;
      }
    }
  }

  if (onlyGenetLen)
    fprintf(out, "\t%lf\n", ibdGenet[1] - ibdGenet[0]);
  else
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
