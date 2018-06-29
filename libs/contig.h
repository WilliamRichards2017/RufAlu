#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "aluHead.h"
#include "denovo.h"
#include "knownAlus.h"
#include "polyATail.h"
#include "util.h"

struct clipCoords;
class denovoEvidence;

struct contigAlignment {
  std::string bamPath;
  std::pair<std::string, int32_t> aluHit;
  std::string chrom;
  clipCoords clipCoords_;
  BamTools::BamAlignment alignedContig;
  BamTools::BamRegion alignedRegion;
  std::vector<polyA> leftBoundTails;
  std::vector<polyA> rightBoundTails;
  std::vector<aluHead> leftBoundHeads;
  std::vector<aluHead> rightBoundHeads;
  bool tailLeftBoundDS = false;
  bool tailRightBoundDS = false;
  bool tailLeftBound = false;
  bool tailRightBound = false;
  bool isDenovo = false; 
  int32_t readsInRegion = 0;
  int32_t longestTail = 0;
  int32_t maxHash = 0;
  std::vector<denovoEvidence> denovoVec_;
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
