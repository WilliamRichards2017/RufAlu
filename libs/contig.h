#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "aluHead.h"
#include "polyATail.h"
#include "util.h"

struct clipCoords;

struct contigAlignment {
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
  int32_t readsInRegion = 0;
  int32_t longestTail = 0;
  int32_t maxHash = 0;
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
