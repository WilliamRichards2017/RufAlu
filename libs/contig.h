#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "polyATail.h"
#include "util.h"

struct clipCoords;

struct contigAlignment {
  std::string aluHit;
  std::string chrom;
  clipCoords clipCoords_;
  BamTools::BamAlignment alignedContig;
  BamTools::BamRegion alignedRegion;
  std::vector<polyA> leftBoundTails;
  std::vector<polyA> rightBoundTails;
  bool leftBoundDS = false;
  bool rightBoundDS = false;
  bool leftBound = false;
  bool rightBound = false;
  int32_t readsInRegion = 0;
  int32_t longestTail = 0;
  int32_t maxHash = 0;
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::string> alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
