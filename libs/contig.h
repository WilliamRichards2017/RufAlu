#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "polyATail.h"
#include "util.h"


struct contigAlignment {
  BamTools::BamAlignment alignedContig;
  BamTools::BamRegion alignedRegion;
  std::vector<polyA> leftBoundTails;
  std::vector<polyA> rightBoundTails;
  bool doubleStranded = false;
  bool leftBound = false;
  bool rightBound = false;
  int readsInRegion = 0;
  uint32_t longestTail = 0;
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::string> alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
