#ifndef CONTIG_H
#define CONTIG_H

#include <vector>
#include <string>

#include "util.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct contigAlignment {
  BamTools::BamAlignment alignedContig;
  BamTools::BamRegion alignedRegion;
  std::vector<BamTools::BamAlignment> supportingReads;
  bool forwardStrand = false;
  bool reverseStrand = false;
  bool doubleStranded = false;
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
