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


class contigAlignment{

 public:
  contigAlignment(std::string, std::vector<std::string>,  std::pair<std::string, int32_t>, BamTools::BamAlignment, std::string, BamTools::BamRegion);
  ~contigAlignment();
  std::string getBamPath();
  std::pair<std::string,int32_t> getAluHit();
  std::string getChrom();
  BamTools::BamAlignment getAlignedContig();
  BamTools::BamRegion getAlignedRegion();
  std::vector<polyA> getTails();
  std::vector<aluHead> getHeads();
  std::vector<polyA> getConsensusTails();
  int32_t getConsensusTailsStartPos();
  clipCoords getClipCoords();
  bool isDenvo();
  int32_t getReadsInRegion();
  int32_t getAltCount();
  int32_t getForwardStrandCount();
  std::pair<int32_t, int32_t> getGenotype();
  int32_t getMaxHash();
  int32_t getLongestTail();
  std::vector<denovoEvidence> getDenovoVec();
  bool isDenovo();
  bool isDoubleStranded();
  

 private:
  int32_t tailSize_ = 10;
  int32_t headSize_ = 10;
  std::string bamPath_;
  std::vector<std::string> parentBamPaths_;
  std::pair<std::string, int32_t> aluHit_;
  std::string chrom_;
  BamTools::BamAlignment alignedContig_;
  BamTools::BamRegion alignedRegion_;
  std::vector<polyA> polyATails_;
  std::vector<aluHead> aluHeads_;
  clipCoords clipCoords_;
  int32_t readsInRegion_;
  int32_t altCount_;
  int32_t forwardStrandCount_;
  std::pair<int32_t, int32_t> genotype_;
  int32_t maxHash_;
  int32_t longestTail_;
  std::vector<denovoEvidence> denovoVec_;

  bool isDenovo_;
  bool doubleStranded_;

  void populateClipCoords();
  void populateMaxHash();
  void populateHeadsAndTails();
  void populateAltCount();
  void populateDenovoEvidence();
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
