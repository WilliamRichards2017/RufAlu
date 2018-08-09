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

/*struct contigAlignment {
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
  int32_t forwardStrands = 0;
  int32_t longestTail = 0;
  int32_t maxHash = 0;
  std::vector<denovoEvidence> denovoVec_;
};
*/

class contigAlignment{
 public:
  contigAlignment(std::string, clipCoords, std::string, BamTools::BamAlignment, std::string, BamTools::BamRegion, int32_t);
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
  bool isDenvo();
  int32_t getReadsInRegion();
  std::pair<int32_t, int32_t> getGenotype();
  int32_t getLongestTail();
  int32_t getMaxHash();		 
  std::vector<denovoEvidence> getDenovoVec();
  

 private:
  std::string bamPath_;
  std::string aluHit_;
  std::string chrom_;
  BamTools::BamAlignment alignedContig_;
  BamTools::BamRegion alignedRegion_;
  std::vector<polyA> polyATails_;
  std::vector<aluHead> aluHeads_;
  clipCoords clipCoords_;
  bool isDenovo_;
  int32_t readsInRegion_;
  std::pair<int32_t, int32_t> genotype_;
  int32_t maxHash_;
  std::vector<denovoEvidence> denovoVec_;
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
