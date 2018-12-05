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

struct genotypeField{
  int32_t RO;
  int32_t AO;
  int32_t DP;
  std::pair<int32_t, int32_t> genotype;
};


class contigAlignment{

 public:
  contigAlignment(std::string, std::vector<std::string>,  std::pair<std::string, int32_t>, BamTools::BamAlignment, std::string, BamTools::BamRegion);
  ~contigAlignment();

  bool isDenovo();
  bool isTailDoubleStranded();
  bool isHeadDoubleStranded();
  bool isReadLeftBound();

  int32_t getReadsInRegion();
  int32_t getRefCount();
  int32_t getAltCount();
  int32_t getForwardStrandCount();
  int32_t getAltForwardStrandCount();
  int32_t getMaxHash();
  int32_t getLongestTail();

  std::string getBamPath();
  std::string getChrom();
  std::string getCigarString();
  std::pair<std::string,int32_t> getAluHit();
  std::pair<int32_t, int32_t> getGenotype();

  BamTools::BamAlignment getAlignedContig();
  BamTools::BamRegion getAlignedRegion();

  clipCoords getClipCoords();
  genotypeField getProbandGT();
  void populateProbandGT();
  std::vector<polyA> getTails();
  std::vector<aluHead> getHeads();
  std::vector<polyA> getConsensusTails();
  std::vector<denovoEvidence> getDenovoVec();
  

 private:
  int32_t tailSize_ = 10;
  int32_t headSize_ = 10;
  std::string bamPath_;
  std::vector<std::string> parentBamPaths_;
  std::pair<std::string, int32_t> aluHit_;
  std::string chrom_;
  std::string cigarString_ = "";
  BamTools::BamAlignment alignedContig_;
  BamTools::BamRegion alignedRegion_;
  std::vector<polyA> polyATails_;
  std::vector<polyA> consensusTails_;
  std::vector<aluHead> aluHeads_;
  clipCoords clipCoords_;
  int32_t readsInRegion_ = 0;
  int32_t altCount_ = 0;
  int32_t forwardStrandCount_ = 0;
  int32_t altForwardStrandCount_ = 0;
  std::pair<int32_t, int32_t> genotype_;
  int32_t maxHash_;
  int32_t longestTail_ = 0;
  std::vector<denovoEvidence> denovoVec_;

  bool isLeftBound_;
  bool isDenovo_ = false;
  bool tailDS_ = false;
  bool headDS_ = false;

  genotypeField probandGT_;
  std::string probandRefPath_;
  std::string probandAltPath_;
  std::vector<std::string> refKmers_;
  std::vector<std::string> altKmers_;
  std::string refSequence_;
  std::string altSequence_;

  std::string referencePath_ = "/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/references/current/human_reference_v37_decoys.fa";


  void populateCigarString();
  void populateClipCoords();
  void populateMaxHash();
  void populateIsLeftBound();
  void populateHeadsAndTails();
  void populateConsensusTails();
  void populateTailDS();
  void populateHeadDS();
  void populateLongestTail();
  void populateAltCount();
  void populateDenovoEvidence();
  void populateAltForwardStrandCount();
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};

#endif
