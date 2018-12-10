#ifndef __LIBS_CONTIG_H__
#define __LIBS_CONTIG_H__

#include <string>
#include <vector>

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
  contigAlignment(std::string, std::vector<std::string>,  std::pair<std::string, int32_t>, 
		  BamTools::BamAlignment, std::string, BamTools::BamRegion, std::string, std::string);
  contigAlignment();
  ~contigAlignment();

  bool isDenovo();
  bool isHeadDoubleStranded();
  bool isReadLeftBound();
  bool isTailDoubleStranded();

  int32_t getAltCount();
  int32_t getForwardStrandCount();
  int32_t getLongestTail();
  int32_t getMaxHash();
  int32_t getReadsInRegion();
  int32_t getRefCount();

  std::string getBamPath();
  std::string getChrom();
  std::string getCigarString();

  std::pair<std::string,int32_t> getAluHit();
  std::pair<int32_t, int32_t> getGenotype();

  BamTools::BamAlignment getAlignedContig();
  BamTools::BamRegion getAlignedRegion();

  clipCoords getClipCoords();
  genotypeField getProbandGT();

  std::vector<polyA> getTails();
  std::vector<aluHead> getHeads();
  std::vector<polyA> getConsensusTails();
  std::vector<denovoEvidence> getDenovoVec();

 private:
  bool headDS_ = false;
  bool isDenovo_ = true;
  bool isLeftBound_;
  bool tailDS_ = false;

  int32_t altCount_ = 0;
  int32_t headSize_ = 10;
  int32_t forwardStrandCount_ = 0;
  int32_t longestTail_ = 0;
  int32_t maxHash_ = 0;
  int32_t readsInRegion_ = 0;
  int32_t tailSize_ = 10;

  std::string altSequence_;
  std::string bamPath_;
  std::string chrom_;
  std::string cigarString_;
  std::string fastaHackPath_;
  std::string probandAltPath_;
  std::string probandRefPath_;
  std::string referencePath_;
  std::string refSequence_;

  std::pair<int32_t, int32_t> genotype_;
  std::pair<std::string, int32_t> aluHit_;

  std::vector<std::string> altKmers_;
  std::vector<std::string> parentBamPaths_;
  std::vector<std::string> refKmers_;

  std::vector<aluHead> aluHeads_;
  std::vector<denovoEvidence> denovoVec_;
  std::vector<polyA> polyATails_;
  std::vector<polyA> consensusTails_;

  BamTools::BamAlignment alignedContig_;
  BamTools::BamRegion alignedRegion_;

  genotypeField probandGT_;
  clipCoords clipCoords_;

  void populateAltCount();
  void populateCigarString();
  void populateClipCoords();
  void populateConsensusTails();
  void populateDenovoEvidence();
  void populateHeadDS();
  void populateHeadsAndTails();
  void populateIsLeftBound();
  void populateLongestTail();
  void populateMaxHash();
  void populateProbandGT();
  void populateTailDS();
};

struct contig {
  std::string name;
  std::string seq;
  std::vector<std::pair<std::string, int32_t> > alusHit;
  std::vector<contigAlignment> contigAlignments;
};



#endif //__LIBS_CONTIG_H__
