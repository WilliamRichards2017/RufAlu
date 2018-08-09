#ifndef __DENOVO_H__
#define __DENOVO_H__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "polyATail.h"
#include "aluHead.h"


class denovoEvidence{
 public:

  denovoEvidence(std::string aluClippedSeq, BamTools::BamRegion, std::string);
  ~denovoEvidence();
  

  bool isDenovo_;

  bool isDenovo();
  void findHeadsAndTails();
  int32_t getRegionCoverage();
  double getStrandBias();
  int32_t getRefCount();
  int32_t getAltCount();
  std::pair<int32_t, int32_t> getGenotype();
  
 private:
  std::pair<int32_t, int32_t> genotype_;
  int32_t regionCoverage_ = 0;
  int32_t forwardCount_ = 0;
  BamTools::BamRegion region_;
  std::vector<aluHead> parentHeads_ = {};
  std::vector<polyA> parentTails_ = {};
  std::string aluClippedSeq_;  
  std::string parentBam_;
  
};

#endif //__DENOVO_H__
