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

  denovoEvidence(const std::string &, const BamTools::BamRegion &, const std::string &, const std::string &, const BamTools::BamAlignment &, const std::vector<std::string> &, const std::vector<std::string> &, const std::string &);
  ~denovoEvidence();
  

  bool isDenovo_;

  bool isDenovo();
  void findHeadsAndTails();
  int32_t getRegionCoverage();
  double getStrandBias();
  int32_t getRefCount();
  int32_t getAltCount();
  std::pair<int32_t, int32_t> getGenotype();

  void populateGTFields();
  int32_t RO_;
  int32_t AO_;
  int32_t DP_;

  
 private:
  std::pair<int32_t, int32_t> genotype_;
  int32_t regionCoverage_ = 0;
  int32_t forwardCount_ = 0;
  BamTools::BamRegion region_;
  std::vector<aluHead> parentHeads_ = {};
  std::vector<polyA> parentTails_ = {};
  std::string aluClippedSeq_;  
  std::string parentBam_;
  std::string probandBam_;
  std::string jellyfishPath_;
  std::vector<std::string> refKmers_;
  std::vector<std::string> altKmers_;
  std::string jhashPath_;

  std::string refSequence_;
  std::string altSequence_;

  BamTools::BamAlignment al_;

  std::string referencePath_ = "/uufs/chpc.utah.edu/common/home/u0991464/d1/home/farrelac/references/current/human_reference_v37_decoys.fa";


  
};

#endif //__DENOVO_H__
