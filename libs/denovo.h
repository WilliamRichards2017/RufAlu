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

  denovoEvidence(std::string aluClippedSeq, BamTools::BamRegion, std::vector<std::string>);
  ~denovoEvidence();
  

  bool isDenovo_;

  bool isDenovo();
  void findHeads();
  void findTails();
  
 private:
  BamTools::BamRegion region_;
  std::vector<aluHead> parentHeads_ = {};
  std::vector<polyA> parentTails_ = {};
  std::string aluClippedSeq_;  
  std::vector<std::string> parentBams_;
  
};

#endif //__DENOVO_H__
