#ifndef __ALU_HEAD_H__
#define __ALU_HEAD_H__

#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class aluHead{
 public:
  aluHead(const std::string &, const BamTools::BamAlignment, const std::string &, const int32_t &);
  ~aluHead();

  BamTools::BamAlignment al_;

  const std::vector<std::pair<std::string, int32_t> > getAluHits();
  const bool isHead();
  const bool mapReadToAluHit();

  

 private:
  bool isHead_;
  int32_t minHeadSize_;
  std::string aluFilePath_;
  std::string aluClippedSeq_;
  std::vector<std::pair<std::string, int32_t> > aluHits_;
  const bool doClipsMapToAlu(const std::vector<std::string> &);
};


#endif // __ALU_HEAD_H__
