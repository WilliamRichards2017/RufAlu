#include <string>
#include <vector>

#include "aluHead.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "minimap.h"
#include "util.h"


const std::vector<std::pair<std::string, int32_t> > aluHead::getAluHits(){
  return aluHits_;
}

const bool aluHead::doClipsMapToAlu(const std::vector<std::string> & readClips){
  for(auto rc : readClips){
    if(rc.substr(0, minHeadSize_).compare(aluClippedSeq_.substr(0, minHeadSize_)) == 0){
      return true;
    }
  }
  return false;
}

const bool aluHead::mapReadToAluHit(){
  const std::vector<std::string> readClipSeqs = util::getClipSeqs(al_);
  return aluHead::doClipsMapToAlu(readClipSeqs);
}

const bool aluHead::isHead(){
  return isHead_;
}

aluHead::aluHead(const std::string & aluClippedSeq, const BamTools::BamAlignment al,  const int32_t & minHeadSize) : aluClippedSeq_(aluClippedSeq), al_(al),  minHeadSize_(minHeadSize){
  isHead_ = aluHead::mapReadToAluHit();
}

aluHead::~aluHead(){
}
