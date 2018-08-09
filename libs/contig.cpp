#include <string>
#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "contig.h"
#include "knownAlus.h"
#include "aluHead.h"
#include "denovo.h"
#include "polyATail.h"
#include "util.h"


std::string contigAlignment::getBamPath(){
  return bamPath_;
}

std::string contigAlignment::getAluHit(){
  return aluHit_;
}

std::string contigAlignment::getChrom(){
  return chrom_;
}

BamTools::BamAlignment contigAlignment::getAlignedContig(){
  return alignedContig_;
}

BamTools::BamRegion contigAlignment::getAlignedRegion(){
  return alignedRegion_;
}

std::vector<polyA> contigAlignment::getTails(){
  return polyATails_;
}

//TODO: implement
std::vector<polyA> contigAligment::getConsensusTails(){
  std::vector<polyA> notImplemented;
  return notImplemented;
}

//TODO: implement
int32_t contigAlignment::getConsensusTailsStartPos(){
  int32_t notImplemented;
  return notImplemented;
}

std::vector<aluHead> contigAlignment::getHeads(){
  return aluHeads_;
}

bool contigAlignment::isDenovo(){
  return isDenovo_;
}

int32_t contigAlignment::getReadsInRegion(){
  return readsInRegion_;
}

std::pair<int32_t, int32_t> contigAlignment::getGenotype(){
  return genotype_;
}

//TODO:: implement
int32_t contigAlignment::getLongestTail(){
  int32_t notImplemented;
  return notImplemented;
}

int32_t contigAlignment::getMaxHash(){
  return maxHash_;
}

std::vector<denovoEvidence> getDenovoVec(){
  return denovoVec_;
}

contigAlignment::contigAlignment(std::string bamPath, clipCoords ccs, std::string aluHit, BamTools::BamAlignment alignedContig, std::string chrom, BamTools::BamRegion alignedRegion, int32_t maxHash) : bamPath_(bamPath), clipCoords_(ccs), aluHit_(aluHit), alignedContig_(alignedContig), chrom_(chrom), alignedRegion_(alignedRegion), maxHash_(maxHash){
  
}

~contigAlignment::contigAlignment(){
}
