#include <algorithm>
#include <string>
#include <vector>
#include <utility>

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

std::pair<std::string, int32_t> contigAlignment::getAluHit(){
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
std::vector<polyA> contigAlignment::getConsensusTails(){
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
  return longestTail_;
}

int32_t contigAlignment::getMaxHash(){
  return maxHash_;
}

int32_t contigAlignment::getAltCount() {
  return contigAlignment::getConsensusTails().size() + contigAlignment::getHeads().size();
}

int32_t contigAlignment::getForwardStrandCount() {
  return forwardStrandCount_;
}

std::vector<denovoEvidence> contigAlignment::getDenovoVec(){
  return denovoVec_;
}

clipCoords contigAlignment::getClipCoords(){
  return clipCoords_;
}



void contigAlignment::populateHeadsAndTails(){
  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  if(!reader.Open(bamPath_)){
    std::cout << "could not open input Bam Path in conigAlignment::populateTails() for " << bamPath_ << std::endl;
    std::cout << "Exiting run with non-zero status..." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }
  
  reader.LocateIndex();
  
  if(!reader.HasIndex()){
    std::cout << "Index for " << bamPath_ << "could not be opened in contigAlignmnet::populateTails()" << std::endl;
    std::cout << "Exiting run with non-zero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  BamTools::BamRegion region = BamTools::BamRegion(alignedContig_.RefID, alignedContig_.Position, alignedContig_.RefID, alignedContig_.GetEndPosition());

  if(!reader.SetRegion(region)){
    std::cout << "could not set region for coords: " << alignedContig_.RefID << ", " << alignedContig_.Position << ", " << alignedContig_.GetEndPosition() << std::endl;
  }

  while(reader.GetNextAlignment(al)){
    ++readsInRegion_;
    if(!al.IsReverseStrand()){
      ++forwardStrandCount_;
    }

    if(readsInRegion_ < 998){
      break;
    }

    polyA tail = {al, tailSize_};
    aluHead head = {util::getClipSeqs(alignedContig_)[0], al, headSize_};
    
    if(tail.isTail()){
      polyATails_.push_back(tail);
    }

    if(head.isHead()){
      aluHeads_.push_back(head);
    }
  }
  reader.Close();
}

void contigAlignment::populateAltCount() {
  altCount_ = contigAlignment::getConsensusTails().size() + contigAlignment::getHeads().size();
}

void contigAlignment::populateDenovoEvidence(){
  for(const auto & pb : parentBamPaths_) {
    denovoEvidence de = {util::getClipSeqs(alignedContig_)[0], alignedRegion_, pb};
    denovoVec_.push_back(de);
    if(de.isDenovo()){
      isDenovo_ = true;
    }
  }
}

void contigAlignment::populateMaxHash(){
  std::vector<int32_t> peakVector = util::getPeakVector(alignedContig_);
  auto it = std::max_element(peakVector.begin(), peakVector.end());
  maxHash_ = peakVector[std::distance(peakVector.begin(), peakVector.end())];
}

void contigAlignment::populateClipCoords(){
  clipCoords_ = util::intersectPeaksAndClips(util::getPeaks(alignedContig_), util::getLocalClipCoords(alignedContig_));
}


contigAlignment::contigAlignment(std::string bamPath, std::vector<std::string> parentBamPaths, std::pair<std::string, int32_t> aluHit, BamTools::BamAlignment alignedContig, std::string chrom, BamTools::BamRegion alignedRegion) : bamPath_(bamPath), parentBamPaths_(parentBamPaths), aluHit_(aluHit), alignedContig_(alignedContig), chrom_(chrom), alignedRegion_(alignedRegion){

  contigAlignment::populateClipCoords();
  contigAlignment::populateMaxHash();
  contigAlignment::populateHeadsAndTails();
  contigAlignment::populateDenovoEvidence();
}

contigAlignment::~contigAlignment(){

}
