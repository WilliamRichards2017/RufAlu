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
void contigAlignment::populateConsensusTails(){

  int modeCount = 0;
  int startPosMode = -1;
  
  std::map<int32_t, int32_t> hashMap;
  for(int i = 0; i < polyATails_.size(); ++i){
    int startPos = polyATails_[i].getGlobalTailStart();
    ++hashMap[startPos];

    modeCount = std::max(modeCount, hashMap[startPos]);
  }
  
  for (auto & startPos : hashMap){
    if (startPos.second == modeCount){
      startPosMode = startPos.first;
    }
  }

  //  std::cout << "~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~" << std::endl;
  //std::cout << "startPosMode is: " << startPosMode << std::endl;
  //std::cout << "size of polyATail_ is: " << polyATails_.size() << std::endl;
  //std::cout << "~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~" << std::endl;

  
  for (auto & t : polyATails_){
    //std::cout << "globalTailStart is: " << t.getGlobalTailStart() << std::endl;
    if(t.getGlobalTailStart() == startPosMode){
      //std::cout << "found consensus tail " << std::endl;
      consensusTails_.push_back(t);
    }
  }
}

std::vector<aluHead> contigAlignment::getHeads(){
  return aluHeads_;
}

bool contigAlignment::isDenovo(){
  return isDenovo_;
}

bool contigAlignment::isTailDoubleStranded(){
  return tailDS_;
}

bool contigAlignment::isHeadDoubleStranded(){
  return headDS_;
}

int32_t contigAlignment::getReadsInRegion(){
  return readsInRegion_;
}

std::pair<int32_t, int32_t> contigAlignment::getGenotype(){
  return genotype_;
}

int32_t contigAlignment::getLongestTail(){
  return longestTail_;
}

int32_t contigAlignment::getMaxHash(){
  return maxHash_;
}

int32_t contigAlignment::getAltCount() {
  return contigAlignment::getConsensusTails().size() + contigAlignment::getHeads().size();
}

int32_t contigAlignment::getRefCount() {
  return readsInRegion_ - altCount_;
}

int32_t contigAlignment::getForwardStrandCount() {
  return forwardStrandCount_;
}

bool contigAlignment::isReadLeftBound(){
  return isLeftBound_;
}

std::vector<denovoEvidence> contigAlignment::getDenovoVec(){
  return denovoVec_;
}

clipCoords contigAlignment::getClipCoords(){
  return clipCoords_;
}

std::vector<polyA> contigAlignment::getConsensusTails(){
  return consensusTails_;
}

std::string contigAlignment::getCigarString(){
  return cigarString_;
}

int32_t contigAlignment::getAltForwardStrandCount(){
  return altForwardStrandCount_;
}

void contigAlignment::populateAltForwardStrandCount(){
  
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

    if(readsInRegion_ > 998){
      break;
    }

    polyA tail = {al, tailSize_};

    aluHead head = {util::getClipSeqs(alignedContig_)[0], al, headSize_};
    
    if(tail.detectPolyTail()){
      //std::cout << "detected polyATail" << std::endl;
      polyATails_.push_back(tail);
    }

    if(head.isHead()){
      //std::cout << "detect aluHead" << std::endl;
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
    denovoEvidence de = {util::getClipSeqs(alignedContig_)[0], alignedRegion_, pb, bamPath_, alignedContig_};
    denovoVec_.push_back(de);
    if(de.isDenovo()){
      isDenovo_ = true;
    }
  }
}

void contigAlignment::populateProbandGT(){

 

  std::vector<std::pair<std::string, int32_t> > refKmerCounts = util::countKmersFromText(probandRefPath_, refKmers_);
  std::vector<std::pair<std::string, int32_t> > altKmerCounts = util::countKmersFromText(probandAltPath_, altKmers_);

  probandGT_.RO = util::countKmerDepth(refKmerCounts);
  std::cout << "RO_ is: " << probandGT_.RO << std::endl;
  probandGT_.AO = util::countKmerDepth(altKmerCounts);
  std::cout << "AO_ is: " << probandGT_.AO << std::endl;
  probandGT_.RO = probandGT_.RO + probandGT_.AO;

  if(probandGT_.RO > 0  and probandGT_.AO > 0){
    probandGT_.genotype = std::make_pair(1, 0);
  }
  else if(probandGT_.RO > 0 and probandGT_.AO == 0){
   probandGT_.genotype = std::make_pair(0, 0);
  }
  else if(probandGT_.AO > 0 and probandGT_.RO == 0){
    probandGT_.genotype = std::make_pair(1, 1);
  }
  else {
    std::cerr << "Error in setGenotype, no ref or alt counts in genotype information" << std::endl;
    std::cerr << "Setting genotype to (0, 0) and proceeding" << std::endl;
    probandGT_.genotype = std::make_pair(0,0);
  }
  
}

genotypeField contigAlignment::getProbandGT(){
  return probandGT_;
}

//TODO:: Implement
void contigAlignment::populateLongestTail(){
  for(auto & t : polyATails_){
    if(t.getLongestTail() > longestTail_){
      longestTail_ = t.getLongestTail();
    }
  }
}

void contigAlignment::populateMaxHash(){
  std::vector<int32_t> peakVector = util::getPeakVector(alignedContig_);
  auto it = std::max_element(peakVector.begin(), peakVector.end());
  maxHash_ = peakVector[std::distance(peakVector.begin(), it)];
}

void contigAlignment::populateClipCoords(){
  clipCoords_ = util::intersectPeaksAndClips(util::getPeaks(alignedContig_), util::getLocalClipCoords(alignedContig_));
}

void contigAlignment::populateCigarString(){
    for(auto it = std::begin(alignedContig_.CigarData); it != std::end(alignedContig_.CigarData); it++){
      cigarString_ += it->Type;
      cigarString_ += std::to_string(it->Length);
    }
}

void contigAlignment::populateTailDS(){
  bool fs = false;
  bool rs = false;

  for(auto & ct : consensusTails_){
    if(ct.al_.IsReverseStrand()){
      rs = true;
    }
    else{
      fs = true;
    }
  }
  tailDS_ = rs && fs;
}

void contigAlignment::populateHeadDS(){
  bool fs = false;
  bool rs = false;

  for(auto & ah : aluHeads_){
    if(ah.al_.IsReverseStrand()){
      rs = true;
    }
    else{
      fs = true;
    }
  }
  headDS_ = rs && fs;
}

void contigAlignment::populateIsLeftBound(){
  isLeftBound_ = util::isReadLeftBound(alignedContig_.CigarData);
  std::cout << "isLeftBound is: " << isLeftBound_ << std::endl;
}


contigAlignment::contigAlignment(std::string bamPath, std::vector<std::string> parentBamPaths, std::pair<std::string, int32_t> aluHit, BamTools::BamAlignment alignedContig, std::string chrom, BamTools::BamRegion alignedRegion) : bamPath_(bamPath), parentBamPaths_(parentBamPaths), aluHit_(aluHit), alignedContig_(alignedContig), chrom_(chrom), alignedRegion_(alignedRegion){


  probandRefPath_ = bamPath_ + ".generator.V2.overlap.asembly.hash.fastq.Ref.sample";
  probandAltPath_ = bamPath_ + ".generator.V2.overlap.asembly.hash.fastq.sample";
  std::pair<int32_t, int32_t> breakpoint;

  if(util::isReadLeftBound(alignedContig_.CigarData)){
    breakpoint = std::make_pair(alignedContig_.RefID, alignedContig_.Position);
  }
  else{
    breakpoint = std::make_pair(alignedContig_.RefID, alignedContig_.GetEndPosition());
  }

  std::vector<BamTools::RefData> refData = util::populateRefData(bamPath);
  refSequence_ = util::pullRefSequenceFromRegion(breakpoint, referencePath_, refData, alignedContig_.QueryBases.size());
  refKmers_ = util::kmerize(refSequence_, 25);
  //altSequence_ = util::getClipSeqs(alignedContig_)[0];
  altSequence_ = alignedContig_.QueryBases;
  altKmers_ = util::kmerize(altSequence_, 25);


  contigAlignment::populateCigarString();
  contigAlignment::populateClipCoords();
  contigAlignment::populateIsLeftBound();
  contigAlignment::populateMaxHash();
  contigAlignment::populateHeadsAndTails();
  contigAlignment::populateConsensusTails();
  contigAlignment::populateTailDS();
  contigAlignment::populateHeadDS();
  contigAlignment::populateAltCount();
  contigAlignment::populateLongestTail();
  contigAlignment::populateProbandGT();
  contigAlignment::populateDenovoEvidence();
}

contigAlignment::~contigAlignment(){

}
