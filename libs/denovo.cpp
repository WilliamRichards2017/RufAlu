#include "aluHead.h"
#include "polyATail.h"
#include "denovo.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include <string>
#include <vector>

void denovoEvidence::findHeadsAndTails(){
  
  BamTools::BamReader reader;

  if (!reader.Open(parentBam_)){
    std::cerr << "Could not open input Bam file" << parentBam_ << std::endl;
    std::cerr << "Existing run with non-sero status.." << std::endl;
    exit (EXIT_FAILURE);
  }
  
  reader.LocateIndex();
  if (!reader.HasIndex()){
    std::cerr << "Index for" << parentBam_ << "could not be opened" << std::endl;
    std::cerr << "Exiting run with non-sero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }
  
  if(!reader.SetRegion(region_)){
    std::cerr << "Could not set region for coords: " << region_.LeftRefID << ", " << region_.LeftPosition << ", " << region_.RightRefID << ", " << region_.RightPosition << std::endl;
    std::cerr << "Exiting run with non-sero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }
  
  std::cout << "setting region coords to be: " << region_.LeftRefID << ", " << region_.LeftPosition << ", " << region_.RightRefID << ", " << region_.RightPosition << std::endl;
  
  
  BamTools::BamAlignment al;
  
  while(reader.GetNextAlignment(al)){
    if(!al.IsReverseStrand()){
      ++forwardCount_;
    }
    ++regionCoverage_;
    aluHead head = {aluClippedSeq_, al, 10};
    if (head.isHead()){
      std::cout << "found head in parent" << std::endl;
      parentHeads_.push_back(head);
    }
    
    polyA tail = {al, 10};
    
    if(tail.isTail()){
      std::cout << "found denovo tail in parent" << std::endl;
      parentTails_.push_back(tail);
    }
  }
  
}

int32_t denovoEvidence::getRegionCoverage() {
  return regionCoverage_;
}

int32_t denovoEvidence::getRefCount(){
  return regionCoverage_ - parentHeads_.size() - parentTails_.size();
}

int32_t denovoEvidence::getAltCount(){
  return parentHeads_.size() + parentTails_.size();
}

double denovoEvidence::getStrandBias(){
  return static_cast<float>(forwardCount_)/(static_cast<float>(regionCoverage_));
}

std::pair<int32_t, int32_t> denovoEvidence::getGenotype(){

  std::cout << "region coverage is: " << regionCoverage_ << std::endl;
  std::cout << "alt counts is: " << parentHeads_.size() + parentTails_.size() << std::endl;

  if(parentHeads_.size() + parentTails_.size() > 0 && regionCoverage_ > parentHeads_.size() + parentTails_.size()){
    return std::make_pair(1,1);
  }
  else if(parentHeads_.size() + parentTails_.size() < 1){
    return std::make_pair(1,0);
  }
  else if(regionCoverage_ <= parentHeads_.size() + parentTails_.size()){
    return std::make_pair(0,1);
  }
  else{
    return std::make_pair(-1,-1);
  }
}


denovoEvidence::denovoEvidence(std::string aluClippedSeq, BamTools::BamRegion region, std::string parentBam) : aluClippedSeq_(aluClippedSeq), region_(region), parentBam_(parentBam){
  denovoEvidence::findHeadsAndTails();
}


denovoEvidence::~denovoEvidence(){
}


bool denovoEvidence::isDenovo(){
  std::cout << "found " << parentHeads_.size() << " heads and " << parentTails_.size() << "tails" << std::endl;
  bool b = parentHeads_.size() > 0 || parentTails_.size() > 0;
  if(b){
    std::cout << "returning isDenovo" << std::endl;
  }

  return !b;
}


