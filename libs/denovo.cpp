#include "aluHead.h"
#include "polyATail.h"
#include "denovo.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "util.h"

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
  
  //std::cout << "setting region coords to be: " << region_.LeftRefID << ", " << region_.LeftPosition << ", " << region_.RightRefID << ", " << region_.RightPosition << std::endl;
  
  
  BamTools::BamAlignment al;
  
  while(reader.GetNextAlignment(al)){
    if(!al.IsReverseStrand()){
      ++forwardCount_;
    }
    ++regionCoverage_;
    aluHead head = {aluClippedSeq_, al, 10};
    if (head.isHead()){
      parentHeads_.push_back(head);
    }
    
    polyA tail = {al, 10};
    
    if(tail.isTail()){
      parentTails_.push_back(tail);
    }
  }
  
}

int32_t denovoEvidence::getRegionCoverage() {
  return regionCoverage_;
}

int32_t denovoEvidence::getRefCount(){
  return RO_;
}

int32_t denovoEvidence::getAltCount(){
  return AO_;
}

double denovoEvidence::getStrandBias(){
  return static_cast<float>(forwardCount_)/(static_cast<float>(regionCoverage_));
}

void denovoEvidence::populateGTFields(){
  std::vector<std::pair<std::string, int32_t> > refKmerCounts = util::countKmersFromText(refPath_, refKmers_);
  std::vector<std::pair<std::string, int32_t> > altKmerCounts = util::countKmersFromText(altPath_, altKmers_);

  RO_ = util::countKmerDepth(refKmerCounts);
  std::cout << "RO_ is: " << RO_ << std::endl;
  AO_ = util::countKmerDepth(altKmerCounts);
  std::cout << "AO_ is: " << AO_ << std::endl;
  DP_ = RO_ + AO_;
  
  if(RO_ > 0  and AO_ > 0){
    genotype_ = std::make_pair(1, 0);
  }
  else if(RO_ > 0 and AO_ == 0){
    genotype_ = std::make_pair(0, 0);
  }
  else if(AO_ > 0 and RO_ == 0){
    genotype_ = std::make_pair(1, 1);
  }
  else {
    std::cerr << "Error in setGenotype, no ref or alt counts in genotype information" << std::endl;
    std::cerr << "Setting genotype to (0, 0) and proceeding" << std::endl;
    genotype_ = std::make_pair(0,0);
  }
}

std::pair<int32_t, int32_t> denovoEvidence::getGenotype(){
  return genotype_;
}


denovoEvidence::denovoEvidence(const std::string & aluClippedSeq, const BamTools::BamRegion & region, const std::string & parentBam, const std::string & probandBam,  const BamTools::BamAlignment & al) : aluClippedSeq_(aluClippedSeq), region_(region), parentBam_(parentBam), probandBam_(probandBam), al_(al){
  
  refPath_ = probandBam + ".generator.V2.overlap.asembly.hash.fastq.Ref." + util::baseName(parentBam_) + ".generator.Jhash";

  //Child.bam.generator.V2.overlap.asembly.hash.fastq.Ref.Father.bam.generator.Jhash
  
  altPath_ = probandBam + ".generator.V2.overlap.asembly.hash.fastq." + util::baseName(parentBam_) + ".generator.Jhash";
    
  std::vector<BamTools::RefData> refData = util::populateRefData(parentBam);
  //TODO: unhardcode ref path

  std::pair<int32_t, int32_t> breakpoint;

  if(util::isReadLeftBound(al.CigarData)){
    breakpoint = std::make_pair(al.RefID, al.Position);
  }
  else{
    breakpoint = std::make_pair(al.RefID, al.GetEndPosition());
  }

  refSequence_ = util::pullRefSequenceFromRegion(breakpoint, referencePath_, refData, al_.QueryBases.size());
  std::cout << "refSequence is: " << refSequence_ << std::endl;
  refKmers_ = util::kmerize(refSequence_, 25);
  std::cout << "refKmers_.size() is: " << refKmers_.size() << std::endl;
  altSequence_ = aluClippedSeq_;
  std::cout << "altSequence is: " << altSequence_ << std::endl;
  altKmers_ = util::kmerize(altSequence_, 25);
  std::cout << "altKmers_.size() is: " << altKmers_.size() << std::endl;
  

  denovoEvidence::findHeadsAndTails();
  denovoEvidence::populateGTFields();
}


denovoEvidence::~denovoEvidence(){
}


bool denovoEvidence::isDenovo(){
  //std::cout << "found " << parentHeads_.size() << " heads and " << parentTails_.size() << "tails" << std::endl;
  bool b = parentHeads_.size() < 1 and parentTails_.size() <  1;
  return b;
}


