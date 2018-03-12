#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "fastqParse.h"
#include "discordantAndChimericReads.h"

DACReads::DACReads(std::string filePath){      
  chimericReads_ = new std::vector<BamTools::BamAlignment*>;
  unmappedReads_ = new std::vector<BamTools::BamAlignment*>;
  DACReads::setAllChimericReads(filePath);
  DACReads::setAllUnmappedReads(filePath);
  
}

DACReads::~DACReads(){

  //  for(unsigned i = 0; i < chimericReads_->size(); ++i){
  //  delete &(chimericReads_)[i];
  // }

  //for(unsigned i = 0; i < unmappedReads_->size(); ++i){
  //  delete &(unmappedReads_)[i];
  // }
    
  chimericReads_->clear();
  unmappedReads_->clear();
}

void DACReads::setAllUnmappedReads(std::string inputFile){
  BamTools::BamReader reader;
  if(!reader.Open(inputFile)){
    std::cout << "Could not open input bam file " << inputFile << std::endl;
    return;
  }

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    if(!al.IsMapped()){
      std::cout << "found unmapped read" << std::endl;
      unmappedReads_->push_back(&al);
    }
  }

  return;
}

void DACReads::setAllChimericReads(std::string inputFile){
  BamTools::BamReader reader;
  if (!reader.Open(inputFile)){
    std::cout << "Could not open input Bam file" << inputFile << std::endl;
    return;
  }

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //Chimeric reads
    if(al.HasTag("SA")) {
      chimericReads_->push_back(&al);
      //std::cout << "read flag field is:" << al.AlignmentFlag << std::endl;
    }
  }
  return;
}

std::vector<BamTools::BamAlignment*> * DACReads::getAllChimericReads(){
  return chimericReads_;
}

std::vector<BamTools::BamAlignment*> * DACReads::getAllUnmappedReads(){
  return unmappedReads_;
}
