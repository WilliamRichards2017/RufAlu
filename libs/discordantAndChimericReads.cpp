#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "fastqParse.h"
#include "discordantAndChimericReads.h"

DACReads::DACReads(std::string filePath){      
  chimericReads_ = new std::vector<BamTools::BamAlignment>;
  unmappedReads_ = new std::vector<BamTools::BamAlignment>;
  DACReads::setAllChimericReads(filePath);
  DACReads::setAllUnmappedReads(filePath);
  
}

DACReads::~DACReads(){

  delete[] chimericReads_;
  delete[] unmappedReads_;
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
      unmappedReads_->push_back(al);
    }
  }
  reader.Close();
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
      chimericReads_->push_back(al);
      //std::cout << "read flag field is:" << al.AlignmentFlag << std::endl;
    }
  }
  reader.Close();
  return;
}

std::vector<BamTools::BamAlignment> * DACReads::getAllChimericReads(){
  return chimericReads_;
}

std::vector<BamTools::BamAlignment> * DACReads::getAllUnmappedReads(){
  return unmappedReads_;
}
