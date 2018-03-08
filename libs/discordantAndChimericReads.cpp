#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "fastqParse.h"
#include "discordantAndChimericReads.h"

DACReads::DACReads(std::string filePath){
  DACReads::setAllChimericReads(filePath);
  DACReads::setAllUnmappedReads(filePath);
  
}

DACReads::~DACReads(){
  delete &chimericReads_;
  delete &unmappedReads_;
}

void DACReads::setAllUnmappedReads(std::string inputFile){
  return;
}

void DACReads::setAllChimericReads(std::string inputFile){
  std::cout << "inside get all chimeric reads" << std::endl;

  BamTools::BamReader reader;
  if (!reader.Open(inputFile)){
    std::cout << "Could not open input Bam file" << std::endl;
    return;
  }

  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //Chimeric reads
    if(al.HasTag("SA")) {
      chimericReads_.push_back(al);
      std::cout << "read flag field is:" << al.AlignmentFlag << std::endl;
    }
  }
}

std::vector<BamTools::BamAlignment> DACReads::getAllChimericReads(){
  return chimericReads_;
}

std::vector<BamTools::BamAlignment> DACReads::getAllUnmappedReads(){
  return unmappedReads_;
}
