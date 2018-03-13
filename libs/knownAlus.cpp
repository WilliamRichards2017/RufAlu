#include "knownAlus.h"
#include "fastqParse.h"

#include <vector>

KnownAlus::KnownAlus(std::string contigFilePath, std::string aluFilePath) : contigFilePath_(contigFilePath), aluFilePath_(aluFilePath){
  KnownAlus::findContigsContainingKnownAlus();
}

KnownAlus::~KnownAlus(){
  contigsContainingKnownAlus_->clear();
}

void KnownAlus::findContigsContainingKnownAlus(){
  contigsContainingKnownAlus_ = new std::vector<std::string*>;
}

std::vector<std::string*> * KnownAlus::getContigsContainingKnownAlus(){
  return contigsContainingKnownAlus_;
}
