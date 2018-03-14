#include "knownAlus.h"
#include "fastqParse.h"
#include "polyATail.h"


#include <vector>

KnownAlus::KnownAlus(std::string contigFilePath, std::string aluFilePath) : contigFilePath_(contigFilePath), aluFilePath_(aluFilePath){
  KnownAlus::findContigsContainingKnownAlus();
}

KnownAlus::~KnownAlus(){
  contigsContainingKnownAlus_->clear();
}

void KnownAlus::findContigsContainingKnownAlus(){
  contigsContainingKnownAlus_ = new std::vector<std::string*>;
  auto fastaSeqs = bioio::read_fasta_seqs("../test_data/primate_non-LTR_Retrotransposon.fasta");
  //auto recordCount = bioio::count_fasta_records("../test_data/primate_non-LTR_Retrotransposon.fasta");

  for(int i=0; i < fastaSeqs.size(); ++i){
    contigsContainingKnownAlus_->push_back(&fastaSeqs[i]);
  }

  std::cout << "found " << contigsContainingKnownAlus_->size() << " reads containing known alu sequences" << std::endl;

}

std::vector<std::string*> * KnownAlus::getContigsContainingKnownAlus(){
  return contigsContainingKnownAlus_;
}
