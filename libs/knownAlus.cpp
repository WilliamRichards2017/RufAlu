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
  std::vector<std::string> fastaSeqs = bioio::read_fasta_seqs("../test_data/primate_non-LTR_Retrotransposon.fasta");
  //auto recordCount = bioio::count_fasta_records("../test_data/primate_non-LTR_Retrotransposon.fasta");

  for(int i=0; i < fastaSeqs.size(); ++i){
    //if(polyA::detectPolyATail(std::string(fastaSeqs[i]))){
    bool b = polyA::detectPolyATail(fastaSeqs[i]);
    if(b){
      std::cout << "found alu w/poly a tail" << std::endl;
      contigsContainingKnownAlus_->push_back(&fastaSeqs[i]);

    }
  }
  
  std::cout << "found " << contigsContainingKnownAlus_->size() << " reads containing known alu sequences" << std::endl;

}

std::vector<std::string*> * KnownAlus::getContigsContainingKnownAlus(){
  return contigsContainingKnownAlus_;
}
