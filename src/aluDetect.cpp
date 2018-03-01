#include "fastqParse.h"
#include <stdio.h>


int main(void){

  std::string fastqPath = "../test_data/test.fastq";
  auto fastqSeqs = bioio::read_fastq_seqs(fastqPath);
  auto recordCount = bioio::count_fastq_records(fastqPath);

  std::cout << fastqSeqs[0];


  return 0;
  
}


bool detectPolyATail(std::string seq) {
  int pos = 0;
  int prop = 0.0;
  if(seq.size() < 10 ){
    return false;
  }
  else{
    char* window = new char[10]
    for(int i = 0; i < 10; ++1){
      char* c = seq[i];
      if(strcmp(seq[i], 'A')==0){
	prop+=0.1;
      }
      window[i] = c;
    }
  }
  
  for(unsigned i = 9; i < seq.size()-10; ++i){
   
  }
  if(prop >= 0.9){
    return true;
  }
  else{
    return false;
  } 
}
