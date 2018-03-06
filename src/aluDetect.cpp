#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>
#include "../libs/polyATail.h"

int main(void){

  std::string fastqPath = "../test_data/test.fastq";
  auto fastqSeqs = bioio::read_fastq_seqs(fastqPath);
  auto recordCount = bioio::count_fastq_records(fastqPath);

  std::cout << fastqSeqs[0];

  for(int i=0; i < fastqSeqs.size(); ++i){
    aluDetect::detectPolyATail(fastqSeqs[i]);
  }

  return 0;
}

