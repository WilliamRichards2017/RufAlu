#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>
#include "../libs/polyATail.h"
#include "../libs/discordantAndChimericReads.h"

int main(void){

  //std::string bamPath = "../test_data/big-test.bam";
  std::string bamPath = "/uufs/chpc.utah.edu/common/home/u0991464/RUFUS.test.set/Family1.mother.bam";
  std::string fastqPath = "../test_data/test.fastq";
  auto fastqSeqs = bioio::read_fastq_seqs(fastqPath);
  auto recordCount = bioio::count_fastq_records(fastqPath);

  std::cout << fastqSeqs[0];

  DACReads dacReads(bamPath);

  for(int i=0; i < fastqSeqs.size(); ++i){
    polyA::detectPolyATail(fastqSeqs[i]);
  }

  return 0;
}

