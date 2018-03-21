#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>
#include "libs/polyATail.h"
#include "libs/discordantAndChimericReads.h"
#include "libs/knownAlus.h"

int main(void){

  const char * bamPath = "../test_data/big-test.bam";
  //std::string bamPath = "/uufs/chpc.utah.edu/common/home/u0991464/RUFUS.test.set/Family1.mother.bam";
  const char * fastqPath = "../test_data/test.fastq";
  auto fastqSeqs = bioio::read_fastq_seqs(fastqPath);
  auto recordCount = bioio::count_fastq_records(fastqPath);



  
  const char * contigFilePath = "../../RUFUS/scripts/Family1.child.bam.generator.V2.overlap.hashcount.fastq";

  // const char * contigFilePath = "../test_data/primate_non-LTR_Retrotransposon.fasta";

  const char * aluFilePath = "../test_data/primate_non-LTR_Retrotransposon.fasta";
  const char * aluIndexPath =  "../test_data/primate_non-LTR_Retrotransposon.fasta.fai";
  const char * refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const char * refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";
  

  //const char * contigFilePath = "../test_data/primate_non-LTR_Retrotransposon.fasta";
  //const char * aluFilePath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  //const char * aluIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";

  DACReads *dacReads = new DACReads(bamPath);
  KnownAlus *knownAlus = new KnownAlus(contigFilePath, aluFilePath, aluIndexPath, refPath, refIndexPath);

  std::cout << "finished building DAC reads\n";
 
  delete dacReads;
  
  std::cout << "destroyed DAC reads\n";
  
  //  for(int i=0; i < fastqSeqs.size(); ++i){
  //  polyA::detectPolyATail(fastqSeqs[i]);
  //}
  
  return 0;
}

