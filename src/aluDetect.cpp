#include "fastqParse.h"
#include <string>
#include <list>
#include <stdio.h>
#include <iostream>
#include <stdexcept>


#include "libs/polyATail.h"
#include "libs/discordantAndChimericReads.h"
#include "libs/knownAlus.h"

int main(int argc, char *argv[] ){
  
  if ( argc != 4 ){ 
    std::cout << "Please provide the proper number of arguments";
    exit (EXIT_FAILURE);

  }

  const char * contigFilePath = argv[1];
  const char * aluFilePath = argv[2];
  const char * aluIndexPath = (std::string(aluFilePath) + ".fai").c_str();
  const char * refPath = argv[3];
  const char * refIndexPath = (std::string(refPath) + ".fai").c_str();;
  
  /*const char * bamPath = "../test_data/big-test.bam";
  const char * contigFilePath = "../../RUFUS/scripts/Family1.child.bam.generator.V2.overlap.hashcount.fastq";
  const char * aluFilePath = "../test_data/primate_non-LTR_Retrotransposon.fasta";
  const char * aluIndexPath =  "../test_data/primate_non-LTR_Retrotransposon.fasta.fai";
  const char * refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const char * refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";*/

  //DACReads *dacReads = new DACReads(bamPath);

  KnownAlus *knownAlus = new KnownAlus(contigFilePath, aluFilePath, aluIndexPath, refPath, refIndexPath);

  //odelete bamPath;
  //delete aluFilePath;
  //delete aluIndexPath;
  //delete refPath;
  //delete refIndexPath;
  

  return 0;
}

