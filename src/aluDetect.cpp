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
  
  if ( argc != 5 ){ 
    std::cout << "Please provide the proper number of arguments";
    exit (EXIT_FAILURE);
  }

  std::string bamPath = std::string(argv[1]);
  std::string contigFilePath = std::string(argv[2]);
  std::string contigBamPath = contigFilePath + ".bam";
  const char * aluFilePath = argv[3];
  const char * aluIndexPath = (std::string(aluFilePath) + ".fai").c_str();
  const char * refPath = argv[4];
  const char * refIndexPath = (std::string(refPath) + ".fai").c_str();

  std::cout << "raw Bam path is " << bamPath << std::endl;
  std::cout << "Contig fastq path is: " << contigFilePath << std::endl;
  std::cout << "Contig bam path is: "<< contigBamPath << std::endl;
  std::cout << "Known Alu list path is: " << aluFilePath << std::endl;
  std::cout << "Reference path is: " << refPath << std::endl;
  std::cout << "Reference index path is: " << refIndexPath << std::endl;
  
  KnownAlus knownAlus = {bamPath, contigFilePath, contigBamPath, aluFilePath, aluIndexPath, refPath, refIndexPath};
  return 0;
}

