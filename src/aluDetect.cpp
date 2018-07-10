#include <string>
#include <list>
#include <stdio.h>
#include <iostream>
#include <stdexcept>

#include "libs/knownAlus.h"

int main(int argc, char *argv[] ){
  
  if ( argc < 5 ){ 
    std::cout << "Please provide the proper number of arguments";
    exit (EXIT_FAILURE);
  }

  std::string bamPath = std::string(argv[1]);
  std::string contigFilePath = std::string(argv[2]);
  std::string contigBamPath = contigFilePath + ".bam";
  std::string aluFilePath = argv[3];
  std::string aluIndexPath = std::string(aluFilePath) + ".fai";
  std::string refPath = argv[4];
  std::string refIndexPath = std::string(refPath) + ".fai";
  std::string vcfOutPath = bamPath + ".V2.overlap.hashcount.fastq.bam.vcf";

  std::vector<std::string> parentBams;

  std::cout << "argc is: " << argc << std::endl;

  std::cout << "argv[5] is : " << argv[5] << std::endl;
  std::cout << "argv[6] is : " << argv[6] << std::endl;

  for(int i = 5; i < argc; ++i){
    parentBams.push_back(std::string(argv[i]));
  }

  for(auto p : parentBams){
    std::cout << "parent bam is: " << p << std::endl;
  }

  std::cout << "raw Bam path is " << bamPath << std::endl;
  std::cout << "Contig fastq path is: " << contigFilePath << std::endl;
  std::cout << "Contig bam path is: "<< contigBamPath << std::endl;
  std::cout << "Known Alu list path is: " << aluFilePath << std::endl;
  std::cout << "Reference path is: " << refPath << std::endl;
  std::cout << "Reference index path is: " << refIndexPath << std::endl;
  
  KnownAlus knownAlus = {bamPath, contigFilePath, contigBamPath, aluFilePath, aluIndexPath, refPath, refIndexPath, vcfOutPath, parentBams};
  
  std::cout << "Exiting Run." << std::endl;
  return 0;
}

