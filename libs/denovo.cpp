#include "aluHead.h"
#include "polyATail.h"
#include "denovo.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include <string>
#include <vector>

void denovoEvidence::findHeads(){
  for(auto pb : parentBams_){
    

    BamTools::BamReader reader;

    if (!reader.Open(pb)){
      std::cerr << "Could not open input Bam file" << pb << std::endl;
      std::cerr << "Existing run with non-sero status.." << std::endl;
      exit (EXIT_FAILURE);
    }

    reader.LocateIndex();
    if (!reader.HasIndex()){
      std::cerr << "Index for" << pb << "could not be opened" << std::endl;
      std::cerr << "Exiting run with non-sero status.." << std::endl;
      reader.Close();
      exit (EXIT_FAILURE);
    }

    if(!reader.SetRegion(region_)){
      std::cerr << "Could not set region for coords: " << region_.LeftRefID << ", " << region_.LeftPosition << ", " << region_.RightRefID << ", " << region_.RightPosition << std::endl;
      std::cerr << "Exiting run with non-sero status.." << std::endl;
      reader.Close();
      exit (EXIT_FAILURE);
    }
    
    std::cout << "setting region coords to be: " << region_.LeftRefID << ", " << region_.LeftPosition << ", " << region_.RightRefID << ", " << region_.RightPosition << std::endl;

    
    BamTools::BamAlignment al;

    while(reader.GetNextAlignment(al)){
      aluHead head = {aluClippedSeq_, al, 10};
      if (head.isHead()){
	std::cout << "found head in parent" << std::endl;
	parentHeads_.push_back(head);
      }

      polyA tail = {al, 10};

      if(tail.isTail()){
	std::cout << "found denovo tail in parent" << std::endl;
	parentTails_.push_back(tail);
      }
    }
    
  }
}

void denovoEvidence::findTails(){
}


denovoEvidence::denovoEvidence(std::string aluClippedSeq, BamTools::BamRegion region, std::vector<std::string> parentBams) : aluClippedSeq_(aluClippedSeq), region_(region), parentBams_(parentBams){
  denovoEvidence::findHeads();
  denovoEvidence::findTails();
  
}

denovoEvidence::~denovoEvidence(){
}

bool denovoEvidence::isDenovo(){
  std::cout << "found " << parentHeads_.size() << " heads and " << parentTails_.size() << "tails" << std::endl;
  bool b = parentHeads_.size() > 0 || parentTails_.size() > 0;
  if(b){
    std::cout << "returning isDenovo" << std::endl;
  }

  return !b;
}


