#include "fastqParse.h"
#include <string.h>
#include <ctype.h>
#include <list>
#include <stdio.h>

#include "polyATail.h"

//TODO: move to util
void printWindow(std::list<const char*> window){
  while(! window.empty()){
    std::cout << *(window.front());
    window.pop_front();
  }
  std::cout << std::endl;
}

std::vector<std::pair<int32_t, int32_t> > polyA::getClipCoords(BamTools::BamAlignment al){
  std::vector<BamTools::CigarOp> cigar = al.CigarData;
  
  std::vector<std::pair<int32_t, int32_t> > clipCoords;
  
  int32_t pos = 0;
  for(auto it = std::begin(cigar); it != std::end(cigar); ++it){
    if (it->Type == 'H' || it->Type == 'S'){
      std::pair<int32_t, int32_t> coordPair = std::make_pair(pos, pos+it->Length);
	clipCoords.push_back(coordPair);
    }
    pos+= it->Length;
  }  
  return clipCoords;

}


void printTailDebug(BamTools::BamAlignment  al, std::string tail){
  std::vector<BamTools::CigarOp> cigar = al.CigarData;
  std::cout << "Name of Read: " << al.Name << std::endl;
  std::cout << "Query Bases: " << al.QueryBases << std::endl;
  //std::cout << "Tail       : " << tail << std::endl;
  std::cout << "Cig        : ";
  for(auto cIt = std::begin(cigar); cIt != std::end(cigar); ++cIt){
    for(uint32_t i = 0;  i < cIt->Length; ++i){
      std::cout << cIt->Type;
    }
    
  }
  std::cout << std::endl;
  
  std::vector<std::pair<int32_t, int32_t> > clips = polyA::getClipCoords(al);
  std::cout << "Clipped coords are: ";
  for(auto cIt = std::begin(clips); cIt != std::end(clips); ++cIt){
    std::cout << cIt->first  << ", " << cIt->second << std::endl;
  }
  std::cout << std::endl;  
} 

std::string detectTailInWindowDebug(BamTools::BamAlignment al, uint32_t pos, uint32_t tailSize){

  std::string tail = "";
  for(uint32_t i = 0; i < pos; ++i){
    tail += "-";
  }
  for(uint32_t j = 0; j < tailSize; ++j){
    tail += "A";
  }
  return tail;
}



bool detectTailInWindow(BamTools::BamAlignment al, std::pair<int32_t, int32_t> coords, int32_t tailSize, const char at){
  std::string seq = al.QueryBases;

  if(coords.first == 0){
    int32_t i = 0;
    if(coords.second-tailSize < 0){
      //std::cout << "coords.second - tailSize is: " << coords.second-tailSize << std::endl;
      return false;
      while(i < tailSize){
	
	const char c = seq[coords.second-i];
	if(c==at){
	  ++i;
	}
	else{
	  return false;
	}
      }
      return true;
    }
    else{
      int32_t i = 0;
      //std::cout << "Coods.first + tailSize is: " << coords.first+tailSize << std::endl;
      if (coords.first+tailSize > al.Length-1){
	return false;
      }
      while(i < tailSize){
	const char c = seq[coords.first+i];
	if(c==at){
	  ++i;
	}
	else{
	  return false;
	}
      }
    }
    return true;
  }
  return false;
}
  
  
bool polyA::detectPolyTailClips(BamTools::BamAlignment al, int32_t tailSize){
  //std::cout << "Inside detect polyA clips" << std::endl;
  std::vector<std::pair<int32_t, int32_t> > clipCoords = polyA::getClipCoords(al);
  //std::cout << "clips.size() is " << clips.size();
  std::string seq = al.QueryBases;
  const char ac = 'A';
  const char tc = 'T';
  for(auto cIt = std::begin(clipCoords); cIt != std::end(clipCoords); ++cIt){
    bool a = detectTailInWindow(al, *cIt, tailSize, ac);
    bool t = detectTailInWindow(al, *cIt, tailSize, tc);
    
    if(a || t){
      /*std::string debugTail = "";
      debugTail = detectTailInWindowDebug(al, cIt->first, tailSize);
      printTailDebug(al, debugTail);*/
      return true;
    }
  }
  return false;
}
