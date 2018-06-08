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


/*void printTailDebug(BamTools::BamAlignment al){
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
  
  std::vector<std::pair<int32_t, int32_t> > clips = polyA::getClipCoords();
  std::cout << "Clipped coords are: ";
  for(auto cIt = std::begin(clips); cIt != std::end(clips); ++cIt){
    std::cout << cIt->first  << ", " << cIt->second << std::endl;
  }
  std::cout << std::endl;  
  }
*/ 


bool polyA::isTail(){
  return isTail_;
}

bool polyA::isTailLeftBound(){
  return isTailLeftBound_;
}

bool polyA::isTailReverseStrand(){
  return isTailReverseStrand_;
}

bool polyA::detectTailInWindow(std::pair<int32_t, int32_t>  coords,  const char at){
  std::string seq = al_.QueryBases;
  
  if(coords.first == 0){
    isTailLeftBound_ = true;

    int32_t i = 0;
    
    if(coords.second-tailSize_ < 0){
      //std::cout << "coords.second - tailSize is: " << coords.second-tailSize << std::endl;
      return false;
    }
    while(i < tailSize_){
      
      const char c = seq[coords.second-i];
      if(c==at){
	  ++i;
      }
      else{
	return false;
      }
    }
    return true;
  } else {
    int32_t i = 0;
    if (coords.first+tailSize_ > al_.Length){
      return false;
    }
    while(i < tailSize_){
      const char c = seq[coords.first+i];
      if(c==at){
	++i;
      }
      else{
	return false;
      }
      
    }
    return true;
  }
  return false;
}

std::vector<std::pair<int32_t, int32_t> > polyA::getLocalClipCoords(){
  std::vector<BamTools::CigarOp> cigar = al_.CigarData;

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
  
void polyA::setGlobalClipCoords(std::pair<int32_t, int32_t>  c){
  if(c.first == 0){
    coords_.clipStart = c.second + al_.Position;
    coords_.clipEnd = c.first + al_.Position;
    coords_.clipDir = rtl;  
  }
  else{
    coords_.clipStart = c.first + al_.Position;
    coords_.clipEnd - c.second + al_.Position;
    coords_.clipDir = ltr;
  }
}
  
bool polyA::detectPolyTail(){
  //std::cout << "Inside detect polyA clips" << std::endl;
  std::vector<std::pair<int32_t, int32_t> > localClipCoords = polyA::getLocalClipCoords();
  //std::cout << "clips.size() is " << clips.size();
  std::string seq = al_.QueryBases;
  const char ac = 'A';
  const char tc = 'T';
  for(auto cIt = std::begin(localClipCoords); cIt != std::end(localClipCoords); ++cIt){
    bool a = polyA::detectTailInWindow(*cIt,  ac);
    bool t = polyA::detectTailInWindow(*cIt,  tc);
    if(a || t){
      polyA::setGlobalClipCoords(*cIt);
      return true;
    }
  }
  return false;
}

polyA::polyA(BamTools::BamAlignment al, int32_t tailSize) : al_(al), tailSize_(tailSize){
  isTail_ = detectPolyTail();
  if(al_.IsReverseStrand()){
    isTailReverseStrand_ = true;
  }
}

polyA::~polyA(){
}
