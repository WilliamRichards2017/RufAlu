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
  return (coords_.clipDir == rtl);
}

bool polyA::isTailReverseStrand(){
  return al_.IsReverseStrand();
}



bool polyA::detectTailInWindow(clipCoords  c,  const char atChar){
  std::string seq = al_.QueryBases;

  //std::cout << "detecting tail for seq: " << seq << std::endl;
  //std::cout << "with coords: " << c.clipStart << ", " << c.clipEnd << ", " << c.clipDir << std::endl;
  
  if(c.clipDir == rtl){
    if(c.clipStart - tailSize_ < 0){
      return false;
    }
    
    int i = 0;
    while(i < tailSize_){
      const char cc = seq.at(c.clipStart-i);
      if(cc==atChar){
	  ++i;
      }
      else{
	return false;
      }
    }
    //std::cout << "found polyTail " << std::endl;
    
    return true;
  } else {
    if (c.clipStart + tailSize_ >= seq.size()-1){
      return false;
    }
    int i = 0;
    while(i < tailSize_){
      const char cc = seq.at(c.clipStart+i);
      if(cc==atChar){
	++i;
      }
      else{
	return false;
      }
    }
    //std::cout << "found polyTail " << std::endl;
    return true;
  }
  return false;
}

std::vector<clipCoords> polyA::getLocalClipCoords() {
  std::vector<clipCoords> coordsVec = {};
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;


  
  al_.GetSoftClips(clipSizes, readPositions, genomePositions);

  clipCoords c = {};
 
  for(int32_t i = 0; i < readPositions.size(); ++i){
        
    if(readPositions[i]-clipSizes[i]==0){
      c.clipDir = rtl;
      c.clipStart = readPositions[i]-1;
      c.clipEnd = c.clipStart - clipSizes[i]+1;
    }
    else{
      c.clipDir = ltr;
      c.clipStart = readPositions[i];
      c.clipEnd = readPositions[i] + clipSizes[i];
    }
    c.index = i;
    coordsVec.push_back(c);
    
    std::cout << "detecting tail for seq: " << al_.QueryBases << std::endl;
    std::cout << "with clip coords " << c.clipStart << ", " << c.clipEnd << ", " << c.clipDir << std::endl;
  }
  return coordsVec;
}


  
void polyA::setGlobalClipCoords(int32_t index){
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al_.GetSoftClips(clipSizes, readPositions, genomePositions);
  
  if(readPositions[index] - clipSizes[index] == 0){
    coords_.clipDir = rtl;
    coords_.clipStart = genomePositions[index] + clipSizes[index];
    coords_.clipEnd = genomePositions[index];
  }
  else{
    coords_.clipDir = ltr;
    coords_.clipStart = genomePositions[index];
    coords_.clipEnd = coords_.clipStart + clipSizes[index];
  }

  //std::cout << "Setting global clip coords to be: " << coords_.clipDir << ", " << coords_.clipStart << ", " << coords_.clipEnd << std::endl;
}

bool polyA::detectPolyTail(){
  //std::cout << "Inside detect polyA clips" << std::endl;
  std::vector<clipCoords > localClipCoords = polyA::getLocalClipCoords();
  //std::cout << "clips.size() is " << clips.size();
  std::string seq = al_.QueryBases;
  const char ac = 'A';
  const char tc = 'T';
  for(auto cIt = std::begin(localClipCoords); cIt != std::end(localClipCoords); ++cIt){
    //std::cout << "clip coords are:" << cIt->first << ", " << cIt->second << std::endl;
    //std::cout << "QueryBases is: " << seq << std::endl;
    bool a = polyA::detectTailInWindow(*cIt,  ac);
    bool t = polyA::detectTailInWindow(*cIt,  tc);
    if(a || t){
      polyA::setGlobalClipCoords(cIt->index);
      return true;
    }
  }
  return false;
}

polyA::polyA(BamTools::BamAlignment al, int32_t tailSize) : al_(al), tailSize_(tailSize){
  isTail_ = detectPolyTail();
}

polyA::~polyA(){
}
