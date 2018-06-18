#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <list>
#include <stdio.h>

#include "polyATail.h"
#include "util.h"

//TODO: move to util
void printWindow(std::list<const char*> window){
  while(! window.empty()){
    std::cout << *(window.front());
    window.pop_front();
  }
  std::cout << std::endl;
}

int32_t polyA::detectLeftTail(const clipCoords & cc, const char & c){
  if (cc.clipStart == -1){
    return -1;
  }
  //std::cerr << "Inside detectLeftTail, substr coords are: " <<cc.clipStart <<", " << cc.clipEnd <<  "for read: " 
  //    << al_.Name << " at coords " << al_.RefID << ", " << al_.Position  <<  std::endl;
  
  
  int32_t leftTail = 0;
  std::string leftClip = al_.QueryBases.substr(0, cc.clipStart);

  

  std::reverse(leftClip.begin(), leftClip.end());
  
  //std::cerr << "left clip string is: " << leftClip << std::endl;
  //std::cerr << "clip starts at " << cc.clipStart << std::endl;
  
  for(auto i : leftClip) {
    //std::cout << "checking if " << i << " == " << c << std::endl;
    if(i==c){
      leftTail += 1;
    }
    else{
      break;
    }
  }
  //std::cerr << "left Tail inside detectLeftTail is: " << leftTail << std::endl;
  return leftTail;
}

int32_t polyA::detectRightTail(const clipCoords & cc, const char & c){


  if(cc.clipStart == -1){
    return -1;
  }

  int32_t rightTail = 0;

  //std::cout << "Inside detectRightTail, substr coords are: " << cc.clipStart << ", " << cc.clipStart + (cc.clipEnd-cc.clipStart) <<  " for read: " << al_.Name << " at coords " << al_.RefID << ", " << al_.Position << std::endl;


  //std::cout << "query bases size is: " << al_.QueryBases.size() << std::endl;
  
  //std::cout << "Query bases are: " << al_.QueryBases << std::endl;
  std::string rightClip = al_.QueryBases.substr(cc.clipStart, cc.clipEnd-cc.clipStart);
  //std::cout << "rightClip is: " << rightClip << std::endl;

  for(auto i : rightClip) {
    //std::cout << "checking if " << i << " == " << c << std::endl;
    if(i==c){
      rightTail += 1;
    }
    else{
      break;
    }
  }
  return rightTail;
}



void polyA::setLongestTail(const std::vector<clipCoords> & clipVec){
  int32_t mla = -1;
  int32_t mlt = -1;
  int32_t mra = -1;
  int32_t mrt = -1;
  
  for(auto c : clipVec) {
    int32_t leftATail = 0;
    int32_t leftTTail = 0;
    int32_t rightATail = 0;
    int32_t rightTTail = 0;
    if(c.clipDir == rtl){
      leftATail = polyA::detectLeftTail(c, 'A');
      leftTTail = polyA::detectLeftTail(c, 'T');
      //std::string leftClip = al_.QueryBases.substr(0, c.clipStart-1);
    }
    else if(c.clipDir == ltr){
      rightATail = polyA::detectRightTail(c, 'A');
      rightTTail = polyA::detectRightTail(c, 'T');
    }
    
    if(leftATail > mla){
      mla = leftATail;
    }
    if(leftTTail > mlt){
      mlt = leftTTail;
    }
    if(rightATail > mra){
      mra = rightATail;
    }
    if(rightTTail > mrt){
      mrt = rightTTail;
    }
  }
  longestTail_ = std::max(mla, std::max(mlt, std::max(mra, mrt)));
}


bool polyA::isTail(){
  return isTail_;
}

bool polyA::isTailLeftBound(){
  return (coords_.clipDir == rtl);
}

bool polyA::isTailReverseStrand(){
  return al_.IsReverseStrand();
}

bool polyA::detectTailInWindow(const clipCoords &  c,  const char & atChar){
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

  const std::vector<int32_t> insertionVec = util::getInsertionVec(al_);

 
  for(int32_t i = 0; i < readPositions.size(); ++i){
        
    if(readPositions[i]-clipSizes[i]==0){
      c.clipDir = rtl;
      c.clipStart = readPositions[i] + insertionVec[i];
      c.clipEnd = 0;
    }
    else{
      c.clipDir = ltr;
      c.clipStart = readPositions[i] + insertionVec[i];
      c.clipEnd = c.clipStart + clipSizes[i];
    }
    c.index = i;
    coordsVec.push_back(c);
    
    //std::cout << "detecting tail for seq: " << al_.QueryBases << std::endl;
    //std::cout << "with clip coords " << c.clipStart << ", " << c.clipEnd << ", " << c.clipDir << std::endl;
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

const int32_t polyA::getLongestTail(){
  return longestTail_;
}

polyA::polyA(BamTools::BamAlignment al, int32_t tailSize) : al_(al), tailSize_(tailSize){
  isTail_ = detectPolyTail();
  polyA::setLongestTail(polyA::getLocalClipCoords());
}

polyA::~polyA(){
}
