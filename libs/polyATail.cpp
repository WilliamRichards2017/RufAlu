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


//TODO: implement
void polyA::setGlobalTailStart(){
  globalTailStart_ = globalClipCoords_.clipStart;
}


int32_t polyA::detectLeftTail(const clipCoords & cc, const char & c){
  //  std::cout << "Inside detect left tail " << std::endl;
  if (cc.clipStart == -1){
    return -1;
  }

  int32_t leftTail = 0;
  std::string leftClip = al_.QueryBases.substr(0, cc.clipStart);

  std::reverse(leftClip.begin(), leftClip.end());
  for(auto i : leftClip) {
    if(i==c){
      leftTail += 1;
    }
    else{
      break;
    }
  }
  return leftTail;
}

int32_t polyA::detectRightTail(const clipCoords & cc, const char & c){
  //std::cout << "inside detect left tail " << std::endl;
  if(cc.clipStart == -1){
    return -1;
  }

  int32_t rightTail = 0;
  std::string rightClip = al_.QueryBases.substr(cc.clipStart, cc.clipEnd-cc.clipStart+1);

  for(auto i : rightClip) {
    if(i==c){
      rightTail += 1;
    }
    else{
      break;
    }
  }
  return rightTail;
}


void polyA::setLongestTail(){
  int32_t mla = 0;
  int32_t mlt = 0;
  int32_t mra = 0;
  int32_t mrt = 0;

  int32_t longestTailIndex = 0;
  int32_t count = 0;
  
  std::pair<int32_t, int32_t> leftATail = std::make_pair(0, 0);
  std::pair<int32_t, int32_t> leftTTail = std::make_pair(0, 0);
  std::pair<int32_t, int32_t> rightATail = std::make_pair(0, 0);
  std::pair<int32_t, int32_t> rightTTail = std::make_pair(0, 0);


  
  
  for(auto c : localClipCoords_) {
        
    int32_t lat = 0;
    int32_t ltt = 0;
    int32_t rat = 0;
    int32_t rtt = 0;

    //std::cout << "Clip coords are: " << c.clipStart << ", " << c.clipEnd << ", " << c.clipDir << std::endl;
    
    if(c.clipDir == rtl){
      lat = polyA::detectLeftTail(c, 'A');
      ltt = polyA::detectLeftTail(c, 'T');
      //std::cout << "lat, ltt: " << lat << ", " << ltt << std::endl;
    }
    else if(c.clipDir == ltr){
      rat = polyA::detectRightTail(c, 'A');
      rtt = polyA::detectRightTail(c, 'T');
      //std::cout << "rat, rtt: " << rat << ", " << rtt << std::endl;
    }
    else {
      std::cout << "Clip dir is uninitialized " << std::endl;
      break;
    }

    if(lat > mla){
      mla = lat;
      leftATail.first = lat;
      leftATail.second = count;
    }
    if(ltt > mlt){
      mlt = ltt;
      leftTTail.first = ltt;
      leftTTail.second = count;
    }
    if(rat > mra){
      mra = rat;
      rightATail.first = rat;
      rightATail.second = count;
    }
    if(rtt > mrt){
      mrt = rtt;
      rightTTail.first = rtt;
      rightTTail.second = count;
    }
    ++count;
    //std::cout << "count is " << count << std::endl;
  }


  std::vector<std::pair<int32_t, int32_t> > allMaxTails;
  allMaxTails.push_back(leftATail);
  allMaxTails.push_back(leftTTail);
  allMaxTails.push_back(rightATail);
  allMaxTails.push_back(rightTTail);


  // yi would be proud
  auto find = std::max_element(begin(allMaxTails), end(allMaxTails), 
  			       [](const std::pair<int32_t, int32_t> & a, const std::pair<int32_t, int32_t> & b)
			       { 
				 return a.first < b.first;
			       }
			       );
  longestTail_ = find->first;
  longestTailIndex_ = find->second; 

  //longestTail_ = 0;	     
  //longestTailIndex_ = 0;
  
}


bool polyA::isTail(){
  return isTail_;
}

bool polyA::isTailReverseStrand(){
  return al_.IsReverseStrand();
}

clipCoords polyA::getGlobalClipCoords(){
  return globalClipCoords_;
}

int32_t polyA::getGlobalTailStart(){
  return globalClipCoords_.clipStart;
}
  

void polyA::setGlobalClipCoords(){
  //std::cout << "index inside setGlobalClipCoords is " << longestTailIndex_ << std::endl;
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  al_.GetSoftClips(clipSizes, readPositions, genomePositions);
  
  const std::vector<int32_t> insertionVec = util::getInsertionVec(al_);
  
  if(clipSizes.size() > 0){
    
    if(util::isReadLeftBound(al_.CigarData)){
      globalClipCoords_.clipDir = rtl;
    }
    else{
      globalClipCoords_.clipDir = ltr;
    }
    
    
    //std::cout << "size of clipSizeVec: " << clipSizes.size() << std::endl;
    //std::cout << "values of clipSizeVec:" << std::endl;
    //for(auto c : clipSizes){
    //  std::cout << c << ",";
    //}
    //std::cout << std::endl;
    //std::cout << "longestTailIndex_ is : " << longestTailIndex_ << std::endl;
    //std::cout << "insertion vec values:" << std::endl;
    //for (auto i : insertionVec){
    //  std::cout << i << ",";
    //}
    //std::cout << std::endl;
    
    globalClipCoords_.clipStart = genomePositions[longestTailIndex_] + insertionVec[longestTailIndex_];
    globalClipCoords_.clipEnd = globalClipCoords_.clipStart + clipSizes[longestTailIndex_] - 1;
    
    //std::cout << "Setting global clip coords to be: " << globalClipCoords_.clipDir << ", " << globalClipCoords_.clipStart << ", " << globalClipCoords_.clipEnd << std::endl;
  }
  else{
    globalClipCoords_.clipStart = 0;
    globalClipCoords_.clipEnd = 0;
    globalClipCoords_.clipDir = ltr;
  }
}

bool polyA::detectPolyTail(){
  //std::cout << "is longest tail of " << longestTail_ << " greater than tailSize of " << tailSize_ << std::endl;
  if(longestTail_ >= tailSize_){
    return true;
  }
  return false;
}

void polyA::setLocalClipCoords(){
  localClipCoords_ = util::getLocalClipCoords(al_);
}

const int32_t polyA::getLongestTail(){
  return longestTail_;
}

polyA::polyA(BamTools::BamAlignment al, int32_t tailSize) : al_(al), tailSize_(tailSize){
  polyA::setLocalClipCoords();
  polyA::setGlobalClipCoords();
  polyA::setLongestTail();
  isTail_ = detectPolyTail();
}

polyA::~polyA(){
}
