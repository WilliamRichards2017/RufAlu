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

std::vector<uint32_t> polyA::getClipStarts(BamTools::BamAlignment al){
  std::vector<BamTools::CigarOp> cigar = al.CigarData;
  std::vector<uint32_t> clipStarts;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  if (al.GetSoftClips(clipSizes, readPositions, genomePositions)){
    //for(auto it = std::begin(readPositions); it != std::end(readPositions); ++it){                                                                                                                   
    for(unsigned i = 0; i < readPositions.size(); ++i) {
      //std::cout << "Pushing back soft clip at pos " << readPositions[i] << std::endl;
      clipStarts.push_back(readPositions[i]);
    }
  }

  uint32_t pos = 0;
  for(auto it = std::begin(cigar); it != std::end(cigar); ++it){
    if (it->Type == 'H'){
      clipStarts.push_back(pos);
      //std::cout << "pushing back hard clip at pos " << pos << std::endl;
    }
    pos+= it->Length;
  }
  return clipStarts;
}
bool detectTailInWindow(BamTools::BamAlignment al, uint32_t pos, uint32_t tailSize, const char at){
  std::string seq = al.QueryBases;
  uint32_t i = 0;
  while(i < tailSize){
    const char c = seq[i+pos];
    if(c==at){
      ++i;
    }
    else{
      return false;
    }
  }
  return true;
}

bool polyA::detectPolyTailClips(BamTools::BamAlignment al, uint32_t tailSize){
  //std::cout << "Inside detect polyA clips" << std::endl;
  std::vector<uint32_t> clips = polyA::getClipStarts(al);
  //std::cout << "clips.size() is " << clips.size();
  std::string seq = al.QueryBases;
  const char ac = 'A';
  const char tc = 'T';
  for(auto it = std::begin(clips); it != std::end(clips); ++it){
    //std::cout << "looping through clips bb " << std::endl;
    if(*it+tailSize <= al.Length){
      //std::cout << "passing bounds check" << std::endl;
      bool a = detectTailInWindow(al, *it, tailSize, ac);
      bool t = detectTailInWindow(al, *it, tailSize, tc);
      if(a || t){
	//std::cout << "Found poly tail supporting evidence" << std::endl;
	return true;
      }
    }
  }
  return false;
}

uint32_t polyA::longestTail(BamTools::BamAlignment al){

  std::string seq = al.QueryBases;  
  uint32_t longestTail = 0;
  uint32_t savedLongest = 0;
  for(unsigned i = 0; i < seq.length(); ++i){
    const char* c = &seq[i];
    if(toupper(*c)=='T'){
      longestTail = 1;
      //for(unsigned j = i; j < seq.length(); ++j){
      while (i < seq.length() and toupper(seq[i])=='T'){
	++longestTail;
	++i;
      }
    } 
    else if(toupper(*c)=='A'){
      longestTail = 1;
      //for(unsigned j = i; j < seq.length(); ++j){                                                                                                                                                                                           
      while (i < seq.length() and toupper(seq[i])=='A'){
	++longestTail;
        ++i;
      } 
      
    }
    if(longestTail > savedLongest){
      savedLongest = longestTail;
      longestTail = 0;
    }
  }
  return savedLongest;
}
