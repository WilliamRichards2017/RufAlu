#include "fastqParse.h"
#include <string.h>
#include <ctype.h>
#include <list>
#include <stdio.h>

#include "polyATail.h"


void printWindow(std::list<const char*> window){
  while(! window.empty()){
    std::cout << *(window.front());
    window.pop_front();
  }
  std::cout << std::endl;
}

std::vector<uint32_t> getClipPositions(BamTools::BamAlignment al){
  std::vector<BamTools::CigarOp> cigar = al.CigarData;

  std::vector<uint32_t> clips;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  if (al.GetSoftClips(clipSizes, readPositions, genomePositions)){
    for(auto it = std::begin(readPositions); it != std::end(readPositions); ++it){
      clips.push_back(*it);
    }
  }

  uint32_t pos = 0;
  for(auto it = std::begin(cigar); it != std::end(cigar); ++it){
    if (it->Type == 'H'){
      clips.push_back(pos);
    }
    pos+= it->Length;
  }
  return clips;
}

bool checkBounds(int32_t pos, std::vector<uint32_t> posVec, std::string readName){
  for(auto it = std::begin(posVec); it != std::end(posVec); ++it){
    if((pos > (*it)-11) and pos < (*it)+11){
      //std::cout << "AYYY check bounds passed for read: " << readName << std::endl;
      //std::cout << "Start of tail detected at position: " << pos << "clipping starts at: " << *it << std::endl;
      return true;
    }
  }
}

bool polyA::detectPolyATail(BamTools::BamAlignment al) {
  std::string seq = al.QueryBases;


  std::vector<uint32_t> clips = getClipPositions(al);

  //std::cout << "detecting polyATail for sequence: " << seq << std::endl;
  std::list<const char*> window;
  int windowSize = 10;
  int pos = 0;
  float prop = 0.0;
  float targetProp = 0.95;
  if(seq.length() < windowSize ){
    return false;
  }
  else{
    for(int i = 0; i < windowSize; ++i){
      const char* c = &seq[i];
      window.push_front(c);
      if(*c=='A'){
	prop+=1.0/(windowSize);
      }
    }
  }
  if(prop >= targetProp and checkBounds(10, clips, al.Name)){
    std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
    window.clear();
    return true;
  }
  
  for(unsigned i = windowSize; i < seq.length(); ++i){
    const char* c = &seq[i];
    window.push_front(c);
    if(toupper(*c)=='A' and toupper(*window.back()) != 'A'){
      prop+=1.0/windowSize;
      //printWindow(window);
      //std::cout << "Prop is: " << prop << std::endl;
      if(prop >= targetProp and checkBounds(i-windowSize, clips, al.Name)){
	std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
	window.clear();
	return true;
      }
    }
    else if(toupper(*c) != 'A' and toupper(*window.back())=='A'){
      //prop-=1.0/(window.size()-1);
      prop = std::max(0.0, prop-(1.0/(windowSize)));
    }
    if(prop >= targetProp and checkBounds(i-windowSize, clips, al.Name)){
      std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
      return true;
    }   
    window.pop_back();
  }
  //std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;
  window.clear();
  return false;
}


bool polyA::detectPolyTTail(BamTools::BamAlignment al) {
  std::string seq = al.QueryBases;
  //std::cout << "detecting polyATail for sequence: " << seq << std::endl;                                                                                                                                                                    
  std::vector<uint32_t> clips = getClipPositions(al);


  std::list<const char*> window;
  int windowSize = 20;
  int pos = 0;
  float prop = 0.0;
  float targetProp = 0.95;

  if(seq.length() < windowSize ){
    return false;
  }
  else{
    for(int i = 0; i < windowSize; ++i){
      const char* c = &seq[i];
      window.push_front(c);
      if(*c=='T'){
        prop+=1.0/(windowSize);
      }
    }
  }
  if(prop >= targetProp and checkBounds(0, clips, al.Name)){
    std::cout << "Detected Poly T tail for sequence " << seq << " with score of " << prop << std::endl;
    window.clear();
    return true;
  }

  for(unsigned i = windowSize; i < seq.length(); ++i){
    const char* c = &seq[i];
    window.push_front(c);
    if(toupper(*c)=='T' and toupper(*window.back()) != 'T'){
      prop+=1.0/windowSize;
      //printWindow(window);                                                                                                                                                                                                                  
      //std::cout << "Prop is: " << prop << std::endl;                                                                                                                                                                                        
      if(prop >= targetProp and checkBounds(i-windowSize, clips, al.Name)){
	std::cout << "Detected Poly t tail for sequence " << seq << " with score of " << prop << std::endl;
        window.clear();
        return true;
      }
    }
    else if(toupper(*c) != 'T' and toupper(*window.back())=='T'){
      //prop-=1.0/(window.size()-1);                                                                                                                                                                                                          
      prop = std::max(0.0, prop-(1.0/(windowSize)));
    }
    if(prop >= targetProp and checkBounds(i-windowSize, clips, al.Name)){
      std::cout << "Detected Poly T tail for sequence " << seq << " with score of " << prop << std::endl;
      return true;
    }
    window.pop_back();
  }
  //std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;                                                                                                                               
  window.clear();
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
