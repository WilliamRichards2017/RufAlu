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

bool polyA::detectPolyATail(std::string seq) {
  //std::cout << "detecting polyATail for sequence: " << seq << std::endl;
  std::list<const char*> window;
  int windowSize = 30;
  int pos = 0;
  float prop = 0.0;
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
  if(prop >= 0.9){
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
      if(prop >= 0.90){
	std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
	window.clear();
	return true;
      }
    }
    else if(toupper(*c) != 'A' and toupper(*window.back())=='A'){
      //prop-=1.0/(window.size()-1);
      prop = std::max(0.0, prop-(1.0/(windowSize)));
    }
    if(prop >= 0.9){
      std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
      return true;
    }   
    window.pop_back();
    }
  //std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;
  window.clear();
  return false;
}


