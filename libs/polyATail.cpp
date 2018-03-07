#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>

#include "polyATail.h"


bool polyA::detectPolyATail(std::string seq) {
  std::list<const char*> window;
  int windowSize=10;
  int pos = 0;
  float prop = 0.0;
  if(seq.size() < 10 ){
    return false;
  }
  else{
    for(int i = 0; i < windowSize; ++i){
      const char* c = &seq[i];
      window.push_front(c);
      if(*c=='A'){
	prop+=0.1;
      }
    }
  }
  if(prop >= 0.9){
    std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
    return true;
  }
  
  for(unsigned i = 10; i < seq.size(); ++i){
    const char* c = &seq[i];
    window.push_front(c);

    if(*c=='A' and *window.back() != 'A'){
      prop+=1.0/(window.size()-1);
      if(prop >= 0.9){
	std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
	return true;
      }
    }
    else if(*c != 'A' and *window.back()=='A'){
      prop-=1.0/(window.size()-1);
    }
    if(prop >= 0.9){
      std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
      return true;
    }   
    window.pop_back();
    }
  std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;
  return false;
}


