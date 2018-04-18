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

std::pair<bool, int> polyA::detectPolyATail(std::string seq) {
  //std::cout << "detecting polyATail for sequence: " << seq << std::endl;
  std::list<const char*> window;
  int windowSize = 10;
  int pos = 0;
  int currentTail = 0;
  int maxTail = 0;
  float prop = 0.0;
  float targetProp = 0.95;
  bool b = false;
  if(seq.length() < windowSize ){
    std::pair<bool, int> r = std::make_pair(false,0);
    return r;
  }
  else{
    for(int i = 0; i < windowSize; ++i){
      const char* c = &seq[i];
      window.push_front(c);
      if(*c=='A'){
	prop+=1.0/(windowSize);
	++currentTail; 
      }
      else {
	currentTail = 0;
      }
      if(currentTail > maxTail){
	maxTail = currentTail;
      }
    }
  }
  
  if(prop >= targetProp){
    std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
    window.clear();
    b = true;
    std::pair<bool, int> ret = std::make_pair(b, maxTail);
    return ret;

  }
  
  for(unsigned i = windowSize; i < seq.length(); ++i){
    const char* c = &seq[i];
    window.push_front(c);
    if(toupper(*c)=='A' and toupper(*window.back()) != 'A'){
      prop+=1.0/windowSize;
      //printWindow(window);
      //std::cout << "Prop is: " << prop << std::endl;
      if(prop >= targetProp){
	std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
	window.clear();
	b = true;
	std::pair<bool, int> ret = std::make_pair(b, maxTail);
	return ret;

      }
      currentTail++;
    }
    else if(toupper(*c) != 'A' and toupper(*window.back())=='A'){
      //prop-=1.0/(window.size()-1);
      prop = std::max(0.0, prop-(1.0/(windowSize)));
      currentTail = 0;
    }
    else {
      currentTail=0;
    }
    if(currentTail > maxTail){
      maxTail = currentTail;
    }
    
    if(prop >= targetProp){
      std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
      b = true;
      std::pair<bool, int> ret = std::make_pair(b, maxTail);
      return ret;

    }   
    window.pop_back();
    }
  //std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;
  window.clear();
  std::cout << "About to return value of: " << b << " and " << maxTail << std::endl;
  std::pair<bool, int> ret = std::make_pair(b, maxTail);
  return ret;
}


std::pair<bool, int> polyA::detectPolyTTail(std::string seq) {
  //std::cout << "detecting polyATail for sequence: " << seq << std::endl;                                                                                                                                                                    
  std::list<const char*> window;
  int windowSize = 20;
  int pos = 0;
  int maxTail = 0;
  int currentTail = 0;
  bool b = false;
  float prop = 0.0;
  float targetProp = 0.95;

  if(seq.length() < windowSize ){
    std::pair<bool, int> r = std::make_pair(false,0);
    return r;
  }
  else{
    for(int i = 0; i < windowSize; ++i){
      const char* c = &seq[i];
      window.push_front(c);
      if(*c=='T'){
        prop+=1.0/(windowSize);
	currentTail++;
      }
      else{
	currentTail = 0;
      }
      if(currentTail > maxTail){
        maxTail = currentTail;
      }
    }
  }
  if(prop >= targetProp){
    std::cout << "Detected Poly T tail for sequence " << seq << " with score of " << prop << std::endl;
    window.clear();
    b = true;
    std::pair<bool, int> ret = std::make_pair(b, maxTail);
    return ret;

  }

  for(unsigned i = windowSize; i < seq.length(); ++i){
    const char* c = &seq[i];
    window.push_front(c);
    if(toupper(*c)=='T' and toupper(*window.back()) != 'T'){
      currentTail++;
      prop+=1.0/windowSize;
      //printWindow(window);                                                                                                                                                                                                                  
      //std::cout << "Prop is: " << prop << std::endl;                                                                                                                                                                                        
      if(prop >= targetProp){
	std::cout << "Detected Poly T tail for sequence " << seq << " with score of " << prop << std::endl;
        window.clear();
        b = true;
	std::pair<bool, int> ret = std::make_pair(b, maxTail);
	return ret;

      }
    }
    else if(toupper(*c) != 'T' and toupper(*window.back())=='T'){
      //prop-=1.0/(window.size()-1);                                                                                                                                                                                                          
      prop = std::max(0.0, prop-(1.0/(windowSize)));
      currentTail = 0;
    }
    else{
      currentTail = 0;
    }
    if(currentTail > maxTail){
      maxTail = currentTail;
    }
    
    if(prop >= targetProp){
      std::cout << "Detected Poly T tail for sequence " << seq << " with score of " << prop << std::endl;
      b = true;
      std::pair<bool, int> ret = std::make_pair(b, maxTail);
      return ret;

    }
    window.pop_back();
  }
  
  //std::cout<< "Failed to detected Poly a tail for sequence "<< seq << " with score of " << prop << std::endl;                                                                                                                               
  window.clear();
  std::pair<bool, int> ret = std::make_pair(b, maxTail);
  return ret;
}
