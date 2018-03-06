#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;

#include "fastqParse.h"
#include <string.h>
#include <list>
#include <stdio.h>

#include "polyATail.h"


bool aluDetect::detectPolyATail(std::string seq) {
  std::list<const char*> window;
  int pos = 0;
  float prop = 0.0;
  if(seq.size() < 10 ){
    return false;
  }
  else{
    for(int i = 0; i < 10; ++i){
      const char* c = &seq[i];
      std::cout << "checking if init character " << *c << " equals " << 'A' << std::endl;
      window.push_front(c);
      if(*c=='A'){
	prop+=0.1;
	std::cout << *c << "should be A" << std::endl;
	std::cout << "prop during initialization is: " << prop << std::endl;
      }
    }
  }

  if(prop >= 0.9){
    std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
    return true;
  }

  std::cout << "prop after initialization is: " << prop << std::endl;
  
  for(unsigned i = 10; i < seq.size(); ++i){
    const char* c = &seq[i];
    window.push_front(c);
    

    if(*c=='A' and *window.back() != 'A'){
      prop+=1.0/(window.size()-1);
      std::cout << "prop after addition is: " << prop << std::endl;
      if(prop >= 0.9){
	std::cout << "Detected Poly a tail for sequence " << seq << " with score of " << prop << std::endl;
	return true;
      }
    }
    else if(*c != 'A' and *window.back()=='A'){
      prop-=1.0/(window.size()-1);
      std::cout << "prop after subtraction is: " << prop << std::endl;
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


