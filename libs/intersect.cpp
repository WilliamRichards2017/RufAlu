#include "intersect.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>




Intersect::~Intersect(){
}

void Intersect::intersectBams(){
  std::cout << "Inside intersectBams " << std::endl;
  //std::vector<BamTools::BamAlignment> * intersection_ = new std::vector<BamTools::BamAlignment>;
  //std::vector<BamTools::BamAlignment>  intersection_;
  std::vector<BamTools::BamRegion> coords;

  BamTools::BamReader reader;
  if(!reader.Open(a_)){
    std::cout << "Could not open the following input bamfile: " << a_ << std::endl;
    exit (EXIT_FAILURE);

  }

  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    //std::cout << "looping through reads in file a" << std::endl;
    
    if(al.RefID != -1){
     //if(al.Name.compare("NODE_1348.bam.generator.V2_336_L191_D11:8:3::MH0")==0) {
    
    // std::cout << "left Ref Id: " << al.RefID << "   left position: " << al.Position << "   right Ref Id: " << al.RefID << "  right position: " << al.GetEndPosition() << std::endl;
      BamTools::BamRegion region = BamTools::BamRegion(al.RefID, al.Position, al.RefID, al.GetEndPosition()); 
      std::cout << "left Ref Id: " << al.RefID << "   left position: " << al.Position << "   right Ref Id: " << al.RefID << "  right position: " << al.GetEndPosition() << std::endl;
      coords.push_back(region);
    }
  }

  reader.Close();
  

  if(!reader.Open(b_)){
    std::cout << "could not open the following input bamfile: " << b_ << std::endl;
    exit (EXIT_FAILURE);
  }


  BamTools::BamAlignment bl;

  for(auto it = std::begin(coords); it != std::end(coords); ++it){
    //std::cout << "looping through coords: " << std::endl;
    reader.SetRegion(*it);

    while(reader.GetNextAlignment(bl)){
      //std::cout << "found overlapping read" << std::endl;
      intersection_.push_back(bl);
    }
  }
}



std::vector<BamTools::BamAlignment>  Intersect::getIntersection(){
  return intersection_;
}

Intersect::Intersect(const char * a, const char * b) : a_(a), b_(b) {
  std::cout << "input bam files for intersection are: "<< a_ << ", " << b_ << std::endl;
  Intersect::intersectBams();
}
