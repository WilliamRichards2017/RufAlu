#include "intersect.h"
#include "util.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>

Intersect::~Intersect(){
}

bool checkForClips(BamTools::BamAlignment al){

  //std::cout << "read name is: " << bl.Name << std::endl;                                                                                                                                                                              
  //std::vector<int> clipSizes;
  //std::vector<int> readPositions;
  //std::vector<int> genomePositions;
  /*if (al.GetSoftClips(clipSizes, readPositions, genomePositions)){
    //std::cout << "found some soft clips BB " << std::endl;
    return true;
  }

  for (auto it = std::begin(al.CigarData); it != std::end(al.CigarData); ++it){
    if (it->Type == 'C' || it->Type == 'H'){
      //std::cout << "Detected soft of hard clip of type: " << it->Type << std::endl;
      return true;
    }
    }*/
  if(al.HasTag("SA")){
    std::cout << "found SA tag" << std::endl;
    return true;
  }
  return false;
}

std::string Intersect::getContigHits(std::string overlapPath, std::string stub){
  BamTools::BamReader overlapReader;
  BamTools::BamReader aluReader;
  std::string out = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus-all-hits.sorted" + stub + ".bam";

  if (!overlapReader.Open(overlapPath)){
    std::cout << "Could not open input Bam overlap file" << overlapPath << std::endl;
    exit (EXIT_FAILURE);
  }

  const BamTools::SamHeader header = overlapReader.GetHeader();
  const BamTools::RefVector references = overlapReader.GetReferenceData();
  BamTools::BamWriter writer;

  if (!writer.Open(out, header, references)){
    std::cout << "could not open bam writer" << std::endl;
  }
  if (!aluReader.Open("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted" + stub + ".bam")){
    std::cout << "Could not open input Bam file contigs-with-alus.sorted" + stub + ".bam" << std::endl;
    exit (EXIT_FAILURE);
  }

  std::vector<std::string> hitNames;

  BamTools::BamAlignment al; 
  while(aluReader.GetNextAlignment(al)){
    hitNames.push_back(al.Name);
  }
  aluReader.Close();

  BamTools::BamAlignment bl;
  while(overlapReader.GetNextAlignment(bl)){
    for(auto it = std::begin(hitNames); it !=std::end(hitNames); ++it){
      if (bl.Name.compare(*it) == 0){
	//std::cout << "found hit for alu\n";
	//std::cout << "read name is: " << bl.Name << std::endl;
	if (checkForClips(bl) and bl.MapQuality > 0) {
	  writer.SaveAlignment(bl);  
	}
      }
    }
  }
  
  //sort and index bam file before returning path
  util::exec(("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/bamtools/src/bamtools_project/bin/bamtools sort -in /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus-all-hits.sorted" + stub + ".bam -out /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus-all-hits.sorted" + stub + ".bam").c_str());
  util::exec(("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/bamtools/src/bamtools_project/bin/bamtools index -in /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus-all-hits.sorted" + stub + ".bam").c_str());

  return out;
}

void Intersect::intersectBams(){
  std::cout << "Inside intersectBams " << std::endl;
  //std::vector<BamTools::BamAlignment> * intersection_ = new std::vector<BamTools::BamAlignment>;
  //std::vector<BamTools::BamAlignment>  intersection_;
  std::vector<std::pair<BamTools::BamRegion, BamTools::BamAlignment> > coords;

  BamTools::BamReader reader;
  if(!reader.Open(a_)){
    std::cout << "Could not open the following input bamfile: " << a_ << std::endl;
    exit (EXIT_FAILURE);

  }

  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    //std::cout << "looping through reads in file " << a_ << std::endl;
    //std::cout << "al.RefID: " << al.RefID << std::endl;
    
    //if(al.RefID != -1){
    //if(al.Name.compare("NODE_1348.bam.generator.V2_539_L424_D15:12:14::MH0")==0) {
    
    // std::cout << "left Ref Id: " << al.RefID << "   left position: " << al.Position << "   right Ref Id: " << al.RefID << "  right position: " << al.GetEndPosition() << std::endl;
    std::pair<BamTools::BamRegion, BamTools::BamAlignment> regionPair;
    BamTools::BamRegion region = BamTools::BamRegion(al.RefID, al.Position, al.RefID, al.GetEndPosition()+5); 
   
    regionPair.first = region;
    regionPair.second = al;
    //std::cout << "left Ref Id: " << al.RefID << "   left position: " << al.Position << "   right Ref Id: " << al.RefID << "  right position: " << al.GetEndPosition() << std::endl;
    coords.push_back(regionPair);
    //}
  }

  reader.Close();
  

  if(!reader.Open(b_)){
    std::cout << "could not open the following input bamfile: " << b_ << std::endl;
    exit (EXIT_FAILURE);
  }


  BamTools::BamAlignment bl;

  for(auto it = std::begin(coords); it != std::end(coords); ++it){
    contigWindow cWindow{};

    cWindow.contig = (*it).second;
    //std::cout << "looping through coords: " << std::endl;
    reader.SetRegion((*it).first);
    
    while(reader.GetNextAlignment(bl)){
      cWindow.window.push_back(bl);
      //std::cout << "found overlapping read" << std::endl;
      //intersection_.push_back(bl);
    }
    intersection_.push_back(cWindow);
    //std::cout << "pushing back cWindow" << std::endl;
    //util::printContigWindow(cWindow);
        
  }
}



std::vector<contigWindow>  Intersect::getIntersection(){
  return intersection_;
}

Intersect::Intersect(std::string a, std::string  b) : a_(a), b_(b) {
  std::cout << "input bam files for intersection are: "<< a_ << ", " << b_ << std::endl;
  Intersect::intersectBams();
}
