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

/*bool checkForClips(BamTools::BamAlignment al){

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
    }
  if(al.HasTag("SA")){
    std::cout << "found SA tag" << std::endl;
    return true;
  }
  return false;
}
*/

/*std::string Intersect::getContigHits(std::string overlapPath, std::string stub){
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
*/

void Intersect::intersectBams(){
  std::cout << "Inside intersectBams " << std::endl;
  std::vector<std::pair<BamTools::BamRegion, BamTools::BamAlignment> > coords;



  for(auto it = std::begin(contigVec_); it != std::end(contigVec_); ++it){
    std::pair<BamTools::BamRegion, BamTools::BamAlignment> regionPair;
    for(auto aIt = std::begin(it->contigAlignments); aIt != std::end(it->contigAlignments); ++aIt){
      BamTools::BamRegion region = BamTools::BamRegion(aIt->RefID, (aIt->Position)-5, aIt->RefID, (aIt->GetEndPosition())+5);
      std::cout << "pushing back contig alignment region" << std::endl;
      it->contigAlignmentRegions.push_back(region);
    }

  }

  BamTools::BamReader reader;

  if(!reader.Open(bamPath_)){
    std::cout << "could not open the following input bamfile: " << bamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  
  
  BamTools::BamAlignment bl;
  
  for(auto it = std::begin(contigVec_); it != std::end(contigVec_); ++it){
    for(auto aIt = std::begin(it->contigAlignmentRegions); aIt != std::end(it->contigAlignmentRegions); ++aIt){
      contigWindow cWindow{};
      reader.SetRegion(*aIt);
      
      while(reader.GetNextAlignment(bl)){
	//std::cout << "found overlaping read " << bl.Name << std::endl;
	cWindow.window.push_back(bl);
      }
      it->overlapingReads.push_back(cWindow);
    }
  }
}

std::vector<contig> Intersect::getContigVec() {
  return contigVec_;
}

Intersect::Intersect(std::vector<contig> contigVec, std::string  bamPath) : contigVec_(contigVec), bamPath_(bamPath) {
  Intersect::intersectBams();
}
