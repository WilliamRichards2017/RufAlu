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

#include "contig.h"
#include "intersect.h"
#include "util.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


Intersect::~Intersect(){
}

bool checkForClips(BamTools::BamAlignment al){

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  if (al.GetSoftClips(clipSizes, readPositions, genomePositions)){
    return true;
  }

  for (auto it = std::begin(al.CigarData); it != std::end(al.CigarData); ++it){
    if (it->Type == 'S' || it->Type == 'H'){
      return true;
    }
  }
  if(al.HasTag("SA")){
    return true;
  }
  return false;
}

void Intersect::intersectBams(){
  
  BamTools::BamReader reader;
  BamTools::BamAlignment al;
  if (!reader.Open(bamPath_)){
    std::cout << "Could not open input Bam file" << bamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }

  std::cout << "Does intersect reader have index? : " << reader.HasIndex() << std::endl;

  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    for(auto caIt = std::begin(cvIt->contigAlignments); caIt != std::end(cvIt->contigAlignments); ++caIt){
      
      const BamTools::BamRegion & region = BamTools::BamRegion(caIt->alignedContig.RefID, caIt->alignedContig.Position,
							       caIt->alignedContig.RefID, caIt->alignedContig.GetEndPosition());
      caIt->alignedRegion = region;
    }
  }
}

std::vector<contig> Intersect::getContigVec() {
  return contigVec_;
}

Intersect::Intersect(std::vector<contig> contigVec, std::string  bamPath) : contigVec_(contigVec), bamPath_(bamPath) {
  Intersect::intersectBams();
}
