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
  if (!reader.Open(bamPath_)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }

  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    for(auto caIt = std::begin(cvIt->contigAlignments); cvIt != std::end(cvIt->contigAlignments); ++caIt){
      BamTools::BamAlignment al;
      BamTools::BamRegion region = BamTools::BamRegion(caIt->RefID, caIt->Position, caIt->RefID, caIt->GetEndPosition());
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
