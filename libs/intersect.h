#ifndef INTERSECT_H
#define INTERSECT_H

#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "util.h"


//forward delcaration
struct contig;

class Intersect{

 public:
  Intersect(std::vector<contig>, std::string);
  ~Intersect();

  std::vector<contig> getContigVec();

    
 private:
  

  std::vector<contig> contigVec_;
  //std::string a_;
  std::string bamPath_;
    void intersectBams();
  

};

#endif // INTERSECT_H
