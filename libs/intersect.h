#ifndef INTERSECT_H
#define INTERSECT_H

#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct contigWindow{
  BamTools::BamAlignment contig;
  std::vector<BamTools::BamAlignment> window;

};

class Intersect{

 public:
  Intersect(std::string, std::string);
  ~Intersect();

  static std::string getContigHits(std::string, std::string);
  std::vector<contigWindow>  getIntersection();

 private:
  std::string a_;
  std::string b_;
  std::vector<contigWindow>  intersection_;
  void intersectBams();
  

};

#endif // INTERSECT_H
