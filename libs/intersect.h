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
  Intersect(const char *, const char *);
  ~Intersect();

  static const char * getContigHits(const char *);
  std::vector<contigWindow>  getIntersection();

 private:
  const char * a_;
  const char * b_;
  std::vector<contigWindow>  intersection_;
  void intersectBams();
  

};

#endif // INTERSECT_H
