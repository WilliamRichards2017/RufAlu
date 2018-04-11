#ifndef INTERSECT_H
#define INTERSECT_H

#include <vector>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class Intersect{

 public:
  Intersect(const char *, const char *);
  ~Intersect();

  static const char * getContigHits(const char *);
  std::vector<BamTools::BamAlignment>  getIntersection();

 private:
  const char * a_;
  const char * b_;
  std::vector<BamTools::BamAlignment>  intersection_;
  void intersectBams();
  

};

#endif // INTERSECT_H
