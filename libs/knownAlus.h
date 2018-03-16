#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
#include <vector>
#include <tuple>

using FastaRecord = std::tuple<std::string, std::string>;

class KnownAlus{
 public:
  KnownAlus(const char *, const char *, const char *);
  ~KnownAlus();
  
  std::vector<const char *> * getContigsContainingKnownAlus();
  

 private:
  const char * contigFilePath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  std::vector<const char *> *contigsContainingKnownAlus_;
  
  void findContigsContainingKnownAlus();
  void alignContigsContainingKnownAlus(const char *);
  void runMiniMap();

};

#endif // KNOWNALUS_H
