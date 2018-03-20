#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
//#include "kseq.h"
//#include "minimap.h"
#include "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/minimap2/src/minimap2_project/minimap.h"
#include <vector>
#include <utility> // std::pair


class KnownAlus{
 public:
  KnownAlus(const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
  std::vector<const char *> * getContigsContainingKnownAlus();
  

 private:
  const char * contigFilePath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  const char * refIndexPath_;

  std::vector<const char *> *contigsContainingKnownAlus_;
  std::vector<mm_reg1_t *> *knownAlus_;
  void findContigsContainingKnownAlus();
  void alignContigsContainingKnownAlus(const char *);
};

#endif // KNOWNALUS_H
