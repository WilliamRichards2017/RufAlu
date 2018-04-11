#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
#include "kseq.h"
#include "minimap.h"
#include <vector>
#include <string>
#include <utility> // std::pair

#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class KnownAlus{
 public:
  KnownAlus(const char *, const char *, const char *, const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
  std::vector<fastqRead> * getContigsContainingKnownAlus();
  

 private:
  const char * contigFilePath_;
  const char * contigBamPath_;
  const char * mutationPath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  const char * refPath_;
  const char * refIndexPath_;

  std::vector<fastqRead> *contigsContainingKnownAlus_;
  const std::vector<BamTools::RefData> *refData_;

  void populateRefData(const char *);
  void findContigsContainingKnownAlus();
  void alignContigsContainingKnownAlus(const char *);
  void mapContigsToRef(const char *);
  void findContigsContainingPolyATails(const char *);
  void findReadsContainingPolyATails(std::vector<BamTools::BamAlignment>, const char *);
  std::string getChromosomeFromRefID(int32_t);

};

#endif // KNOWNALUS_H
