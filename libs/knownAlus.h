#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
#include "kseq.h"
#include "minimap.h"
#include <vector>
#include <string>
#include <iostream>
#include <utility> // std::pair

#include "contig.h"
#include "intersect.h" //contigWindow struct
#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


//Forward decleration of required struct
struct contig;

class KnownAlus{
 public:
  KnownAlus(std::string std::string, std::string, std::string, const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
 private:
  std::string bamPath_;
  std::string contigFilePath_;
  std::string contigBamPath_;
  std::string  mutationPath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  const char * refPath_;
  const char * refIndexPath_;

  std::string stub_;

  std::vector<contig> contigVec_;
  std::vector<BamTools::RefData> refData_;

  std::string getChromosomeFromRefID(int32_t);
  void populateRefData(std::string);
  void findContigsContainingKnownAlus();
  void mapContigsToRef(const char *);
  void pullNamesWithHits(std::vector<contig>, std::string);
  void findReadsContainingPolyTails(std::string, uint32_t);
  void writeBedPEHeader(std::ofstream &)
  void writeContigVecToBedPE(std::ofstream &);

};

#endif // KNOWNALUS_H
