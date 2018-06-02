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
#include "intersect.h"
#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


class KnownAlus{
 public:
  KnownAlus(std::string, std::string, std::string, std::string, const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
 private:
  std::string rawReadsPath_;
  std::string contigFastqPath_;
  std::string contigBamPath_;
  std::string  mutationPath_;
  const char * aluPath_;
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
  void pullNamesWithHits(std::string);
  void findReadsContainingPolyTails(uint32_t);
  void writeBedPEHeader(std::ofstream &);
  void writeContigVecToBedPE(std::ofstream &);

};

#endif // KNOWNALUS_H
