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

struct contigAlignment;

class KnownAlus{
 public:
  KnownAlus(std::string, std::string, std::string, std::string, std::string, std::string, std::string);
  ~KnownAlus();
  
 private:
  std::string rawBamPath_;
  std::string contigFastqPath_;
  std::string contigBamPath_;
  std::string aluFastaPath_;
  std::string aluIndexPath_;
  std::string refPath_;
  std::string refIndexPath_;

  std::string stub_;
  std::vector<contig> contigVec_;
  std::vector<BamTools::RefData> refData_;
  std::string getChromosomeFromRefID(int32_t);

  void populateRefData();
  void findContigsContainingKnownAlus();
  void pullContigAlignments();
  void findReadsContainingPolyTails(int32_t);

  bool bedFilter(contigAlignment &);
  void writeBedPEHeader(std::ofstream &);
  void writeContigVecToBedPE(std::ofstream &);
  
  void printContigVec();

};

#endif // KNOWNALUS_H
