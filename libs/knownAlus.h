#ifndef KNOWNALUS_H
#define KNOWNALUS_H


#include <vector>
#include <string>
#include <iostream>
#include <utility> // std::pair

#include "contig.h"
#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "kseq.h"
#include "minimap.h"


struct contig;
struct contigAlignment;

class KnownAlus{
 public:
  KnownAlus(std::string, std::string, std::string, std::string, std::string, std::string, std::string);
  ~KnownAlus();
  
 private:
  const std::string rawBamPath_;
  const std::string contigFastqPath_;
  const std::string contigBamPath_;
  const std::string aluFastaPath_;
  const std::string aluIndexPath_;
  const std::string refPath_;
  const std::string refIndexPath_;
  const std::string stub_;
  const std::string prefix_ = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/out/";

  std::vector<contig> contigVec_;
  std::vector<BamTools::RefData> refData_;
  std::string getChromosomeFromRefID(int32_t);

  void populateRefData();
  void findContigsContainingKnownAlus();
  void pullContigAlignments();
  void findReadsContainingPolyTails(int32_t);
  void findReadsContainingHeads();

  bool bedFilter(contigAlignment &);
  void writeBedPEHeader(std::ofstream &);
  void writeContigVecToBedPE(std::ofstream &);

  void writeToVCF(std::string &);
  void writeContigVecToVCF(std::ofstream &);
  
  void writeToBed(std::string &);
  
  void printContigVec();

};

#endif // KNOWNALUS_H
