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
#include "denovo.h"


struct contig;
struct contigAlignment;

class KnownAlus{
 public:
  KnownAlus(std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::vector<std::string>, std::string);
  ~KnownAlus();

  const std::vector<contig> getContigVec();
  
 private:
  const std::string rawBamPath_;
  const std::string contigFastqPath_;
  const std::string contigBamPath_;
  const std::string aluFastaPath_;
  const std::string aluIndexPath_;
  const std::string refPath_;
  const std::string refIndexPath_;
  const std::string stub_;
  const std::string vcfOutPath_;
  const std::string fastaHackPath_;

  std::vector<contig> contigVec_;
  std::vector<std::string> parentBams_;
  std::vector<contigAlignment> parentContigAlignments_;
  std::vector<BamTools::RefData> refData_;
  std::string getChromosomeFromRefID(int32_t);

  void populateRefData();
  void findContigsContainingKnownAlus();
  void pullContigAlignments();
  void writeContigVecToBedPE(std::ofstream &);

  void writeToVCF(std::string &);
  void writeContigVecToVCF(std::fstream &);

  std::vector<contigAlignment> findParentContigAlignments(const BamTools::BamAlignment &, const BamTools::BamRegion &, const std::vector<std::string> &);
  std::vector<contigAlignment> populateParentContigAlignments(std::vector<contigAlignment>);

  const bool isDenovo(const std::vector<contigAlignment> &);
  
  void writeToBed(std::string &);
  void printContigVec();

};

#endif // KNOWNALUS_H
