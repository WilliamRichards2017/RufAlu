#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
#include "kseq.h"
#include "minimap.h"
#include <vector>
#include <string>
#include <iostream>
#include <utility> // std::pair

#include "intersect.h" //contigWindow struct
#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct bedPELine {
  int32_t chrom1;
  int32_t chrom1Start;
  int32_t chrom1End;
  int32_t chrom2;
  int32_t chrom2Start;
  int32_t chrom2End;
  std::string name_rufus_contig;
  std::string name_alu_hit;
  int32_t score_numHits;
  int32_t longestTail;
};

class KnownAlus{
 public:
  KnownAlus(std::string, std::string, std::string, const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
  std::vector<fastqRead> * getContigsContainingKnownAlus();
  

 private:
  std::string contigFilePath_;
  std::string contigBamPath_;
  std::string  mutationPath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  const char * refPath_;
  const char * refIndexPath_;

  std::string stub_;

  

  std::vector<fastqRead> *contigsContainingKnownAlus_;
  const std::vector<BamTools::RefData> *refData_;

  void populateRefData(std::string);
  void writeHitToBed(std::ofstream&, bedPELine *);
  void findContigsContainingKnownAlus();
  void alignContigsContainingKnownAlus(const char *);
  void mapContigsToRef(const char *);
  void findContigsContainingPolyATails(const char *);
  void findReadsContainingPolyATails(std::vector<contigWindow>, std::string);
  std::string getChromosomeFromRefID(int32_t);

};

#endif // KNOWNALUS_H
