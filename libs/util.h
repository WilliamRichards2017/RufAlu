#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <utility> //std::pair

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "knownAlus.h"
#include "intersect.h"

struct coords {
  int leftPos;
  int RightPos;
  const char * chrom;
  
};

struct contigWindow{
  std::vector<BamTools::BamAlignment> window;
};

struct samRead {
  std::string qname;
  short flag;
  int rname;
  long pos;
  short mapq;
  std::string cigar;
  std::string mrnm;
  long mpos;
  int isize;
  std::string seq;
  std::string qual;
  std::string tags;
};

struct contig {
  std::string name;
  std::string seq;
  std::string qual;
  std::vector<std::string> alusHit;
  std::vector<BamTools::BamAlignment> contigAlignments;
  std::vector<BamTools::BamRegion> contigAlignmentRegions;
  std::vector<contigWindow> overlapingReads;
  std::vector<BamTools::BamAlignment> supportingReads;
  uint32_t longestTail;
};

struct contigWindow;

class util{

 public:
  static void exec(char const*);
  static bool overlap(std::pair<int, int>, std::vector<std::pair<int, int> >);
  static std::vector<BamTools::BamAlignment> intersect(const char *, const char *);
  static std::vector<BamTools::BamAlignment> intersectBams(const char *, const char *);
  static std::string contigsToFastq(std::vector<contig> , std::string);
  static std::string baseName(std::string);
  static void printContigWindow(contigWindow c);


  static const std::vector<std::string> Split(const std::string&, const char);
  static const char * getRootDirectory(std::string);

 private:

};

#endif // UTIL_H
