#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <utility> //std::pair

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct coords {
  int leftPos;
  int RightPos;
  const char * chrom;
  
};

struct fastqRead {
  std::string name;
  std::string seq;
  std::string qual;
fastqRead(std::string n, std::string s, std::string q) : name(n), seq(s), qual(q) {}
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

class util{

 public:
  static void exec(char const*);
  static bool overlap(std::pair<int, int>, std::vector<std::pair<int, int> >);
  static std::vector<BamTools::BamAlignment> intersect(const char *, const char *);
  static std::vector<BamTools::BamAlignment> intersectBams(const char *, const char *);
  static const char * contigsToFastq(std::vector<fastqRead>*, const char *);
  static std::string baseName(std::string);

  static const std::vector<std::string> Split(const std::string&, const char);
  static const char * getRootDirectory(std::string);

 private:

};

#endif // UTIL_H
