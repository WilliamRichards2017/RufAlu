#ifndef KNOWNALUS_H
#define KNOWNALUS_H

#include "fastqParse.h"
#include "kseq.h"
#include "minimap.h"
#include <vector>
#include <string>
#include <utility> // std::pair
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"



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

class KnownAlus{
 public:
  KnownAlus(const char *, const char *, const char *, const char *, const char *);
  ~KnownAlus();
  
  std::vector<fastqRead *> * getContigsContainingKnownAlus();
  

 private:
  const char * contigFilePath_;
  const char * aluFilePath_;
  const char * aluIndexPath_;
  const char * refPath_;
  const char * refIndexPath_;

  std::vector<fastqRead *> *contigsContainingKnownAlus_;
  void findContigsContainingKnownAlus();
  void alignContigsContainingKnownAlus(const char *);
  void mapContigsToRef(const char *);
  void recoverPolyATails();
  void findContigsContainingPolyATails(const char *);
  void findReadsContainingPolyATails(std::vector<BamTools::BamAlignment>);

};

#endif // KNOWNALUS_H
