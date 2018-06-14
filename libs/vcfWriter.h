#ifndef __VCF_WRITER_H__
#define __VCF_WRITER_H__

#include "contig.h"

#include <iostream>
#include <time.h>
#include <string>

struct vcfLine{
  
  std::string chrom= ".";
  int32_t pos = -1;
  std::string id = ".";
  std::string ref = ".";
  std::string alt = ".";
  std::string info = ".";
  std::string contigName= ".";
  std::string cigarString = ".";
  int16_t qual= -1;

};

class vcfWriter{
 public:
  vcfWriter(contigAlignment &, std::ofstream &, std::string);
  ~vcfWriter();

  void writeVCFHeader();
  
 private:
  vcfLine vcfLine_ = {};
  void populateVCFLine();
  const contigAlignment & ca_;
  std::ofstream & vcfStream_;
  std::string stub_;
  
};

#endif //__VCF_WRITER_H__
