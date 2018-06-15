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
  std::string cigarString = "";
  int16_t qual= -1;

};

class vcfWriter{
 public:
  vcfWriter(const contigAlignment &, std::ofstream &, const std::string &);
  ~vcfWriter();

  const bool vcfFilter();
  void writeVCFLine();

  static void writeVCFHeader(std::ofstream &, const std::string &);


  
 private:
  vcfLine vcfLine_ = {};
  std::ofstream & vcfStream_;
  
  const contigAlignment & ca_;
  const std::string stub_;
  
  void populateVCFLine();  

};

#endif //__VCF_WRITER_H__
