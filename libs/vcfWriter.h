#ifndef __VCF_WRITER_H__
#define __VCF_WRITER_H__

#include "contig.h"
#include "denovo.h"

#include <iostream>
#include <time.h>
#include <string>

  
struct filterField{
  bool DS; //FILTER=<ID=DS,Description="polyA tails detected on both forward and reverse strand">
  bool SB; //FILTER=<ID=SB,Description="polyA tails on only 1 side of the events"> 
  int32_t polyA; //FILTER=<ID=SB,Description="atleast 3 polyA tails that start at the same point"> 
};
struct infoField {
  int32_t NR  =  -1; // INFO=<ID=NR,Number=1,Type=Integer,Description="Number of total reads in target region">
  int32_t NT  =  -1; // INFO=<ID=NT,Number=1,Type=Integer,Description="Number of polyA tails in target region">
  std::string TB;  // INFO=<ID=TB,Number=1,Type=Integer,Description="Is tail left bound, right bound, or double bound"> 
  int32_t NH = -1; // INFO=<ID=NH,Number=1,Type=Integer,Description="Number of alu heads in target region"> 
  int32_t LT  =  -1; // INFO=<ID=LT,Number=1,Type=Integer,Description="Longest polyA tail in target region"> 
  std::string SVTYPE; // INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV detected">
  int32_t SVLEN = -1; // INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV detected"> 
  int32_t END = -1; // INFO=<ID=END,Number=1,Type=Integer,Description="END of SV detected">
  std::string RN; // INFO=<ID=RN,Number=1,Type=String,Description="Name of contig that produced the call">
  int16_t MQ = -1; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::string cigar; 
  std::string CVT = "ME"; //Compressed variant type
  std::vector<int32_t> HD; // hashcount for kmers overlapping variant
  std::pair<bool, bool> GT; //Genotype information
  int32_t DP = -1; // Total kmer depth
  int32_t RO = -1; // reference kmer count
  int32_t AO = -1; // Altername kmer count
  int32_t LP = -1; //number of low coverage bases
  int32_t PC = -1; // Mode of parent coverage
  float SB = -1.0; // Strand Bias

};

struct formatField {
  int32_t DP = -1; // FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">
  int32_t RO = -1; // FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\"
  int32_t AO = -1; // FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">
  int32_t LP = -1; // FORMAT=<ID=LP,Number=1,Type=Integer,Description=\"Number of lowcoverage parent bases\"
  int32_t PC = -1; // FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Mode of parents coverage\">
  float SB = -1.0; // FORMAT=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">
};

struct vcfLine {
  std::string CHROM = ".";
  int32_t POS = -1;
  std::string ID = ".";
  std::string REF = ".";
  std::string ALT = ".";
  int32_t QUAL = -1;
  filterField FILTER = {};
  infoField INFO = {};
  formatField FORMAT = {};
};




class vcfWriter{
 public:
  vcfWriter(const contigAlignment &, std::fstream &, const std::string &);
  ~vcfWriter();

  const bool vcfFilter();
  void writeVCFLine();

  static void writeVCFHeader(std::fstream &, const std::string &);

  
 private:
  
  vcfLine vcfLine_ = {};
  std::fstream & vcfStream_;
  
  contigAlignment ca_;
  const std::string stub_;
  
  void populateVCFLine();  

  void writeFilter();
  void writeInfo();

};

#endif //__VCF_WRITER_H__
