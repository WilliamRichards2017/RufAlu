#ifndef __VCF_WRITER_H__
#define __VCF_WRITER_H__

#include "contig.h"
#include "denovo.h"

#include <iostream>
#include <time.h>
#include <string>

struct genotypeField{
  std::pair<bool, bool> GT = std::make_pair(-1,-1); //Genotype information
  int32_t DP = -1; // Total kmer depth
  int32_t RO = -1; // reference kmer count
  int32_t AO = -1; // Altername kmer count
};
  
struct filterField{
  bool HDS;
  bool TDS; //FILTER=<ID=DS,Description="polyA tails detected on both forward and reverse strand">
  int32_t polyA; //FILTER=<ID=SB,Description="atleast 3 polyA tails that start at the same point"> 
};
struct infoField {
  int32_t NR  =  -1; // INFO=<ID=NR,Number=1,Type=Integer,Description="Number of total reads in target region">
  int32_t NT  =  -1; // INFO=<ID=NT,Number=1,Type=Integer,Description="Number of polyA tails in target region">
  int32_t NH = -1; // INFO=<ID=NH,Number=1,Type=Integer,Description="Number of alu heads in target region"> 
  int32_t LT  =  -1; // INFO=<ID=LT,Number=1,Type=Integer,Description="Longest polyA tail in target region"> 
  std::string SVTYPE; // INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV detected">
  int32_t SVLEN = -1; // INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV detected"> 
  int32_t END = -1; // INFO=<ID=END,Number=1,Type=Integer,Description="END of SV detected">
  std::string RN; // INFO=<ID=RN,Number=1,Type=String,Description="Name of contig that produced the call">
  int16_t MQ = -1; // INFO=<ID=MQ,Number=1,Type=Integer,Description="Mapping quality of the contig that created the call">
  std::string cigar = ""; 
  std::string CVT = "ME"; //Compressed variant type
  double SB; // Strand Bias

  //TODO - REMOVE HARD CODING WHEN HD IS POPULATED
  std::vector<int32_t> HD = {-1,-1}; // hashcount for kmers overlapping variant
  //END TODO

  genotypeField probandGT;
  std::vector<genotypeField> parentGTs;

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
  vcfWriter(contigAlignment &, std::fstream &, const std::string &);
  ~vcfWriter();

  const bool vcfFilter();
  void writeVCFLine();

  static void writeVCFHeader(std::fstream &, std::string probandBam);

  
 private:
  
  vcfLine vcfLine_ = {};
  std::fstream & vcfStream_;
  contigAlignment ca_;
  const std::string & probandBam_;


  void printContigAlignment();
  void populateVCFLine();  
  void populateGenotypes();
  void populateParentGenotypes();

  void writeGenotypes();
  void writeFilter();
  void writeInfo();

};

#endif //__VCF_WRITER_H__
