#include <stdexcept>
#include <string>
#include <time.h>

#include "contig.h"
#include "vcfWriter.h"


//TODO: Implement vcf filter
const bool vcfWriter::vcfFilter(){
  return true;
}

void vcfWriter::populateVCFLine(){
  vcfLine_.CHROM = ca_.chrom;
  vcfLine_.POS = ca_.clipCoords_.clipStart+ca_.alignedContig.Position;
  //TODO: figure out how to report as denovo or inherited
  vcfLine_.ID = "denvo/inherited";
  //TODO: //write function to get nucleotide at alu head start pos
  vcfLine_.REF = "N";
  vcfLine_.ALT = "INS:ME:"+ca_.aluHit;  
  vcfLine_.QUAL = ca_.alignedContig.MapQuality;

  vcfLine_.FILTER.SB = (ca_.leftBound) ^ (ca_.rightBound); // ^ = XOR
  vcfLine_.FILTER.DS = (ca_.leftBoundDS ^ ca_.rightBoundDS) and vcfLine_.FILTER.SB;
  
  vcfLine_.INFO.SVTYPE = "INS";
  vcfLine_.INFO.SVLEN = std::abs(ca_.clipCoords_.clipStart - ca_.clipCoords_.clipEnd);
  vcfLine_.INFO.END = ca_.clipCoords_.clipEnd + ca_.alignedContig.Position;
  vcfLine_.INFO.HD = ca_.maxHash;
  vcfLine_.INFO.MQ = ca_.alignedContig.MapQuality;
  vcfLine_.INFO.RN = ca_.alignedContig.Name;
  if(ca_.leftBound){
    vcfLine_.INFO.NT = ca_.leftBoundTails.size();
  }
  else if (ca_.rightBound){
    vcfLine_.INFO.NT = ca_.leftBoundTails.size();
  }
  vcfLine_.INFO.NR = ca_.readsInRegion;
  vcfLine_.INFO.LT = util::getLongestTail(ca_.leftBoundTails, ca_.rightBoundTails);


  for(auto it = std::begin(ca_.alignedContig.CigarData); it != std::end(ca_.alignedContig.CigarData); ++it){
    vcfLine_.INFO.cigar += it->Type;
    vcfLine_.INFO.cigar += std::to_string(it->Length);
  }

}

vcfWriter::vcfWriter( const contigAlignment & ca, std::ofstream & vcfStream, const std::string & stub) : ca_(ca), vcfStream_(vcfStream), stub_(stub) {
  if(!vcfStream_.is_open()){
    std::cerr << "vcfStream is not open, exiting run with non-zero exit status " << std::endl;
    exit (EXIT_FAILURE);
  }
  vcfWriter::populateVCFLine();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeFilter(){
  if (vcfLine_.FILTER.DS and vcfLine_.FILTER.SB){
    vcfStream_ << "PASS\t";
  }
  else if((!vcfLine_.FILTER.DS) and (!vcfLine_.FILTER.SB)) {
    vcfStream_ << "DS:SB\t";
  }
  else if (!vcfLine_.FILTER.SB){
    vcfStream_ << "SB\t";
  }
  else if (!vcfLine_.FILTER.DS){
    vcfStream_ << "DS\t";
  }
  else{
    vcfStream_ << '\t';
  }
}

void vcfWriter::writeInfo(){
  vcfStream_ << "NR=" << vcfLine_.INFO.NR << ";NT=" << vcfLine_.INFO.NT << ";LT=" << vcfLine_.INFO.LT <<  ";SVTYPE=" << vcfLine_.INFO.SVTYPE 
	     << ";SVLEN=" << vcfLine_.INFO.SVLEN << ";END=" << vcfLine_.INFO.END << ";HD=" << vcfLine_.INFO.HD << ";RN=" << vcfLine_.INFO.RN << ";cigar=" << vcfLine_.INFO.cigar << ';';
}

void vcfWriter::writeVCFLine(){
  vcfStream_ << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t'  << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT 
	     << '\t' << vcfLine_.QUAL << '\t';
  
  vcfWriter::writeFilter();
  vcfWriter::writeInfo();

  vcfStream_ << std::endl;


}

void vcfWriter::writeVCFHeader(std::ofstream & vcfStream, const std::string & stub){
  vcfStream << "##fileformat=VCFv4.1" << std::endl;
  vcfStream << "##fileDate=" << time(0) << std::endl;
  vcfStream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  vcfStream << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << std::endl;
  vcfStream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << std::endl;
  vcfStream << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << std::endl;
  vcfStream << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << std::endl;
  vcfStream << "##FORMAT=<ID=LP,Number=1,Type=Integer,Description=\"Number of lowcoverage parent bases\">" << std::endl;
  vcfStream << "##FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Mode of parents coverage\">" << std::endl;
  vcfStream << "##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">" << std::endl;
  vcfStream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << std::endl;
  vcfStream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << std::endl; 
  vcfStream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << std::endl; 
  vcfStream << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">" << std::endl;
  vcfStream << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">" << std::endl;
  vcfStream << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">" << std::endl;
  vcfStream << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">" << std::endl;
  vcfStream << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">" << std::endl;
  vcfStream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">" << std::endl;
  vcfStream << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">" << std::endl;
  vcfStream << "##INFO=<ID=NT,Number=1,Type=Integer,Description=\"Number of polyA tails in target region\">" << std::endl; 
  vcfStream << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of total reads in target region\">" << std::endl;
  vcfStream << "##INFO=<ID=LT,Number=1,Type=Integer,Description=\"Longest polyA tail in target region\">" << std::endl;
  vcfStream << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
  vcfStream << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;
  vcfStream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << stub << std::endl;
}

