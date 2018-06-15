
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
  vcfLine_.chrom = ca_.chrom;
  vcfLine_.pos = ca_.clipPeak;
  //TODO: figure out how to report as denovo or inherited
  vcfLine_.id = "denvo/inherited";
  //TODO: //write function to get nucleotide at alu head start pos
  vcfLine_.ref = "N";

  vcfLine_.alt = "INS:ME:"+ca_.aluHit;
  vcfLine_.info = ".";
  vcfLine_.contigName = ca_.alignedContig.Name;
  //TOOD: write function to parse cigar data into string"
  for(auto it = std::begin(ca_.alignedContig.CigarData); it != std::end(ca_.alignedContig.CigarData); ++it){
    vcfLine_.cigarString += it->Type;
    vcfLine_.cigarString += std::to_string(it->Length);
  }
  vcfLine_.qual = ca_.alignedContig.MapQuality;
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

void vcfWriter::writeVCFLine(){
  vcfStream_ << vcfLine_.chrom << '\t' << vcfLine_.pos << '\t'  << vcfLine_.id << '\t' << vcfLine_.ref << '\t' << vcfLine_.alt 
	     << '\t' << vcfLine_.qual << '\t' <<vcfLine_.info << '\t' << vcfLine_.contigName << '\t' << vcfLine_.cigarString 
	     << '\t' << std::endl;
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
  vcfStream << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< std::endl;
  vcfStream << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">"<< std::endl;
  vcfStream << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">"<< std::endl;
  vcfStream << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">"<< std::endl;
  vcfStream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">"<< std::endl;
  vcfStream << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">"<< std::endl;
  vcfStream << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
  vcfStream << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;
  vcfStream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  vcfStream << stub << std::endl;
  
}
