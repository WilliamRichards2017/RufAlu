#include "contig.h"
#include "vcfWriter.h"


void vcfWriter::populateVCFLine(){
  vcfLine_.chrom = ca_.chrom;
  vcfLine_.pos = ca_.alignedContig.Position;
  //TODO: figure out how to report as denovo or inherited
  vcfLine_.id = "Unknown";
  //TODO: //write function to get nucleotide at alu head start pos
  vcfLine_.ref = "Unknown";

  //TODO: write function to get type of mobile element (parse out aluHit list)
  vcfLine_.alt = "INS:ME:ALU/INS:ME:L1";
  vcfLine_.info = ".";
  vcfLine_.contigName = ca_.alignedContig.Name;
  //TOOD: write function to parse cigar data into string"
  vcfLine_.cigarString = "CIGAR";
  vcfLine_.qual = ca_.MapQuality;
}

vcfWriter::vcfWriter(contig & ca, std::ofstream & vcfStream, std::string & stub) : ca_(ca), vcfStream_(vcfStream), stub_(stub) {
  if(!vcfStream.is_open()){
    std::string bs = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/out/" + stub_ + ".vcf";
    vcfStream_.open(bs);
  }
  vcfWriter::writeVCFLine();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeVCFLine(){
  std::cout << vcfLine_.chrom << '\t' << vcfLine_.id << '\t' << vcfLine_.ref << '\t' << vcfLine_.alt 
	    << '\t' << vcfLine_.info << '\t' << vcfLine_.contigName << '\t' << vcfLine_.cigarString 
	    << '\t' << vcfLine_.qual << std::endl;
}

void vcfWriter::writeVCFHeader(){
  vcfStream_ << "##fileformat=VCFv4.1" << std::endl;
  vcfStream_ << "##fileDate=" << std::time(0) << std::endl;
  vcfStream_ << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=LP,Number=1,Type=Integer,Description=\"Number of lowcoverage parent bases\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=PC,Number=1,Type=Integer,Description=\"Mode of parents coverage\">" << std::endl;
  vcfStream_ << "##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">" << std::endl;
  vcfStream_ << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << std::endl;
  vcfStream_ << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << std::endl; 
  vcfStream_ << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << std::endl; 
  vcfStream_ << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">" << std::endl;
  vcfStream_ << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">"<< std::endl;
  vcfStream_ << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">"<< std::endl;
  vcfStream_ << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">"<< std::endl;
  vcfStream_ << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">"<< std::endl;
  vcfStream_ << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">"<< std::endl;
  vcfStream_ << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">"<< std::endl;
  vcfStream_ << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
  vcfStream_ << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;
  vcfStream_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
  vcfStream_ << stub_ << endl;
  
}
