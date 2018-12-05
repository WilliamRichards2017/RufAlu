#include <stdexcept>
#include <string>
#include <time.h>

#include "contig.h"
#include "vcfWriter.h"


//TODO: Implement vcf filter
const bool vcfWriter::vcfFilter(){
  return true;
}

void vcfWriter::printContigAlignment(){
  
}

void vcfWriter::populateGenotypes(){
  vcfLine_.INFO.probandGT = ca_.getProbandGT();
  for(auto & pDE : ca_.getDenovoVec()){
    genotypeField gt = {};
    gt.genotype = pDE.getGenotype();
    gt.DP = pDE.DP_;
    gt.RO = pDE.RO_;
    gt.AO = pDE.AO_;
    vcfLine_.INFO.parentGTs.push_back(gt);
  }
}

void vcfWriter::populateVCFLine(){
  vcfLine_.CHROM = ca_.getChrom();
  //vcfLine_.POS = ca_.getClipCoords().clipStart+ca_.getAlignedContig().Position;
  if(ca_.isReadLeftBound()){
    vcfLine_.POS = ca_.getAlignedContig().Position + 1;
  }
  else{
    vcfLine_.POS = ca_.getAlignedContig().GetEndPosition();
  }
  
  if(ca_.isDenovo()){
    vcfLine_.ID = "ME-DeNovo";
  }
  else{
    vcfLine_.ID = "Inherited";
  }

  //TODO: //write function to get nucleotide at alu head start pos
  vcfLine_.REF = "N";
  vcfLine_.ALT = "INS:ME:"+ca_.getAluHit().first;  
  vcfLine_.QUAL = ca_.getAlignedContig().MapQuality;

  //vcfLine_.FILTER.SB = (ca_.tailLeftBound) ^ (ca_.tailRightBound); // ^ = XOR
  vcfLine_.FILTER.TDS = ca_.isTailDoubleStranded();
  vcfLine_.FILTER.HDS = ca_.isHeadDoubleStranded();

  if(!vcfLine_.FILTER.TDS or !vcfLine_.FILTER.HDS){
    vcfLine_.ID = "ME-StrandBias";
  }
  
  vcfLine_.INFO.SVTYPE = "INS";
  vcfLine_.INFO.SVLEN = std::abs(ca_.getClipCoords().clipStart - ca_.getClipCoords().clipEnd);
  vcfLine_.INFO.END = ca_.getClipCoords().clipEnd + ca_.getAlignedContig().Position;
  vcfLine_.INFO.MQ = ca_.getAlignedContig().MapQuality;
  vcfLine_.INFO.RN = ca_.getAlignedContig().Name;
  vcfLine_.INFO.NT = ca_.getConsensusTails().size();
  
  vcfLine_.INFO.SB = static_cast<double>(ca_.getForwardStrandCount())/static_cast<double>(ca_.getReadsInRegion());

  vcfLine_.INFO.NH = ca_.getHeads().size(); 
  vcfLine_.INFO.NR = ca_.getReadsInRegion();
  vcfLine_.INFO.LT = util::getLongestTail(ca_.getConsensusTails());

  vcfLine_.INFO.cigar = ca_.getCigarString();
  
  
  vcfWriter::populateGenotypes();

}

vcfWriter::vcfWriter(contigAlignment & ca, std::fstream & vcfStream, const std::string & probandBam) : ca_(ca), vcfStream_(vcfStream), probandBam_(probandBam) { 
  if(!vcfStream_.is_open()){
    std::cerr << "vcfStream is not open, exiting run with non-zero exit status " << std::endl;
    exit (EXIT_FAILURE);
  }

  vcfWriter::populateVCFLine();
}

vcfWriter::~vcfWriter(){
}

void vcfWriter::writeFilter(){
  
  std::cout << "\n\n WRITING FILTER" << std::endl;
  std::cout << "tds is : " << vcfLine_.FILTER.TDS << std::endl;
  std::cout << "hds is : " << vcfLine_.FILTER.HDS << std::endl << std::endl;
  if (vcfLine_.FILTER.TDS && vcfLine_.FILTER.HDS){
    vcfStream_ << "PASS\t";
  }
  else if(!vcfLine_.FILTER.TDS && vcfLine_.FILTER.HDS) {
    vcfStream_ << "TDS\t";
  }
  else if (!vcfLine_.FILTER.HDS && vcfLine_.FILTER.TDS) {
      vcfStream_ << "HDS\t";
  }
  else{
    vcfStream_ << "TDS,HDS\t";
  }
}

void vcfWriter::writeGenotypes(){
  vcfStream_ << "\tGT:DP:RO:AO " << vcfLine_.INFO.probandGT.genotype.first << '/' << vcfLine_.INFO.probandGT.genotype.second << ':' << vcfLine_.INFO.probandGT.DP << ':' << vcfLine_.INFO.probandGT.RO << ':' << vcfLine_.INFO.probandGT.AO << '\t';

  for(auto pGT : vcfLine_.INFO.parentGTs){
    vcfStream_ << pGT.genotype.first << '/' << pGT.genotype.second << ':' << pGT.DP << ':' << pGT.RO << ':' << pGT.AO << '\t';
  }
  vcfStream_ << std::endl;

}

void vcfWriter::writeInfo(){
  vcfStream_ << "NR=" << vcfLine_.INFO.NR << ";NT=" << vcfLine_.INFO.NT <<  ";NH=" << vcfLine_.INFO.NH << ";LT=" << vcfLine_.INFO.LT <<  ";SVTYPE=" << vcfLine_.INFO.SVTYPE << ";SVLEN=" << vcfLine_.INFO.SVLEN << ";END=" << vcfLine_.INFO.END << ";RN=" << vcfLine_.INFO.RN << ";cigar=" << vcfLine_.INFO.cigar << ";SB=" << vcfLine_.INFO.SB <<  ";CVT=" << vcfLine_.INFO.CVT << ";HD=";
  for(auto h : vcfLine_.INFO.HD){
    vcfStream_ << h << '_';
  }
  vcfStream_ << ";AO=" << vcfLine_.INFO.probandGT.AO << ";VT=" << vcfLine_.INFO.CVT << ';';

  vcfWriter::writeGenotypes();
  
}

void vcfWriter::writeVCFLine(){

  if(vcfWriter::vcfFilter()){
  vcfStream_ << vcfLine_.CHROM << '\t' << vcfLine_.POS << '\t'  << vcfLine_.ID << '\t' << vcfLine_.REF << '\t' << vcfLine_.ALT 
	     << '\t' << vcfLine_.QUAL << '\t';
  
  vcfWriter::writeFilter();
  vcfWriter::writeInfo();
  
  vcfStream_ << std::endl;
  }
}

void vcfWriter::writeVCFHeader(std::fstream & vcfStream, std::string probandBam){
  vcfStream << "##fileformat=VCFv4.1" << std::endl;
  vcfStream << "##fileDate=" << time(0) << std::endl;
  vcfStream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  vcfStream << "##FORMAT=<ID=AK,Number=1,Type=Integer,Description=\"Alternte Kmer Count\">" << std::endl;
  vcfStream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Kmer depth across the variant\">" << std::endl;
  vcfStream << "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Mode of reference kmer counts\">" << std::endl;
  vcfStream << "##FORMAT=<ID=AO,Number=1,Type=Integer,Description=\"Mode of alt kmer counts\">" << std::endl;
  vcfStream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV detected\">" << std::endl;
  vcfStream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV detected\">" << std::endl; 
  vcfStream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"END of SV detected\">" << std::endl; 
  vcfStream << "##INFO=<ID=AO,Number=1,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">" << std::endl;
  vcfStream << "##INFO=<ID=HD,Number=.,Type=String,Description=\"Hash counts for each k-mer overlapping the vareint, -1 indicates no info\">" << std::endl;
  vcfStream << "##INFO=<ID=RN,Number=1,Type=String,Description=\"Name of contig that produced the call\">" << std::endl;
  vcfStream << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the contig that created the call\">" << std::endl;
  vcfStream << "##INFO=<ID=cigar,Number=1,Type=String,Description=\"Cigar string for the contig that created the call\">" << std::endl;
  vcfStream << "##INFO=<ID=SB,Number=1,Type=Float,Description=\"StrandBias\">" << std::endl;
  vcfStream << "##INFO=<ID=VT,Number=1,Type=String,Description=\"Varient Type\">" << std::endl;
  vcfStream << "##INFO=<ID=CVT,Number=1,Type=String,Description=\"Compressed Varient Type\">" << std::endl;
  vcfStream << "##INFO=<ID=NT,Number=1,Type=Integer,Description=\"Number of polyA tails in target region\">" << std::endl; 
  vcfStream << "##INFO=<ID=TB,Number=1,Type=Integer,Description=\"Is tail left bound, right bound, or double bound\">" << std::endl;
  vcfStream << "##INFO=<ID=NH,Number=1,Type=Integer,Description=\"Number of alu heads in target region\">" << std::endl;
  vcfStream << "##INFO=<ID=NR,Number=1,Type=Integer,Description=\"Number of total reads in target region\">" << std::endl;
  vcfStream << "##INFO=<ID=LT,Number=1,Type=Integer,Description=\"Longest polyA tail in target region\">" << std::endl;
  vcfStream << "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">" << std::endl;
  vcfStream << "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">" << std::endl;
  vcfStream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << probandBam << std::endl;
}

