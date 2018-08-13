#include <algorithm>
#include <assert.h>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <zlib.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include "aluHead.h"
#include "contig.h"
#include "knownAlus.h"
#include "denovo.h"
#include "kseq.h"
#include "minimap.h"
#include "polyATail.h"
#include "util.h"
#include "vcfWriter.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

KSEQ_INIT(gzFile, gzread)

bool debugPrintFilter(contigAlignment & ca){
  return true;
  if(ca.getAlignedRegion().LeftPosition != -1 and ca.getAlignedRegion().RightPosition != -1){
    return true;
  }
  return false;
}

void printContigAlignment(contigAlignment & ca){
  //  if(ca.getAlignedContig().Name.compare("NODE_8800.bam.generator.V2_372_L151_D6:6:0::MH0") == 0){
    std::cout << std::endl;
    std::cout << "printing contig alignment: " << ca.getAlignedContig().Name;
    std::cout << "polyAails.size(): " << ca.getConsensusTails().size();
    std::cout << "aluHeads.size(): " << ca.getHeads().size();
    std::cout << std::endl;
    //}
}

void printContig(contig & c){
  for(auto it = std::begin(c.contigAlignments); it != std::end(c.contigAlignments); ++it){
    printContigAlignment(*it);
  }
  
}


void KnownAlus::printContigVec(){
  for(auto it = std::begin(contigVec_); it != std::end(contigVec_); ++it){
    printContig(*it);
  }
}

void KnownAlus::populateRefData(){
  BamTools::BamReader reader;
  if (!reader.Open(rawBamPath_)){
    std::cerr << "Could not open raw reads Bam file" << rawBamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
  reader.Close();
}


std::string KnownAlus::getChromosomeFromRefID(int32_t id){
  std::string ret = "";
  
  if(id == -1) {
    ret = "unmapped";
    return ret;
  }
  ret = refData_[id].RefName;
  return ret;
}

void KnownAlus::writeToVCF(std::string & vcfFile){

  std::fstream vcfStream;
  vcfStream.open(vcfFile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
  //vcfWriter::writeVCFHeader(vcfStream, stub_);
  KnownAlus::writeContigVecToVCF(vcfStream);
  vcfStream.close();
}

void KnownAlus::writeContigVecToVCF(std::fstream & vcf){
  //for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
  for(auto  c : contigVec_){
    for(auto  ca : c.contigAlignments){
      vcfWriter writer = {ca, vcf, stub_};
	if(ca.getDenovoVec().size() > 0){
	  std::cout << "CA is denovo in KnownAlus " << std::endl;
	}
	if(writer.vcfFilter()){
	  writer.writeVCFLine();
      }
    }
  }
}
  

void KnownAlus::findContigsContainingKnownAlus()
{
   mm_idxopt_t iopt;
   mm_mapopt_t mopt;
   int n_threads = 3;

   //std::cout << "reading in contig fasta file " << contigFilePath_ << std::endl;


   mm_verbose = 3; // print to std out
   mm_set_opt(0, &iopt, &mopt); //initialize alignment parameters to default
   mopt.flag |= MM_F_CIGAR; // perform alignment                                                                                                                                                                                               

   gzFile f = gzopen(contigFastqPath_.c_str(), "r");
   assert(f);
   kseq_t *ks = kseq_init(f);

   // open index reader
   mm_idx_reader_t *r = mm_idx_reader_open(aluFastaPath_.c_str(), &iopt, 0);
   mm_idx_t *mi;
   while ((mi = mm_idx_reader_read(r, n_threads)) != 0) { // traverse each part of the index
     mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence
     mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread                                                                                                                              

     while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence                                                                                                                                                           
       //std::cout << "looping through kseq_reads" << std::endl; 

       mm_reg1_t *reg;
       int j, i, n_reg;

       reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0); // get all hits for the query

       contig c = {};
       c.name = std::string(ks->name.s);
       c.seq = ks->seq.s;
       
       //std::cout << "Checking contig for alu Hits: " << ks->name.s << std::endl;
       
       for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	 mm_reg1_t *r = &reg[j];
	 c.alusHit.push_back(std::make_pair(mi->seq[r->rid].name, int(r->mapq)));
	 free(r->p);
	 //std::cout << "found alu hit for contig: " << ks->name.s << std::endl;
       }
       
       if(c.alusHit.size() > 0){
	 // std::cout << "found contig containing knonw alu: " <<  c.name << std::endl;
	 contigVec_.push_back(c);
       }
       free(reg);
     }
     mm_tbuf_destroy(tbuf);
     mm_idx_destroy(mi);
   }
   mm_idx_reader_close(r); 
   kseq_destroy(ks); // close the query file                                                                                                                                                                                                   
   gzclose(f);
 }


KnownAlus::KnownAlus(std::string rawBamPath, std::string contigFastqPath, std::string contigBamPath, std::string aluFastaPath, std::string aluIndexPath, std::string refPath, std::string refIndexPath, std::string vcfOutPath, std::vector<std::string> parentBams) :  rawBamPath_(rawBamPath), contigFastqPath_(contigFastqPath), contigBamPath_(contigBamPath), aluFastaPath_(aluFastaPath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath), stub_(util::baseName(rawBamPath)), vcfOutPath_(vcfOutPath), parentBams_(parentBams){
   
  contigVec_ = {};
  refData_ = {};

  std::cerr << "[1/7]  Populating reference data for " << stub_ << std::endl;
  KnownAlus::populateRefData();

  std::cerr << "[2/7]  Finding contigs containing known alus for " << stub_ << std::endl;
  KnownAlus::findContigsContainingKnownAlus();

  std::cerr << "[3/7]  Pulling contig hit alignments for " << stub_ <<  std::endl;
  KnownAlus::pullContigAlignments();

  //  std::cerr << "[4/7]  Finding reads containing polyA tails for " << stub_ << std::endl;
  //KnownAlus::findReadsContainingPolyTails(9);

  //std::cerr << "[5/7]  Finding reads containing alu heads tails for " << stub_ << std::endl;
  // KnownAlus::findReadsContainingHeads();

  //std::cerr << "[6/7] Flaging denovos for " << stub_ << std::endl; 
  //KnownAlus::findDenovoEvidence();

 

  //std::cout << "[5/5] Writing out results to bed file " << stub_  << ".bed" << std::endl;
  //KnownAlus::writeToBed(prefix_ + stub_ + "bed");

  std::string tempy = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/testy.vcf";

  KnownAlus::printContigVec();


  std::cerr << "[7/7] Writing out results to vcf file " << vcfOutPath << std::endl;
  //KnownAlus::writeToVCF(vcfOutPath);
  KnownAlus::writeToVCF(tempy);

 }



KnownAlus::~KnownAlus(){
}

void KnownAlus::pullContigAlignments(){

  BamTools::BamReader reader;
  if (!reader.Open(contigBamPath_)){
    std::cout << "Could not open input Bam file " << contigBamPath_ << std::endl;
    exit (EXIT_FAILURE);
  }

  BamTools::BamAlignment al;
  
  for(auto  cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    while(reader.GetNextAlignment(al)){
      if(cvIt->name.compare(al.Name)==0 and al.HasTag("SA")){
	//std::cout << "Checking for peak and clip coord intersection for read: " << al.Name << std::endl;
	//util::printCigar(al.CigarData);
	//std::cout << al.Qualities << std::endl;
	clipCoords clipPeak = util::intersectPeaksAndClips(util::getPeaks(al), util::getLocalClipCoords(al));
	if(clipPeak.clipStart != -1){
	  //std::cout << "Found intersection between peak and clips for " << cvIt->name << std::endl;
	  std::pair<std::string, int32_t> hqAlu = util::getHighestQualityAluHit(cvIt->alusHit);
	  std::string chrom = getChromosomeFromRefID(al.RefID);
	  BamTools::BamRegion region = {al.RefID, al.Position, al.RefID, al.GetEndPosition()};
	  contigAlignment ca = {rawBamPath_, parentBams_, hqAlu, al, chrom, region };
	  cvIt->contigAlignments.push_back(ca);
	}
      }
    }
    reader.Rewind();
  }
  reader.Close();
}
