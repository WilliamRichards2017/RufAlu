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

#include "contig.h"
#include "fastqParse.h"
#include "knownAlus.h"
#include "kseq.h"
#include "minimap.h"
#include "polyATail.h"
#include "intersect.h"
#include "util.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

KSEQ_INIT(gzFile, gzread)

bool debugPrintFilter(contigAlignment & ca){
  if(ca.leftBound || ca.rightBound and ca.readsInRegion < 201){
    return true;
  }
  return false;
}

void printContigAlignment(contigAlignment & ca){
  if(debugPrintFilter(ca)){
    std::cout << "Aligned Contig name: " << ca.alignedContig.Name << std::endl;
    std::cout << "Quality string: " << ca.alignedContig.Qualities << std::endl;
  }
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
    std::cout << "Could not open raw reads Bam file" << rawBamPath_ << std::endl;
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

void KnownAlus::writeBedPEHeader(std::ofstream &bed){
  bed << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' << "chrom" << '\t'<< "contig_name" << '\t' << "alu_hit" << '\t' << "primary_alignment" << '\t' << "numReadsInRegion" << '\t' << "numPolyATails" <<  '\t' << "both_strands" << '\t' << "longest_tail"  << std::endl;
}

bool KnownAlus::bedFilter(contigAlignment & ca) {
  if(ca.readsInRegion < 200 and ca.doubleStranded and (!(ca.leftBound and ca.rightBound))){
    auto peaks = util::getPeaks(ca.alignedContig);
    
      if(ca.leftBound){
	if (ca.rightBoundTails.size() > 1){
	  return false;
	}
	else {
	  return true;
	}
      }
      else{
	if(ca.leftBoundTails.size() > 1){
	  return false;
	}
	else{
	  return true;
	}
      }
    }
  return false;
}

void KnownAlus::writeContigVecToBedPE(std::ofstream &bed){
  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    for(auto caIt = std::begin(cvIt->contigAlignments); caIt != std::end(cvIt->contigAlignments); ++caIt){
      if(KnownAlus::bedFilter(*caIt)) {
	bed << getChromosomeFromRefID(caIt->alignedContig.RefID) << '\t'<< caIt->alignedContig.Position << '\t' << caIt->alignedContig.GetEndPosition() 
	    << '\t'<< cvIt->name << '\t' << cvIt->alusHit[0] << '\t' << caIt->readsInRegion << '\t' << caIt->alignedContig.IsPrimaryAlignment();
	
	if(caIt->leftBound){
	  bed  << '\t' << caIt->leftBoundTails.size() << '\t' << caIt->doubleStranded << std::endl;
	}
	else if(caIt->rightBound){
	  bed << '\t' << caIt->rightBoundTails.size() << '\t' << caIt->doubleStranded << std::endl;
	}
	else {
	  bed << '\t' << "tail not left or right bound..." << std::endl;
	}
      }
    }
  }
}


void KnownAlus::findReadsContainingPolyTails(int32_t tailSize){

  BamTools::BamReader reader;
  BamTools::BamAlignment al;
  if (!reader.Open(rawBamPath_)){
    std::cout << "Could not open input Bam file" << rawBamPath_ << std::endl;
    std::cout << "Existing run with non-sero status.." << std::endl;
    exit (EXIT_FAILURE);
  }

  reader.LocateIndex();
  if (!reader.HasIndex()){
    std::cout << "Index for" << rawBamPath_ << "could not be opened" << std::endl;
    std::cout << "Exiting run with non-sero status.." << std::endl;
    reader.Close();
    exit (EXIT_FAILURE);
  }

  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){    
    for(auto caIt = std::begin(cvIt->contigAlignments); caIt != std::end(cvIt->contigAlignments); ++caIt){
      BamTools::BamRegion region = BamTools::BamRegion(caIt->alignedContig.RefID, caIt->alignedContig.Position, caIt->alignedContig.RefID, caIt->alignedContig.GetEndPosition());
   
      std::cout << "setting region for coords : " << caIt->alignedContig.RefID << ", " <<  caIt->alignedContig.Position << ", " << caIt->alignedContig.RefID << ", " << caIt->alignedContig.GetEndPosition() << std::endl;
   
      if(!reader.SetRegion(region)) {
	std::cout << "could not set region for coords : " << caIt->alignedContig.RefID << ", " <<  caIt->alignedContig.Position << ", " << caIt->alignedContig.RefID << ", " << caIt->alignedContig.GetEndPosition() << std::endl;
	
      }
      
      while(reader.GetNextAlignment(al)){
	caIt->readsInRegion += 1;
	if(caIt->readsInRegion > 998){
	  break;
	}
	//std::cout << "count of reads in region is: " << caIt->readsInRegion << std::endl;
	polyA tail = {al, tailSize};
	if(tail.isTail()){
	  if(tail.isTailLeftBound()){
	    caIt->leftBoundTails.push_back(tail);
	  }
	  else {
	    caIt->rightBoundTails.push_back(tail);
	  }
	}
      }
      if(util::checkDoubleStranded(caIt->leftBoundTails)){
	caIt->doubleStranded = true;
	caIt->leftBound = true;
      }
      if(util::checkDoubleStranded(caIt->rightBoundTails)){
	caIt->doubleStranded = true;
	caIt->rightBound = true;
      }
    }
  }
  reader.Close();
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
	 c.alusHit.push_back(mi->seq[r->rid].name);
	 free(r->p);
	 // std::cout << "found alu hit for contig: " << ks->name.s << std::endl;
       }

       if(c.alusHit.size() > 0){
	 //std::cout << "found contig containing knonw alu" << std::endl;
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


KnownAlus::KnownAlus(std::string rawBamPath, std::string contigFastqPath, std::string contigBamPath, std::string aluFastaPath, std::string aluIndexPath, std::string refPath, std::string refIndexPath) :  rawBamPath_(rawBamPath), contigFastqPath_(contigFastqPath), contigBamPath_(contigBamPath), aluFastaPath_(aluFastaPath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath), stub_(util::baseName(rawBamPath)){
   
  contigVec_ = {};
  refData_ = {};

  std::cout << "[1/5]  Populating reference data for " << stub_ << std::endl;
  KnownAlus::populateRefData();
  std::cout << "[2/5]  Finding contigs containing known alus for " << stub_ << std::endl;
  KnownAlus::findContigsContainingKnownAlus();
  std::cout << "[3/5]  Pulling contig hit alignments for " << stub_ <<  std::endl;
  KnownAlus::pullContigAlignments();

  std::cout << "[4/5]  Finding reads containing polyA tails for " << stub_ << std::endl;
  KnownAlus::findReadsContainingPolyTails(10);


  KnownAlus::printContigVec();

  std::ofstream bed;
  std::string bs = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/out/" + stub_ + ".bed";
  bed.open(bs);

  writeBedPEHeader(bed);
  std::cout << "[5/5] Writing out results to bed file " << stub_  << ".bed" << std::endl;
  writeContigVecToBedPE(bed);

  bed.close();

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
  
  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    while(reader.GetNextAlignment(al)){
      if(cvIt->name.compare(al.Name)==0 and al.HasTag("SA")){
	std::cout << "Checking for peak and clip coord intersection for read: " << al.Name << std::endl;
	if(util::intersectPeaksAndClips(util::getPeaks(al), util::getLocalClipCoords(al))){
	  std::cout << "Found intersection between peak and clips" << std::endl;
	  
	  contigAlignment ca = {};
	  ca.alignedContig = al;
	  //auto peakVector = util::getPeaks(al);
	  cvIt->contigAlignments.push_back(ca);
	}
      }
    }
    reader.Rewind();
  }
  reader.Close();
}

