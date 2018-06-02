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


void printContig(contig c){

}

void printContigVec(std::vector<contig> contigVec){
  for(auto it = std::begin(contigVec); it != std::end(contigVec); ++it){
    printContig(*it);
  }
}

void KnownAlus::populateRefData(std::string bamPath){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = reader.GetReferenceData();
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
  bed << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' 
      << "contig_name" << '\t' << "alu_hit" << '\t' << "num_hits" << '\t' << "longest_tail" <<  '\t' << "both_strands" <<std::endl;
}

void KnownAlus::writeContigVecToBedPE(std::ofstream &bed){
  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    for(auto caIt = std::begin(cvIt->contigAlignments); caIt != std::end(cvIt->contigAlignments); ++caIt){
      if(caIt->primaryAlignment) {
	bed << getChromosomeFromRefID(caIt->alignedContig.RefID) << '\t'<< caIt->alignedContig.Position << '\t' << caIt->alignedContig.GetEndPosition() 
	    << '\t' << '-' << '\t' << '-' << '\t' << '-' << '\t'<< cvIt->name << '\t' << cvIt->alusHit[0] << '\t' << caIt->supportingReads.size() << '\t' 
	    << caIt->longestTail << '\t' << caIt->doubleStranded << std::endl;
      } else {
        bed << '-' << '\t' << '-' << '\t' << '-' << getChromosomeFromRefID(caIt->alignedContig.RefID) << '\t'<< caIt->alignedContig.Position<< '\t'
		<< caIt->alignedContig.GetEndPosition() << '\t' << '-' << '\t' << '-' << '\t' << '-' << '\t'<< cvIt->name << '\t' << cvIt->alusHit[0] << '\t' << caIt->supportingReads.size() << '\t'
		<< caIt->longestTail << '\t' << caIt->doubleStranded << std::endl;
      }
    }
  }
}

void KnownAlus::findReadsContainingPolyTails(uint32_t tailSize){
  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){    
    for(auto caIt = std::begin(cvIt->contigAlignments); cvIt != std::end(cvIt->contigAlignments); ++caIt){
      BamTools::BamReader reader;
      if (!reader.Open(bamPath_)){
	std::cout << "Could not open input Bam file" << bamPath << std::endl;
	exit (EXIT_FAILURE);
      }
      reader.setRegion(caIt->alignedRegion);
      
      BamTools::BamAlignment al;
      while(reader.GetNextAlignment(al)){
	bool tail = polyA::detectPolyTailClips(al, tailSize);
	if(tail){
	  if(al->IsReverseStrand()){
	    caIt->reverseStrand = true;
	  }
	  else {
	    caIt->forwardStrand = true;
	  }
	}
      }
      if(caIt->forwardStrand and caIt->reverseStrand){
	caIt->doubleStranded = true;
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

   gzFile f = gzopen(contigFilePath_.c_str(), "r");
   assert(f);
   kseq_t *ks = kseq_init(f);

   // open index reader
   mm_idx_reader_t *r = mm_idx_reader_open(aluFilePath_, &iopt, 0);
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
       
       for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	 mm_reg1_t *r = &reg[j];
	 c.alusHit.push_back(mi->seq[r->rid].name);
	 free(r->p);
       }

       if(c.alusHit.size() > 0){
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


 KnownAlus::KnownAlus(std::string bamPath, std::string contigFilePath, std::string contigBamPath, std::string mutationPath, const char * aluFilePath, const char * aluIndexPath, const char * refPath, const char * refIndexPath) :  bamPath_(bamPath), contigFilePath_(contigFilePath), contigBamPath_(contigBamPath), mutationPath_(mutationPath), aluFilePath_(aluFilePath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath), stub_("." + util::baseName(contigBamPath)){
   
  std::string contigsWithAlus = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted" + stub_ + ".bam";

  contigVec_ = {};
  refData_ = {};
  
  KnownAlus::populateRefData(mutationPath_);
  KnownAlus::findContigsContainingKnownAlus();
  KnownAlus::pullNamesWithHits(contigBamPath_);
  Intersect intersect{contigVec_, mutationPath_};
  contigVec_ = intersect.getContigVec();

  findReadsContainingPolyTails(10);

  std::ofstream bed;
  std::string bs = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/out/" + util::baseName(inputFile) + ".bed";
  bed.open(bs);

  writeBedPEHeader(bed);
  writeContigVecToBedPE(bed);

  bed.close();
  
  printContigVec(contigVec_);



}

KnownAlus::~KnownAlus(){
  delete aluFilePath_;
  delete aluIndexPath_;
  delete refPath_;
  delete refIndexPath_;
}

void KnownAlus::pullNamesWithHits(std::string bamPath){

  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file " << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }

  for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
    BamTools::BamReader reader;

    if (!reader.Open(bamPath)){
      std::cout << "Could not open input Bam file " << bamPath << std::endl;
      exit (EXIT_FAILURE);
    }

    BamTools::BamAlignment al;
    while(reader.GetNextAlignment(al)){
      if(cvIt->name.compare(al.Name)==0){
	cvIt->contigAlignments.push_back(al);
      }
    }
  }
}
  
