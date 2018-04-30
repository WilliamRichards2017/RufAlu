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
   if(c.supportingReads.size() > 0) {
   std::cout << "##########__PRINTING_CONTIG__############";
   std::cout << "name: " << c.name << std::endl;
   std::cout << "seq: " << c.seq << std::endl;
   std::cout << "qual: " << c.qual << std::endl;
   std::cout << "**** ALUS HIT ****" << std::endl;

   for (auto it = std::begin(c.alusHit); it != std::end(c.alusHit); ++it){
     std::cout << "    " << *it << std::endl;
   }

   std::cout << "**** contigAlignment ****" << std::endl;
   for (auto it2 = std::begin(c.contigAlignments); it2 != std::end(c.contigAlignments); ++it2){
     std::cout << "    " << it2->Name << std::endl;
   }

   std::cout << "**** contigAlignmentRegions ****" << std::endl;
   for(auto it3 = std::begin(c.contigAlignmentRegions); it3 != std::end(c.contigAlignmentRegions); ++it3){
     std::cout << "    left position: " << it3->LeftPosition << "   right position: " << it3->RightPosition << std::endl;
   }

   std::cout << "**** overlapingReads ****" << std::endl;
   for(auto it4 = std::begin(c.overlapingReads); it4 != std::end(c.overlapingReads); ++it4){
     for(auto it5 = std::begin(it4->window); it5 != std::end(it4->window); ++it5){
       std::cout << "   " <<  it5->Name << std::endl;
     }
   }

   std::cout << "**** supportingReads ****" << std::endl;
   for(auto it6 = std::begin(c.supportingReads); it6 != std::end(c.supportingReads); ++it6){
     std::cout << "    " <<  it6->Name << std::endl;
   }

   std::cout << "longest tail: " << c.longestTail << std::endl;
   }
 }

 void printContigVec(std::vector<contig> contigVec){
   for(auto it = std::begin(contigVec); it != std::end(contigVec); ++it){
     //TODO, use const_iterator to avoid temp
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

   //if(id > 23 or id < 0){
   //ret = "fuck";
     //return ret;
     //}
   std::cout << "id is " << id << std::endl;
   if(id == -1) {
     ret = "unmapped";
     return ret;
   }
   
   std::cout << "tying to look up ref id: " << id << std::endl;
   ret = refData_[id].RefName;
   return ret;
 }

 void writeBedPEHeader(std::ofstream &bed){
   bed << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' << "contig_name" << '\t' << "alu_hit" << '\t' << "num_hits" << '\t' << "longest_tail" << std::endl;
 }


 void KnownAlus::writeHitToBed(std::ofstream &bed, bedPELine * b){
   //std::cout << "Inside writeTobed" << std::endl;
   bed << b->chrom1 << '\t' << b->chrom1Start << '\t' << b->chrom1End << '\t' << b->chrom2 << '\t' << b->chrom2Start << '\t' << b->chrom2End << '\t'<< b->name_rufus_contig << '\t' << b->name_alu_hit << '\t' << b->score_numHits << '\t' << b->longestTail << std::endl;
 }

 void KnownAlus::findReadsContainingPolyATails(std::vector<contig>  contigs, std::string inputFile){
   //std::cout << "inside findReadsCOntainingPolyATails()" << std::endl;
   //std::cout << "number of contigs: " << contigs.size();
   std::ofstream bed;
   std::string bs = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/out/" + util::baseName(inputFile) + ".bed";
   bed.open(bs);
   writeBedPEHeader(bed);                                 

   for(auto it = std::begin(contigs); it != std::end(contigs); ++it){
     bedPELine * b = new bedPELine{};
     //std::cout << "number of overlaping reads in contig: " << it->overlapingReads.size();
     //std::cout << "number of alus hit: " << it->alusHit.size();
     b->name_rufus_contig = it->name;
     b->name_alu_hit = it->alusHit[0];
     for(auto ait = std::begin(it->alusHit); ait != std::end(it->alusHit); ++ait){
       //std::cout << "contig hit alu: " << *ait << std::endl;
     }
     uint32_t count = 0;
     uint32_t maxTail = 0;


     //std::cout << "size of overlaping reads vec is: " << it->overlapingReads.size();
     for(auto rIt = std::begin(it->overlapingReads); rIt != std::end(it->overlapingReads); ++rIt){
       //std::cout << "looping through overlaping reads vector" << std::endl;
       for(auto qIt = std::begin(rIt->window); qIt != std::end(rIt->window); ++qIt ) {
	 //std::cout << "checking for poly tail for read: " << qIt->Name << std::endl;
	 bool a = polyA::detectPolyATail(*qIt);
	 bool t = polyA::detectPolyTTail(*qIt);
	 uint32_t longestTail = polyA::longestTail(*qIt);

	 if(longestTail > maxTail){
	   maxTail=longestTail;
	 }

	 if(a || t){
	   b->score_numHits++;

	   if(count == 0){
	     b->chrom1 = getChromosomeFromRefID(qIt->RefID);
	     b->chrom1Start = qIt->Position;
	     b->chrom1End = qIt->GetEndPosition();
	   }

	   if(count == 1){
	     b->chrom2 = getChromosomeFromRefID(qIt->RefID);
	     b->chrom2Start = qIt->Position;
	     b->chrom2End = qIt->GetEndPosition();
	   }

	   ++count;

	   //std::cout << "found poly A/T Tail evidence supporting alu" << std::endl;
	   //std::cout <<  qIt->QueryBases << std::endl;
	   //std::string chrom = getChromosomeFromRefID(it->RefID);                                                                                                                      
	   //std::string chrom = "Null";
	   //std::cout << "supporting read was found at region: " << qIt->Position << ", " << qIt->GetEndPosition() << " RefID: " << qIt->RefID << " " << "Name: " << qIt->Name << std::endl;
	 }

       }

       b->longestTail = maxTail;
       it->longestTail = maxTail;
       if (b->score_numHits > 1){
	 writeHitToBed(bed, b);
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


       std::string qual = "";
       for(unsigned i = 0; i < std::string(ks->seq.s).length(); ++i){
	 qual+= '!';
       }

       //zero-initialize struct
       // zero initialization of dao                                                                                                                                                    
       bedPELine * b = new bedPELine{};
       contig  c = {};

       c.name = std::string(ks->name.s);
       c.seq = ks->seq.s;
       c.qual = qual;

       //std::cout << "nreg is: " << n_reg << std::endl;
       for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	 //std::cout << "traversing hits" << std::endl;
	 mm_reg1_t *r = &reg[j];

	 c.alusHit.push_back(mi->seq[r->rid].name);
	 //std::cout << "checking if aluHit[j] is populated: " << c.alusHit[j] << std::endl;

	 //assert(r->p); // with MM_F_CIGAR, this should not be NULL
	 //printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
	 //printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
	 //std::cout << "checking if alu name is: " << mi->seq[r->rid].name << std::endl;
	 for (i = 0; i < r->p->n_cigar; ++i) { // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
	   //printf("%d%c", r->p->cigar[i]>>4, "MIDSHN"[r->p->cigar[i]&0xf]);
	 }
	 //putchar('\n');
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
   mm_idx_reader_close(r); // close the index reader                                                                                                                                                                                           
   kseq_destroy(ks); // close the query file                                                                                                                                                                                                   
   gzclose(f);
 }


KnownAlus::KnownAlus(std::string contigFilePath, std::string contigBamPath, std::string mutationPath, const char * aluFilePath, const char * aluIndexPath, const char * refPath, const char * refIndexPath) : contigFilePath_(contigFilePath), contigBamPath_(contigBamPath), mutationPath_(mutationPath), aluFilePath_(aluFilePath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath), stub_("." + util::baseName(contigBamPath)){
   
  std::string contigsWithAlus = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted" + stub_ + ".bam";

  contigVec_ = {};
  refData_ = {};
  
  
  KnownAlus::populateRefData(mutationPath_);
  KnownAlus::findContigsContainingKnownAlus();
  
    //KnownAlus::alignContigsContainingKnownAlus(refIndexPath_);
  contigVec_ = KnownAlus::pullNamesWithHits(contigVec_, contigBamPath_);


  Intersect intersect{contigVec_, mutationPath_};
  contigVec_ = intersect.getContigVec();


  //std::cout << "finished intersection, now finding supporting reads with polyA tail" << std::endl;
  findReadsContainingPolyATails(contigVec_, mutationPath_);

  printContigVec(contigVec_);



}

KnownAlus::~KnownAlus(){
  //contigsContainingKnownAlus_->clear();
  //delete[] refData_;
  // delete contigFilePath_;
  //delete aluFilePath_;
  //delete aluIndexPath_;
  //delete refPath_;
  //delete refIndexPath_;
  //delete mutationPath_;
}

bool KnownAlus::checkIfNameInContigVec(BamTools::BamAlignment al){
  for(unsigned u = 0; u < contigVec_.size(); ++u){

    std::cout << "checking if " << al.Name << " == " << contigVec_[u].name << std::endl;
    if(contigVec_[u].name.compare(al.Name)){
      contigVec_[u].contigAlignments.push_back(al);
      std::cout << "pushing back" << al.Name << " to contigVec[" << u << "]" << std::endl;
      return true;
    }
  }
  return false;
}

std::vector<contig> KnownAlus::pullNamesWithHits(std::vector<contig> contigVec, std::string bamPath){
  std::cout << "inside pull name with hits" << std::endl;
  BamTools::BamReader reader;

  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }

  BamTools::BamAlignment al;
  while (reader.GetNextAlignment(al)){
    for(auto it = std::begin(contigVec); it != std::end(contigVec); ++it){
      if(it->name.compare(al.Name)==0){
	it->contigAlignments.push_back(al);
      }
    }
  }
  return contigVec;
}
