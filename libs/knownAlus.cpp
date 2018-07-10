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
  if(ca.alignedRegion.LeftPosition != -1 and ca.alignedRegion.RightPosition != -1){
    return true;
  }
  return false;
}

void printContigAlignment(contigAlignment & ca){
  if(ca.alignedContig.Name.compare("NODE_8800.bam.generator.V2_372_L151_D6:6:0::MH0") == 0){
    std::cout << std::endl;
    std::cout << "printing contig alignment: " << ca.alignedContig.Name;
    std::cout << "ca.leftBoundHeads.size(): " << ca.leftBoundHeads.size();
    for(auto head : ca.leftBoundHeads){
      std::cout << "head name is: " << head.al_.Name << std::endl;
      std::cout << "head seq is: " << head.al_.QueryBases << std::endl;
    }
    for(auto tail : ca.leftBoundTails){
      std::cout << "tail name is: " << tail.al_.Name << std::endl;
      std::cout << "tail seq is: " << tail.al_.QueryBases << std::endl;
      
    }
    std::cout << std::endl;
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

void KnownAlus::writeBedPEHeader(std::ofstream &bed){
  bed << "chrom" << '\t' << "chromStart" << '\t' << "chromEnd" << '\t' << "chrom" << '\t'<< "contig_name" << '\t' << "alu_hit" << '\t' << "primary_alignment" << '\t' << "numReadsInRegion" << '\t' << "numPolyATails" <<  '\t' << "both_strands" << '\t' << "longest_tail"  << std::endl;
}

bool KnownAlus::bedFilter(contigAlignment & ca) {
  if(ca.readsInRegion < 400 and std::max(ca.rightBoundTails.size(), ca.leftBoundTails.size()) > 1){
    return true;
      if(ca.tailLeftBound){
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


void KnownAlus::writeToVCF(std::string & vcfFile){

  std::ofstream vcfStream;
  vcfStream.open(vcfFile, std::ios::app);
  vcfWriter::writeVCFHeader(vcfStream, stub_);
  KnownAlus::writeContigVecToVCF(vcfStream);
  vcfStream.close();
}

void KnownAlus::writeContigVecToVCF(std::ofstream & vcf){
  //for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
  for(auto c : contigVec_){
    for(auto ca : c.contigAlignments){
      vcfWriter writer = {ca, vcf, stub_};
      if(ca.denovoVec_.size() > 0){
	std::cout << "CA is denovo in KnownAlus " << std::endl;
      }
      if(writer.vcfFilter()){
	writer.writeVCFLine();
      }
    }
  }
}

void KnownAlus::writeContigVecToBedPE(std::ofstream &bed){
  /*for(auto cvIt = std::begin(contigVec_); cvIt != std::end(contigVec_); ++cvIt){
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
	  bed << '\t' << "tail not left or right bound..." << '\t' << std::max(caIt->rightBoundTails.size(), caIt->leftBoundTails.size()) << std::endl;
	}
      }
    }
    }*/
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
      //std::cout << "setting region for coords : " << caIt->alignedContig.RefID << ", " <<  caIt->alignedContig.Position << ", " << caIt->alignedContig.RefID << ", " << caIt->alignedContig.GetEndPosition() << std::endl;
   
      if(!reader.SetRegion(region)) {
	std::cout << "could not set region for coords : " << caIt->alignedContig.RefID << ", " <<  caIt->alignedContig.Position << ", " << caIt->alignedContig.RefID << ", " << caIt->alignedContig.GetEndPosition() << std::endl;
	
      }
      while(reader.GetNextAlignment(al)){
	caIt->readsInRegion += 1;
	if(caIt->readsInRegion > 200){
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
	caIt->tailLeftBoundDS = true;
      }
      if(util::checkDoubleStranded(caIt->rightBoundTails)){
	caIt->tailRightBoundDS = true;
      }
      if(caIt->leftBoundTails.size() > 1){
	caIt->tailLeftBound = true;
      }
      if(caIt->rightBoundTails.size() > 1){
	caIt->tailRightBound = true;
      }
      
    }
  }
  reader.Close();
}


void KnownAlus::findDenovoEvidence(){
  
  for(auto & c : contigVec_){
    for(auto & ca : c.contigAlignments){
      denovoEvidence de = {util::getClipSeqs(ca.alignedContig)[0], ca.alignedRegion,  parentBams_};

      ca.denovoVec_.push_back(de);
      std::cout << "is de denovo ? : " << de.isDenovo() << std::endl;
      ca.isDenovo = de.isDenovo();
      
    }
  }
}

void KnownAlus::findReadsContainingHeads(){
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

  for(auto cIt = std::begin(contigVec_); cIt != std::end(contigVec_); ++cIt){
    for(auto caIt = std::begin(cIt->contigAlignments); caIt != std::end(cIt->contigAlignments); ++caIt){
      if(caIt->readsInRegion < 201){
	
	BamTools::BamRegion region = BamTools::BamRegion(caIt->alignedContig.RefID, caIt->alignedContig.Position, caIt->alignedContig.RefID, caIt->alignedContig.GetEndPosition());
	if(!reader.SetRegion(region)) {
	  std::cout << "could not set region for coords : " << caIt->alignedContig.RefID << ", " <<  caIt->alignedContig.Position << ", " 
		    << caIt->alignedContig.RefID << ", " << caIt->alignedContig.GetEndPosition() << std::endl;
	}
	BamTools::BamAlignment al;
	while(reader.GetNextAlignment(al)){
	  aluHead head = {util::getClipSeqs(caIt->alignedContig)[0], al, 10};
	  if(head.isHead()){
	    caIt->leftBoundHeads.push_back(head);
	  }
	}
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

void KnownAlus::writeToBed(std::string & bedFile){
  std::ofstream bed;
  bed.open(bedFile);
  writeBedPEHeader(bed);

  writeContigVecToBedPE(bed);

  bed.close();
}


std::vector<contigAlignment>  KnownAlus::findParentContigAlignments(const BamTools::BamAlignment & al, const BamTools::BamRegion & region, const std::vector<std::string> & parentBamPaths){
  
  std::vector<contigAlignment> parentContigAlignments;
  for(auto bp : parentBamPaths){
    contigAlignment c;
    c.bamPath = bp;
    c.alignedRegion = region;
    c.alignedContig = al;
    parentContigAlignments.push_back(c);
  }
  return parentContigAlignments;
}


std::vector<contigAlignment> KnownAlus::populateParentContigAlignments(std::vector<contigAlignment>  parentContigAlignments){
  for (auto pCA : parentContigAlignments){

    std::cout << "populating parent contig" << std::endl;

    BamTools::BamReader reader;

    if (!reader.Open(pCA.bamPath)){
      std::cerr << "Could not open input Bam file" << pCA.bamPath << std::endl;
      std::cerr << "Existing run with non-sero status.." << std::endl;
      exit (EXIT_FAILURE);
    }

    reader.LocateIndex();
    if (!reader.HasIndex()){
      std::cerr << "Index for" << pCA.bamPath << "could not be opened" << std::endl;
      std::cerr << "Exiting run with non-sero status.." << std::endl;
      reader.Close();
      exit (EXIT_FAILURE);
    }

    if(!reader.SetRegion(pCA.alignedRegion)){
      std::cerr << "Could not set region for coords: " << pCA.alignedRegion.LeftRefID << ", " << pCA.alignedRegion.LeftPosition << ", " << pCA.alignedRegion.RightRefID << ", " << pCA.alignedRegion.RightPosition << std::endl;
      std::cerr << "Exiting run with non-sero status.." << std::endl;
      reader.Close();
      exit (EXIT_FAILURE);
    }


    std::cout << "setting region coords to be: " << pCA.alignedRegion.LeftRefID << ", " << pCA.alignedRegion.LeftPosition << ", " << pCA.alignedRegion.RightRefID << ", " << pCA.alignedRegion.RightPosition << std::endl;

    BamTools::BamAlignment al;
    while(reader.GetNextAlignment(al)){
      std::cout << "checking alignment for heads and tails" << std::endl;
      aluHead head = {(util::getClipSeqs(pCA.alignedContig))[0], al,10};
      if(head.isHead()){
	std::cout << "found head in parent" << std::endl;
	pCA.leftBoundHeads.push_back(head);
      }
      polyA tail = {al, 10};
      if(tail.isTail()){
	if(tail.isTailLeftBound()){
	std::cout << "found tail in parent" << std::endl;
	  pCA.leftBoundTails.push_back(tail);
	}
	else {
	std::cout << "found tail in parent" << std::endl;
	  pCA.rightBoundTails.push_back(tail);
	}
      }
    }
    parentContigAlignments.push_back(pCA);
  reader.Close();
  }
  return parentContigAlignments;
}


void KnownAlus::flagAllDenovos(const std::vector<std::string> & parentBams){
  for(auto c :contigVec_){
    for(auto & ca : c.contigAlignments){
      if(debugPrintFilter(ca)){
	
	
	std::cout << "checking if alu hit  is denovo" << std::endl;
	std::vector<contigAlignment> caVec = findParentContigAlignments(ca.alignedContig, ca.alignedRegion, parentBams);
	caVec = KnownAlus::populateParentContigAlignments(caVec);
	ca.isDenovo = KnownAlus::isDenovo(caVec);
      }
    }
  }
}



const bool KnownAlus::isDenovo(const std::vector<contigAlignment> & parentContigAlignments){
  for(auto pCA : parentContigAlignments){
    std::cout << "sizes are " << pCA.leftBoundHeads.size() << ", " << pCA.leftBoundTails.size() << ", " << pCA.rightBoundTails.size();
    if(pCA.leftBoundHeads.size() > 0 || pCA.leftBoundTails.size() > 0 || pCA.rightBoundTails.size() > 0){
      std::cout << "returning non-denovo" << std::endl;
      return false;
    }
  }
  std::cout << "returning inherited" << std::endl;
  return true;
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

  std::cerr << "[4/7]  Finding reads containing polyA tails for " << stub_ << std::endl;
  KnownAlus::findReadsContainingPolyTails(9);

  std::cerr << "[5/7]  Finding reads containing alu heads tails for " << stub_ << std::endl;
  KnownAlus::findReadsContainingHeads();

  std::cerr << "[6/7] Flaging denovos for " << stub_ << std::endl; 
  KnownAlus::findDenovoEvidence();

 

  //std::cout << "[5/5] Writing out results to bed file " << stub_  << ".bed" << std::endl;
  //KnownAlus::writeToBed(prefix_ + stub_ + "bed");


  std::cerr << "[7/7] Writing out results to vcf file " << vcfOutPath << std::endl;
  KnownAlus::writeToVCF(vcfOutPath);

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
	//std::cout << "Checking for peak and clip coord intersection for read: " << al.Name << std::endl;
	//util::printCigar(al.CigarData);
	//std::cout << al.Qualities << std::endl;
	clipCoords clipPeak = util::intersectPeaksAndClips(util::getPeaks(al), util::getLocalClipCoords(al));
	if(clipPeak.clipStart != -1){
	  //std::cout << "Found intersection between peak and clips for " << cvIt->name << std::endl;
	  contigAlignment ca;
	  ca.clipCoords_ = clipPeak;
	  //ca.aluHit = cvIt->alusHit[0].first;
	  ca.aluHit = util::getHighestQualityAluHit(cvIt->alusHit);
	  ca.alignedContig = al;
	  ca.chrom = getChromosomeFromRefID(ca.alignedContig.RefID);
	  
	  ca.alignedRegion = {ca.alignedContig.RefID, ca.alignedContig.Position, ca.alignedContig.RefID, ca.alignedContig.GetEndPosition()};
	  
	  std::vector<int32_t> peakVector = util::getPeakVector(al);
	  auto it =  std::max_element(peakVector.begin(), peakVector.end());
	  ca.maxHash = peakVector[std::distance(peakVector.begin(), it)];
	  cvIt->contigAlignments.push_back(ca);
	}
      }
    }
    reader.Rewind();
  }
  reader.Close();
}
