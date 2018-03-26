#include <assert.h>
#include <stdexcept>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <zlib.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>

#include "fastqParse.h"
#include "knownAlus.h"
#include "kseq.h"
#include "minimap.h"
#include "polyATail.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

KSEQ_INIT(gzFile, gzread)

void exec(char const* cmd) {
  char buffer[512];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 128, pipe) != NULL)
	result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return;
}

void KnownAlus::findContigsContainingPolyATails(const char * inputFile){

  BamTools::BamReader reader;
  if (!reader.Open(inputFile)){
    std::cout << "Could not open input Bam file" << inputFile << std::endl;
    return;
  }

  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    bool b = polyA::detectPolyATail(al.QueryBases);
    if(b){
      std::cout << "found alu w/ poly a tail" << std::endl;
      std::cout <<  al.QueryBases << std::endl;
    }
    
  }
}

bool overlap(std::pair<int, int> a, std::vector<std::pair<int, int> > b){
  for (auto it = std::begin(b); it != std::end(b); ++it){
    if((a.first >= it->first && a.first <= it->second) || (a.second >= it->first && a.second <= it->second)){
      return true;
    }
  }
  return false;
}

std::vector<BamTools::BamAlignment> intersectBams(const char * a, const char * b){
  std::vector<BamTools::BamAlignment> intersection;


  std::vector<std::pair<int, int> > aCoords;
  BamTools::BamReader reader;
  if (!reader.Open(a)){
    std::cout << "Could not open the following input bamfile: " << a << std::endl;
    return intersection;
  }

  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    std::pair<int,int> coords = std::make_pair(al.Position, al.GetEndPosition());
    aCoords.push_back(coords);
    //std::cout << "pushing back coords " << coords.first << ", " << coords.second << std::endl;
  }

  if (!reader.Open(b)){
    std::cout << "Could not open the following input bamfile: " << b << std::endl;
    return intersection;
  }
  
  BamTools::BamAlignment bl;
  
  while(reader.GetNextAlignment(bl)){
    std::pair<int, int> coords = std::make_pair(bl.Position, bl.GetEndPosition());
    if (overlap(coords, aCoords)){
      intersection.push_back(bl);
      std::cout << "found overlap at coords " << coords.first << ", " << coords.second << std::endl;
    }
  }

  return intersection;
}


const char *contigsToFastq(std::vector<fastqRead *> *contigs, const char * outFile){

  std::ofstream out;
  out.open(outFile);

  for(auto it = std::begin(*contigs); it != std::end(*contigs); ++it){

    //out << (*it)->id << std::endl;
    out << '@' << (*it)->name << std::endl;
    out << (*it)->seq << std::endl;
    out << '+' << std::endl;
    out << (*it)->qual << std::endl;
    std::cout << "name: " << (*it)->name <<  std::endl;
    std::cout << "sequence: " << (*it)->seq <<  std::endl;
    std::cout << "quality: " << (*it)->qual <<  std::endl;
        
  }  
  
  return outFile;
}

void KnownAlus::findContigsContainingKnownAlus()
{

  mm_idxopt_t iopt;
  mm_mapopt_t mopt;
  int n_threads = 3;


  std::cout << "reading in contig fasta file " << contigFilePath_ << std::endl;


  mm_verbose = 3; // print to std out
  mm_set_opt(0, &iopt, &mopt); //initialize alignment parameters to default
  mopt.flag |= MM_F_CIGAR; // perform alignment                                                                                                                                                                                               

  gzFile f = gzopen(contigFilePath_, "r");
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

      fastqRead *f = new fastqRead(std::string(ks->name.s), ks->seq.s, qual);

      contigsContainingKnownAlus_->push_back(f);


      
      //std::cout << "nreg is: " << n_reg << std::endl;
      for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
	//std::cout << "traversing hits" << std::endl;
        mm_reg1_t *r = &reg[j];
        assert(r->p); // with MM_F_CIGAR, this should not be NULL

	printf("%s\t%d\t%d\t%d\t%c\t", ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
        printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:", mi->seq[r->rid].name, mi->seq[r->rid].len, r->rs, r->re, r->mlen, r->blen, r->mapq);
        for (i = 0; i < r->p->n_cigar; ++i) { // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
	  printf("%d%c", r->p->cigar[i]>>4, "MIDSHN"[r->p->cigar[i]&0xf]);
	}
        putchar('\n');
        free(r->p);
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


KnownAlus::KnownAlus(const char * contigFilePath, const char * aluFilePath, const char * aluIndexPath, const char * refPath, const char * refIndexPath) : contigFilePath_(contigFilePath), aluFilePath_(aluFilePath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath){
  contigsContainingKnownAlus_ = new std::vector<fastqRead *>;
  //KnownAlus::findContigsContainingKnownAlus();
  //KnownAlus::alignContigsContainingKnownAlus(refIndexPath_);
  //KnownAlus::recoverPolyATails();
  //findContigsContainingPolyATails("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/alu_intersect.bam");
  intersectBams("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted.bam", "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/Family1.child.bam.generator.Mutations.fastq.bam");
}

KnownAlus::~KnownAlus(){
  contigsContainingKnownAlus_->clear();
}

void KnownAlus::mapContigsToRef(const char * contigs){
  std::string cmd = "";
  cmd+= "cd /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/externals/minimap2/src/minimap2_project ; ./minimap2 -a ";
  cmd+= refPath_;
  cmd += " /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/";
  cmd+= contigs;
  cmd+= " > /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sam";
  std::cout << "executing command " << cmd << std::endl;
  exec(cmd.c_str());
}

void KnownAlus::recoverPolyATails(){

  std::string cmd = "";
  cmd+= "bedtools intersect -a /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/Family1.child.bam.generator.Mutations.fastq.bam -b /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted.bam ";
  cmd += " > /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/alu_intersect.bam";
  
  std::cout << "executing command " << cmd << std::endl;
  exec(cmd.c_str());
  

}

std::vector<fastqRead *> * KnownAlus::getContigsContainingKnownAlus(){
  return contigsContainingKnownAlus_;
}

void KnownAlus::alignContigsContainingKnownAlus(const char * refPath){

  const char * fastq = contigsToFastq(contigsContainingKnownAlus_, "contigs.fastq"); 
  KnownAlus::mapContigsToRef(fastq);
}
