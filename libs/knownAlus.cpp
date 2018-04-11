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

void KnownAlus::populateRefData(const char * bamPath){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }
  refData_ = &(reader.GetReferenceData());
}

std::string KnownAlus::getChromosomeFromRefID(int32_t id){
  std::cout << "id is " << id << std::endl;
  if(id == -1) {
    return "unmapped";
  }
  return (*refData_)[id].RefName;
}

void KnownAlus::findReadsContainingPolyATails(std::vector<BamTools::BamAlignment>  contigs, const char * inputFile){
  //std::cout << "inside findReadsContainingPlolyATails()" << std::endl;
  for(auto it = std::begin(contigs); it != std::end(contigs); ++it){
    //std::cout << "looping through contigs" << std::endl;
    bool b = polyA::detectPolyATail(it->QueryBases);
    bool t = polyA::detectPolyTTail(it->QueryBases);
    if(b || t){
      std::cout << "found poly A/T Tail evidence supporting alu" << std::endl;
      std::cout <<  it->QueryBases << std::endl;
      //std::string chrom = getChromosomeFromRefID(it->RefID);
      std::string chrom = "Null";
      std::cout << "supporting read was found at region: " << it->Position << ", " << it->GetEndPosition() << " RefID: " << it->RefID << " " << "Name: " << it->Name << std::endl;
    }
  }
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
      contigsContainingKnownAlus_->push_back(*f);
      
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


KnownAlus::KnownAlus(const char * contigFilePath, const char * contigBamPath, const char * mutationPath, const char * aluFilePath, const char * aluIndexPath, const char * refPath, const char * refIndexPath) : contigFilePath_(contigFilePath), contigBamPath_(contigBamPath), mutationPath_(mutationPath), aluFilePath_(aluFilePath), aluIndexPath_(aluIndexPath), refPath_(refPath), refIndexPath_(refIndexPath){
  contigsContainingKnownAlus_ = new std::vector<fastqRead>;

  //const char * rootDir = util::getRootDirectory(std::string(aluFilePath));
  //std::cout << "RUFUS root path is: " << rootDir << std::endl;
  const char * contigsWithAlus = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted.bam";

  std::cout << "Contig bam path is: " << contigBamPath_ << std::endl;

  KnownAlus::populateRefData(mutationPath_);
  //KnownAlus::findContigsContainingKnownAlus();
  //KnownAlus::alignContigsContainingKnownAlus(refIndexPath_);
  std::cout << "finished aligning contings to known alus, now intersecting bams" << std::endl;

  const char * contigsWithAluHits = Intersect::getContigHits(contigBamPath_);
 
  Intersect intersect{contigsWithAluHits, mutationPath_};
  std::vector<BamTools::BamAlignment> reads = intersect.getIntersection();

  for(auto it = std::begin(reads); it != std::end(reads); ++it){
    std::cout << it->QueryBases << std::endl;
  }

  std::cout << "size of intersection is: " << reads.size() << std::endl;

  std::cout << "finished intersection, now finding supporting reads with polyA tail" << std::endl;
  findReadsContainingPolyATails(reads, mutationPath_);

}

KnownAlus::~KnownAlus(){
  //contigsContainingKnownAlus_->clear();
  delete[] contigsContainingKnownAlus_;
  delete[] refData_;
  delete contigFilePath_;
  delete aluFilePath_;
  delete aluIndexPath_;
  delete refPath_;
  delete refIndexPath_;
  delete mutationPath_;
}

void KnownAlus::mapContigsToRef(const char * contigs){

  std::string basepath  = util::baseName(std::string(contigs));
  std::cout << "Base path is: " << basepath << std::endl;

  std::string cmd = "";
  std::string sort = "";
  util::exec("chmod +x /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/bamtools/src/bamtools_project/bin/bamtools");
  cmd+= "cd /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/externals/minimap2/src/minimap2_project ; ./minimap2 -ax sr ";
  cmd+= refPath_;
  cmd += " /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/";
  cmd+= contigs;
  cmd+= " > /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sam";
  //std::cout << "executing command " << cmd << std::endl;
  util::exec(cmd.c_str());
  //convert sam to bam
  util::exec("samtools view -Sb /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sam > /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.bam");
  //Sort bam file by position
  util::exec("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/bamtools/src/bamtools_project/bin/bamtools sort -in /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.bam -out /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted.bam");
  //Index sorted bamfile
  util::exec("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/externals/bamtools/src/bamtools_project/bin/bamtools index -in /uufs/chpc.utah.edu/common/home/u0401321/RufAlu/data/contigs-with-alus.sorted.bam");
  
}
std::vector<fastqRead> * KnownAlus::getContigsContainingKnownAlus(){
  return contigsContainingKnownAlus_;
}

void KnownAlus::alignContigsContainingKnownAlus(const char * refPath){

  const char * fastq = util::contigsToFastq(contigsContainingKnownAlus_, "contigs.fastq"); 
  KnownAlus::mapContigsToRef(fastq);

}
