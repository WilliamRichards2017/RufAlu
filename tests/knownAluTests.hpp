#ifndef RUFALU_TESTS_KNOWN_ALUS_HPP
#define RUFALU_TESTS_KNOWN_ALUS_HPP

#include "knownAlus.h"
#include "fastqParse.h"
#include "kseq.h"


#include "api/BamMultiReader.h"
#include "api/BamWriter.h"


TEST(KnownAluTests, populateRefData){
  BamTools::BamReader reader;
  const char * bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam";
  ASSERT_TRUE(reader.Open(bam));

  BamTools::RefVector refVec = (reader.GetReferenceData());
  ASSERT_EQ(refVec.size(), 86);
}

TEST(KnownAluTests, knownAlusIntegrationTest){

  static const char * contigFilePath = std::string("/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.V2.overlap.hashcount.fastq").c_str();
  static const char * mutationPath = std::string("/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam").c_str();
  static const char * aluFilePath = std::string("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/test_data/primate_non-LTR_Retrotransposon.fasta").c_str();
  static const char * aluIndexPath = std::string("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/test_data/primate_non-LTR_Retrotransposon.fasta.fai").c_str();
  static const char * refPath = std::string("/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa").c_str();
  static const char * refIndexPath = std::string("/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai").c_str();
  static const char * contigBamPath = std::string("/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.V2.overlap.hashcount.fastq.bam").c_str();

  KnownAlus *knownAlus = new KnownAlus(contigFilePath, contigBamPath, mutationPath, aluFilePath, aluIndexPath, refPath, refIndexPath);
  std::cout << "found " << knownAlus->getContigsContainingKnownAlus()->size() << " contig hits" << std::endl;

}


#endif // RUFALU_TESTS_KNOWN_ALUS_HPP
