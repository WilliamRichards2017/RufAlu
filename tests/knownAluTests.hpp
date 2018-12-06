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

  const std::string bamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam";
  const std::string contigFilePath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.V2.overlap.hashcount.fastq";
  const std::string mutationPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam";
  const std::string aluFilePath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta";
  const std::string aluIndexPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta.fai";
  const std::string refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const std::string refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";
  const std::string contigBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.V2.overlap.hashcount.fastq.bam";
  const std::string fatherBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1333.bam";
  const std::string motherBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1334.bam";
  const std::string fastaHackPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/bin/externals/fastahack/src/fastahack_project/bin/tools/fastahack";
  const std::string vcfOutPath = "~/RufAlu/bin/gtest.out";

  std::vector<std::string> parentBamPaths = {fatherBamPath, motherBamPath};

  auto knownAlus = new KnownAlus(bamPath, contigFilePath, contigBamPath, aluFilePath, aluIndexPath, refPath, refIndexPath, vcfOutPath, parentBamPaths, fastaHackPath);



  auto contigVec = knownAlus->getContigVec();

  std::cerr << "contigVec.size() is: " << contigVec.size() << std::endl;
  
  ASSERT_GT(contigVec.size(), 0);

  delete knownAlus;
}


#endif // RUFALU_TESTS_KNOWN_ALUS_HPP

