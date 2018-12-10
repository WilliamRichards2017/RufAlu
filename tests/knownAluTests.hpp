#ifndef RUFALU_TESTS_KNOWN_ALUS_HPP
#define RUFALU_TESTS_KNOWN_ALUS_HPP

#include "fastqParse.h"
#include "knownAlus.h"
#include "kseq.h"
#include "vcfWriter.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct input{
  input(std::string proband, std::string mother, std::string father) : proband_(proband), mother_(mother), father_(father){} 
  
  std::string proband_;
  std::string mother_;
  std::string father_;

  const std::string bamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + proband_ + ".bam";
  const std::string contigFilePath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + proband_ + ".bam.generator.V2.overlap.hashcount.fastq";
  const std::string mutationPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + proband_ + ".bam.generator.Mutations.fastq.bam";
  const std::string aluFilePath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta";
  const std::string aluIndexPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta.fai";
  const std::string refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const std::string refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";
  const std::string contigBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + proband_ + ".bam.generator.V2.overlap.hashcount.fastq.bam";
  const std::string fatherBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + mother_ + ".bam";
  const std::string motherBamPath = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/" + father_ + ".bam";
  const std::string fastaHackPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/bin/externals/fastahack/src/fastahack_project/bin/tools/fastahack";
  const std::string vcfOutPath = "/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/gtest.out";
  const std::vector<std::string> parentBams = {fatherBamPath, motherBamPath};
};

TEST(KnownAluTests, populateRefData){
  BamTools::BamReader reader;
  const char * bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam";
  ASSERT_TRUE(reader.Open(bam));

  BamTools::RefVector refVec = (reader.GetReferenceData());
  ASSERT_EQ(refVec.size(), 86);
}


TEST(KnownAluTests, p1348){

  input i = {std::to_string(1348), std::to_string(1333), std::to_string(1334)};

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  ASSERT_TRUE(reader.Open(i.bamPath));
  ASSERT_TRUE(reader.LocateIndex());

  reader.SetRegion(6, 110305760, 6, 110305762);

  while(reader.GetNextAlignment(al)){
    if(al.Name.compare("/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam") == 0){
      std::cerr << "Found contig: " << al.Name << std::endl;

      //contigAlignmnet ca = contigAlignment::contigAlignment(i.bamPath, i.parentBams, )
    }
  }

 

}


#endif // RUFALU_TESTS_KNOWN_ALUS_HPP

