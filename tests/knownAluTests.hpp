#ifndef RUFALU_TESTS_KNOWN_ALUS_HPP
#define RUFALU_TESTS_KNOWN_ALUS_HPP

#include "fastqParse.h"
#include "knownAlus.h"
#include "kseq.h"
#include "vcfWriter.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

struct input{
  input(std::string prefix, std::string proband, std::string mother, std::string father) : prefix_(prefix), proband_(proband), mother_(mother), father_(father){} 
  
  std::string proband_;
  std::string mother_;
  std::string father_;
  std::string prefix_;

  const std::string bamPath = prefix_ + proband_ + ".bam";
  const std::string contigFilePath = prefix_ + proband_ + ".bam.generator.V2.overlap.hashcount.fastq";
  const std::string mutationPath = prefix_ + proband_ + ".bam.generator.Mutations.fastq.bam";
  const std::string aluFilePath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta";
  const std::string aluIndexPath = "/uufs/chpc.utah.edu/common/home/u0401321/RUFUS/resources/primate_non-LTR_Retrotransposon.fasta.fai";
  const std::string refPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
  const std::string refIndexPath = "/uufs/chpc.utah.edu/common/home/marth-ucgdstor/resources/references/human/GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.fai";
  const std::string contigBamPath = prefix_ + proband_ + ".bam.generator.V2.overlap.hashcount.fastq.bam";
  const std::string fatherBamPath = prefix_ + mother_ + ".bam";
  const std::string motherBamPath = prefix_ + father_ + ".bam";
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

  input i = {"/scratch/ucgd/lustre/u0991464/Projects/CEPH.1kg.cut0.5.v2/", std::to_string(1348), std::to_string(1333), std::to_string(1334)};

  BamTools::BamReader reader;
  BamTools::BamAlignment al;
  std::pair<std::string, int32_t> aluHit = std::make_pair("AluYa5", 60);

  ASSERT_TRUE(reader.Open(i.contigBamPath));
  ASSERT_TRUE(reader.LocateIndex());

  BamTools::BamRegion region = {6, 110305000, 6, 110306000};

  std::fstream streamy;

  streamy.open("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/testy.vcf");

  reader.SetRegion(region);


  bool contigFound = false;

  while(reader.GetNextAlignment(al)){
    std::cerr << "Found contig: " << al.Name << std::endl;

    if(al.Name.compare("NODE_1348.bam.generator.V2_102_L185_D8:6:2::MH0") == 0){
      
      contigFound = true;

      contigAlignment ca = contigAlignment(i.bamPath, i.parentBams, aluHit, al, "7", region, i.refPath, i.fastaHackPath);
      vcfWriter w = vcfWriter(ca, streamy, i.bamPath);
      vcfLine l = w.getVCFLine();

      ASSERT_EQ(l.CHROM, "7");
      ASSERT_EQ(l.POS, 110305761);
      ASSERT_STREQ(l.ID.c_str(), "ME-DeNovo");
      ASSERT_STREQ(l.REF.c_str(), "N");
      ASSERT_STREQ(l.ALT.c_str(), "INS:ME:AluYa5");
      ASSERT_EQ(l.FILTER.HDS and l.FILTER.TDS, true);

      EXPECT_EQ(l.QUAL, 45);
      EXPECT_EQ(l.INFO.NR, 84);
      EXPECT_EQ(l.INFO.NT, 6);
      EXPECT_EQ(l.INFO.NH, 11);
      
      ASSERT_EQ(l.INFO.LT, 72);
      ASSERT_STREQ(l.INFO.cigar.c_str(), "M140S45");

      ASSERT_EQ(l.INFO.probandGT.DP, 50);
      ASSERT_EQ(l.INFO.probandGT.RO, 20);
      ASSERT_EQ(l.INFO.probandGT.AO, 30);
      ASSERT_EQ(l.INFO.probandGT.genotype, std::make_pair(1,0));

      ASSERT_EQ(l.INFO.parentGTs[0].DP, 30);
      ASSERT_EQ(l.INFO.parentGTs[0].RO, 30);
      ASSERT_EQ(l.INFO.parentGTs[0].AO, 0);
      ASSERT_EQ(l.INFO.parentGTs[0].genotype, std::make_pair(0,0));

      ASSERT_EQ(l.INFO.parentGTs[1].DP, 26);
      ASSERT_EQ(l.INFO.parentGTs[1].RO, 26);
      ASSERT_EQ(l.INFO.parentGTs[1].AO, 0);
      ASSERT_EQ(l.INFO.parentGTs[1].genotype, std::make_pair(0,0));

    }
  }

  ASSERT_EQ(contigFound, true);
  streamy.close();

}



TEST(KnownAluTests, p2788){

  input i = {"/scratch/ucgd/lustre/u0991464/Projects/CEPH.1kg.cut0.5.v2/", std::to_string(2788), std::to_string(8125), std::to_string(8126)};

  BamTools::BamReader reader;
  BamTools::BamAlignment al;

  std::pair<std::string, int32_t> aluHit = std::make_pair("SVA_B", 0);

  ASSERT_TRUE(reader.Open(i.contigBamPath));
  ASSERT_TRUE(reader.LocateIndex());

  BamTools::BamRegion region = {10, 9072830, 6, 9072840};

  std::fstream streamy;

  streamy.open("/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/bin/testy.vcf");

  reader.SetRegion(region);


  bool contigFound = false;

  while(reader.GetNextAlignment(al)){
    std::cerr << "Found contig: " << al.Name << std::endl;

    if(al.Name.compare("NODE_2788.bam.generator.V2_736_L261_D14:6:8::MH0") == 0){
      
      contigFound = true;

      contigAlignment ca = contigAlignment(i.bamPath, i.parentBams, aluHit, al, "10", region, i.refPath, i.fastaHackPath);
      vcfWriter w = vcfWriter(ca, streamy, i.bamPath);
      vcfLine l = w.getVCFLine();

      ASSERT_EQ(l.CHROM, "11");
      ASSERT_EQ(l.POS, 110305761);
      ASSERT_STREQ(l.ID.c_str(), "ME-DeNovo");
      ASSERT_STREQ(l.REF.c_str(), "N");
      ASSERT_STREQ(l.ALT.c_str(), "INS:ME:SVA_B");
      ASSERT_EQ(l.FILTER.HDS and l.FILTER.TDS, true);

      EXPECT_EQ(l.QUAL, 0);
      EXPECT_EQ(l.INFO.NR, 74);
      EXPECT_EQ(l.INFO.NT, 5);
      EXPECT_EQ(l.INFO.NH, 12);
      
      ASSERT_EQ(l.INFO.LT, 60);
      ASSERT_STREQ(l.INFO.cigar.c_str(), "M24D4M102S135");

      ASSERT_EQ(l.INFO.probandGT.DP, 37);
      ASSERT_EQ(l.INFO.probandGT.RO, 19);
      ASSERT_EQ(l.INFO.probandGT.AO, 18);
      ASSERT_EQ(l.INFO.probandGT.genotype, std::make_pair(1,0));

      ASSERT_EQ(l.INFO.parentGTs[0].DP, 34);
      ASSERT_EQ(l.INFO.parentGTs[0].RO, 34);
      ASSERT_EQ(l.INFO.parentGTs[0].AO, 0);
      ASSERT_EQ(l.INFO.parentGTs[0].genotype, std::make_pair(0,0));

      ASSERT_EQ(l.INFO.parentGTs[1].DP, 22);
      ASSERT_EQ(l.INFO.parentGTs[1].RO, 22);
      ASSERT_EQ(l.INFO.parentGTs[1].AO, 0);
      ASSERT_EQ(l.INFO.parentGTs[1].genotype, std::make_pair(0,0));

    }
  }

  ASSERT_EQ(contigFound, true);
  streamy.close();

}

#endif // RUFALU_TESTS_KNOWN_ALUS_HPP
