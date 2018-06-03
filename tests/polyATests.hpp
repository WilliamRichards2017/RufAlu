#ifndef RUFALU_TESTS_POLYA_HPP
#define RUFALU_TESTS_POLYA_HPP

#include <string>

#include "polyATail.h"

TEST(PolyATests, p1348){
    BamTools::BamReader reader;
    std::string bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam.generator.Mutations.fastq.bam";
    ASSERT_TRUE(reader.Open(bam));
    ASSERT_TRUE(reader.LocateIndex());
    reader.SetRegion()
    polyA::detectPolyTailClips(BamTools::BamAlignment, uint32_t);
}

TEST(PolyATests, FalsePolyA){
  std::string contig = "AAAAAAACAAAAAAAAAAACAAAAAAAAACAAAAAAAAAACAAAAAAAAAAACAAAAAAAACAAAAAATCAGGAGGCACAAAATATCGGTGCATCTCATTGTTGTGATGCCAAGTTTGATCACTTAAGGCAATATTGACTAGATTTCTCCATCATAA";
  std::pair<bool, int> check = polyA::detectPolyATail(contig);

  ASSERT_FALSE(check.first);
}

TEST(PolyATests, TruePolyT){
  std::string contig = "AGATGCACCGATATTTTGTGCCTCCTGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGTTTTTTTTTTTTTTTTTTTTTTGTTGTTTGTTTTTTTTTTTGTTTGTTTTCTT";
  std::pair<bool, int> check = polyA::detectPolyTTail(contig);
  ASSERT_TRUE(check.first);

}

TEST(PolyATests, FalsePolyT){
  std::string contig = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAGGAGGCACAAAATATCGGTGCATCTCATTGTTGTGATGCCAAGTTTGATCACTTAAGGCAATATTGACTAGATTTCTCCATCATAA";
  std::pair<bool, int> check = polyA::detectPolyTTail(contig);
  ASSERT_FALSE(check.first);
}

#endif //RUFALU_TESTS_POLYA_HPP
