#ifndef RUFALU_TESTS_POLYA_HPP
#define RUFALU_TESTS_POLYA_HPP

#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "polyATail.h"
#include "util.h"


TEST(PolyATests, longestTailEdgeCase_p1288){
  BamTools::BamReader reader;
  std::string bam = "/scratch/ucgd/lustre/u0991464/Projects/CEPH.25/1288.bam";
  ASSERT_TRUE(reader.Open(bam));             
  ASSERT_TRUE(reader.LocateIndex());
  reader.SetRegion(3, 157422200, 3, 157422300);
  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    if(al.Name.compare("ST-E00264:200:HT3HFCCXX:4:2111:18477:52063")==0){

      std::cerr << "Read Name : " << al.Name << " is:" << std::endl;
      std::cerr << "Seq is : " << al.QueryBases << " at position: " << al.Position << std::endl;
      std::cerr << "Cigar is: ";
      util::printCigar(al.CigarData);
      polyA tail = {al, 10};
    }
  }
  ASSERT_TRUE(reader.Close());
} 


TEST(PolyATests, trueLeftClipTails_p2788){
  BamTools::BamReader reader;
  std::string bam = "/scratch/ucgd/lustre/u0991464/Projects/CEPH.25/2788.bam";
  ASSERT_TRUE(reader.Open(bam));
  ASSERT_TRUE(reader.LocateIndex());
  reader.SetRegion(10, 9072650, 10, 9072900);
  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //std::cerr << "Name is: " << al.Name << std::endl;                                                                                                                                                                                           
    if(al.Name.compare("E00286:175:HNWHVCCXX:8:1105:2392:29332")==0){
      polyA tail = {al, 10};
      std::cerr << "Longest Tail for read : " << al.Name << " is:" << std::endl;
      std::cerr << tail.getLongestTail();
      std::cerr << "cigar is: ";
      util::printCigar(al.CigarData);
      ASSERT_TRUE( tail.getLongestTail() == 22);
    }
    else if(al.Name.compare("E00324:137:HNVKKCCXX:7:1203:9404:54102")==0){
      polyA tail = {al, 10};
      std::cerr << "Longest Tail for read : " << al.Name << " is:" << std::endl;
      std::cerr << tail.getLongestTail() << std::endl;
      std::cerr << "cigar is: ";
      util::printCigar(al.CigarData);
      ASSERT_TRUE( tail.getLongestTail() == 55);
    }
  }
}


TEST(PolyATests, longestTailLeftClip_p1348){
    BamTools::BamReader reader;
    std::string bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1348.bam";
    ASSERT_TRUE(reader.Open(bam));
    ASSERT_TRUE(reader.LocateIndex());
    reader.SetRegion(6, 110305600, 6, 110305800);
    BamTools::BamAlignment al;
    while(reader.GetNextAlignment(al)){
      //std::cerr << "Name is: " << al.Name << std::endl;
      if(al.Name.compare("ST-E00312:155:HNW77CCXX:1:2222:6999:26589")==0){
	polyA tail = {al, 10};
	
	std::cerr << "Longest Tail for read : " << al.Name << " is:" << std::endl;
	std::cerr << tail.getLongestTail() << std::endl;
	ASSERT_TRUE( tail.getLongestTail() == 59);
      }
    }
}

TEST(PolyATests, longestTailRightClip_p8210){
  BamTools::BamReader reader;
  std::string bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/8210.bam";
  ASSERT_TRUE(reader.Open(bam));
  ASSERT_TRUE(reader.LocateIndex());
  reader.SetRegion(9, 89415300, 9, 89415400);
  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    //std::cerr << "Name is: " << al.Name << std::endl;                                                                                                                          
    if(al.Name.compare("ST-E00267:203:HNYHTCCXX:2:1105:10967:49039")==0){
      polyA tail = {al, 10};
      
      std::cerr << "Longest Tail for read : " << al.Name << " is:" << std::endl;
      std::cerr << tail.getLongestTail();
      ASSERT_TRUE( tail.getLongestTail() == 10);                                                                                                                                 
    }
  }
  ASSERT_TRUE(reader.Close());
}


TEST(PolyATests, longestTailEdgeCase_p1219){
  BamTools::BamReader reader;
  std::string bam = "/scratch/ucgd/lustre/u0691312/analysis/A414_CEPH/alu_samples/1219.bam";
  ASSERT_TRUE(reader.Open(bam));
  ASSERT_TRUE(reader.LocateIndex());
  reader.SetRegion(5, 8538900, 5, 8538999);
  BamTools::BamAlignment al;
  while(reader.GetNextAlignment(al)){
    if(al.Name.compare("ST-E00313:182:HT2VKCCXX:7:1119:6238:66056")==0){
      std::cerr << "Name : " << al.Name << " is:" << std::endl;
      std::cerr << "Seq is : " << al.QueryBases << " at position: " << al.Position << std::endl;
      std::cerr << "Cigar is: ";
      util::printCigar(al.CigarData);
      std::cerr << "QueryBases.size() is: " << al.QueryBases.size() << std::endl;

      polyA tail = {al, 10};		     
      std::cerr << tail.getLongestTail();
            
    }
  }
  ASSERT_TRUE(reader.Close());
}

#endif //RUFALU_TESTS_POLYA_HPP
