#include "../libs/polyATail.h"
#include "../libs/discordantAndChimericReads.h"
#include "../libs/fastqParse.h"
#include <string>


int testsPassed = 0;
int totalTests = 0;

std::string testBam ="/uufs/chpc.utah.edu/common/home/u0401321/RufAlu/test_data/big-test.bam";
//std::string testBam = "/uufs/chpc.utah.edu/common/home/u0991464/RUFUS.test.set/Family1.mother.bam";

void detectPositivePolyATail() {
  bool returnValue = polyA::detectPolyATail("ggccgggcgcggtggctcacgcctgtaatcccagcactttgggaggccgaggcgggaggattgcttgagcccaggagttcgagaccagcctgggcaacatagcgagaccccgtctctacaaaaaatacaaaaattagccgggcgtggtggcgcgcgcctgtagtcccagctactcgggaggctgaggcaggaggatcgcttgagcccaggagttcgaggctgcagtgagctatgatcgcgccactgcactccagcctgggcgacagagcgagaccctgtctcaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
  if(returnValue){
    std::cout << "detectPositivePolyA test passed!" << std::endl;
    ++testsPassed;
  }
  else{
    std::cout << "ERROR: detectPositivePolyA test failed" << std::endl;
  }
  ++totalTests;
}

void detectNegativePolyATail() {
  bool returnValue = polyA::detectPolyATail("cccccccaaaaaaaagccgggggttttttttt");
  if(!returnValue){
    std::cout << "detectNegativePolyA test passed!" << std::endl;
    ++testsPassed;
  }
  else{
    std::cout << "ERROR: detectNegativePolyA test failed" << std::endl;
  }
  ++totalTests;
}

void buildDACReads(std::string filePath){
  //DACReads dacReads(filePath);
  DACReads *dacReads = new DACReads(filePath);
  ++totalTests;
  ++testsPassed;
  delete dacReads;
}

void runAllTests() {
  detectPositivePolyATail();
  detectNegativePolyATail();
  buildDACReads(testBam);
}

int main(int argc, char* argv[]){
  runAllTests();
  std::cout << testsPassed << "/" << totalTests << " tests passed" << std::endl;
}
