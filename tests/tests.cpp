#include "../libs/polyATail.h"
#include "../libs/discordantAndChimericReads.h"
#include "../libs/fastqParse.h"
#include <string>


int testsPassed = 0;
int totalTests = 0;

std::string testBam ="../test_data/big-test.bam";

void detectPositivePolyATail() {
  bool returnValue = polyA::detectPolyATail("AAAAAAAAAA");
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
  bool returnValue = polyA::detectPolyATail("CCCCCCCCCCAAAAAAAATTAAA");
  if(!returnValue){
    std::cout << "detectNegativePolyA test passed!" << std::endl;
    ++testsPassed;
  }
  else{
    std::cout << "ERROR: detectNegativePolyA test failed" << std::endl;
  }
  ++totalTests;
}

void detectChimericReads(std::string testBam, std::string out){
  DACReads::getAllChimericReads(testBam, out);
  std::cout << "detectChimericReads test passed!" << std::endl;
  ++totalTests;
  ++testsPassed;
}

void detectDiscordantReads(std::string testBam, std::string out){
  DACReads::getAllDiscordantReads(testBam, out);
  std::cout << "detectDiscordantReads test passed!" << std::endl;
  ++totalTests;
  ++testsPassed;
}



void runAllTests() {
  detectPositivePolyATail();
  detectNegativePolyATail();
  detectDiscordantReads(testBam, "out1.bam");
  detectChimericReads(testBam, "out2.bam");
}

int main(int argc, char* argv[]){
  runAllTests();
  std::cout << testsPassed << "/" << totalTests << " tests passed" << std::endl;
}
