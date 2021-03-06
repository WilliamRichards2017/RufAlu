#include "../src/aluDetect.h"
#include "../src/fastqParse.h"
#include <string>


int testsPassed = 0;
int totalTests = 0;

void detectPositivePolyATail() {
  bool returnValue = aluDetect::detectPolyATail("AAAAAAAAAA");
  if(returnValue){
    std::cout << "detectPositivePolyA test passed!" << std::endl;
    ++testsPassed;
  }
  else{
    std::cout << "ERROR: detectPositivePolyA test failed" << std::endl;
  }
  ++totalTests;
}

void runAllTests() {
  detectPositivePolyATail();
}

int main(int argc, char* argv[]){
  runAllTests();
  std::cout << testsPassed << "/" << totalTests << " tests passed" << std::endl;
}
