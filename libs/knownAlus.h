#include "fastqParse.h"
#include <vector>
#include <tuple>

using FastaRecord = std::tuple<std::string, std::string>;

class KnownAlus{
 public:
  KnownAlus(std::string, std::string);
  ~KnownAlus();
  
  std::vector<std::string*> * getContigsContainingKnownAlus();
  

 private:
  std::string contigFilePath_;
  std::string aluFilePath_;
  std::vector<std::string*> *contigsContainingKnownAlus_;
  
  void findContigsContainingKnownAlus();

};
