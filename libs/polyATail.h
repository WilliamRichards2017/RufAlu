#ifndef POLYA_H
#define POLYA_H

#include <string>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

typedef enum {ltr, rtl} direction;

struct clipCoords{
  int32_t clipStart;
  int32_t clipEnd;
  direction clipDir;
  int32_t index;
};

class polyA{
  
 public:
  polyA(BamTools::BamAlignment, int32_t);
  ~polyA();

  BamTools::BamAlignment al_;
  clipCoords coords_ = {0,0,rtl};
  std::vector<clipCoords> getLocalClipCoords();
  

  
  bool isTail();
  bool isTailLeftBound();
  bool isTailReverseStrand();
  const int32_t getLongestTail();
  int32_t getGlobalTailStart();

 private:

  bool isTail_ = false;
  int32_t tailSize_;
  int32_t longestTail_ = 0;
  int32_t globalTailStart_;
  int32_t detectRightTail(const clipCoords &, const char &);
  int32_t detectLeftTail(const clipCoords &, const char &);
  void setLongestTail(const std::vector<clipCoords> & c);
  bool detectPolyTail();
  bool detectTailInWindow(const clipCoords &, const char &);
  void printTailDebug();
  void setGlobalClipCoords(int32_t);
  void printClipsAndSeq(BamTools::BamAlignment);

};

#endif // POLYA_H
