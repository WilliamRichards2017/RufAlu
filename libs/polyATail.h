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

  clipCoords coords_ = {0,0,rtl};

  
  bool isTail();
  bool isTailLeftBound();
  bool isTailReverseStrand();
 private:

  BamTools::BamAlignment al_;
  bool isTail_ = false;
  int32_t tailSize_;
  int32_t longestTail_ = 0;
  
  bool detectPolyTail();
  bool detectTailInWindow(clipCoords, const char);
  void printTailDebug();
  std::vector<clipCoords> getLocalClipCoords();
  void setGlobalClipCoords(int32_t);
  void printClipsAndSeq(BamTools::BamAlignment);

};

#endif // POLYA_H
