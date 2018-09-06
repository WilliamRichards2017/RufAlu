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
  bool isTail();
  bool isTailReverseStrand();
  bool detectPolyTail();


  const int32_t getLongestTail();
  clipCoords getGlobalClipCoords();
  int32_t getGlobalTailStart();

 private:
  
  std::vector<clipCoords> localClipCoords_;
  clipCoords globalClipCoords_;
  bool isTail_ = false;
  int32_t tailSize_;
  int32_t longestTail_ = 0;
  int32_t globalTailStart_;
  int32_t longestTailIndex_ = 0;
  void  setLongestTail();
  void setLocalClipCoords();
  void setGlobalTailStart();
  void setGlobalClipCoords();
  void printClipsAndSeq(BamTools::BamAlignment);
  int32_t detectRightTail(const clipCoords &, const char &);
  int32_t detectLeftTail(const clipCoords &, const char &);

};

#endif // POLYA_H
