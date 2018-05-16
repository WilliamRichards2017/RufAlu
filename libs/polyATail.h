#ifndef POLYA_H
#define POLYA_H

#include <string>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class polyA{
  
 public:
  static bool detectPolyATail(BamTools::BamAlignment);  
  static bool detectPolyTTail(BamTools::BamAlignment);
  static uint32_t longestTail(BamTools::BamAlignment);
  static bool detectPolyTailClips(BamTools::BamAlignment, uint32_t);
  
 private:
  static std::vector<uint32_t> getClipStarts(BamTools::BamAlignment);
};

#endif // POLYA_H
