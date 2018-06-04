#ifndef POLYA_H
#define POLYA_H

#include <string>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

class polyA{
  
 public:
  static std::vector<int32_t> getClipStarts(BamTools::BamAlignment);
  static std::vector<std::pair<int32_t, int32_t> > getClipCoords(BamTools::BamAlignment al);
  static void printClipsAndSeq(BamTools::BamAlignment);
  static bool detectPolyATail(BamTools::BamAlignment);  
  static bool detectPolyTTail(BamTools::BamAlignment);
  static bool detectPolyTailClips(BamTools::BamAlignment, int32_t);
  
 private:

};

#endif // POLYA_H
