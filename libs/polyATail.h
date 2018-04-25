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
 private:
};

#endif // POLYA_H
