#ifndef __DEBUG_H__
#define __DEBUG_H__



#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include <string>
#include <vector>


class debug{
  
public:
  
  static void printRead(BamTools::BamAlignment al, std::string bamFile);
  static void printReadInRegion(BamTools::BamAlignment al, BamTools::BamRegion region, std::string bamFile);
  
private:
  
};

#endif //__DEBUG_H__deb

