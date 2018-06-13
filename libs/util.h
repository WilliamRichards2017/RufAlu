#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <utility> //std::pair

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "intersect.h"
#include "knownAlus.h"
#include "polyATail.h"


class polyA;

struct clipCoords;

struct digitWindow{
  int32_t startPos = -1;
  std::string digits = "";
  int32_t peak = -1;
  int32_t peakPos = -1;
};

class util{

 public:
  static std::vector<clipCoords> getLocalClipCoords(const BamTools::BamAlignment &);
  static void exec(char const*);
  static std::vector<BamTools::BamAlignment> intersect(const char *, const char *);
  static std::string baseName(std::string);
  static bool checkDoubleStranded(std::vector<polyA>);

  static const std::vector<std::string> Split(const std::string &, const char);
  static const char * getRootDirectory(std::string);
  static const std::vector<int32_t> getPeakVector(const BamTools::BamAlignment & al);
  static const std::vector<std::pair<int32_t, int32_t> > getPeaks(const BamTools::BamAlignment &);


 private:
  static const std::vector<int32_t> getPeakVector(BamTools::BamAlignment &);
  static const bool anyOverlap(std::vector<int32_t> const &, std::vector<int32_t> const &);
};

#endif // UTIL_H


