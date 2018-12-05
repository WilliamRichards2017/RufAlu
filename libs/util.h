#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <string>
#include <utility> //std::pair

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "knownAlus.h"
#include "polyATail.h"


class polyA;

struct clipCoords;

class util{

public:
  static std::string exec(char const*);
  static std::vector<BamTools::BamAlignment> intersect(const char *, const char *);
  static std::string baseName(std::string);
  static bool checkDoubleStranded(std::vector<polyA>);
  static const std::vector<std::string> getClipSeqs(const BamTools::BamAlignment &);
  static const int32_t getLongestTail(const std::vector<polyA> &);
  static const std::vector<std::string> Split(const std::string &, const char);
  static const char * getRootDirectory(std::string);
  static const std::vector<std::pair<int32_t, int32_t> > getPeaks(const BamTools::BamAlignment &);
  static const clipCoords intersectPeaksAndClips(const std::vector<std::pair<int32_t, int32_t> > &, const std::vector<clipCoords> &);
  static const std::vector<int32_t> getInsertionVec(const BamTools::BamAlignment &);
  static void printCigar(const std::vector<BamTools::CigarOp> &);
  static const std::vector<int32_t> getPeakVector(const BamTools::BamAlignment &);
  static const std::pair<std::string, int32_t> getHighestQualityAluHit(const std::vector<std::pair<std::string, int32_t > > &);
  static const bool isReadLeftBound(const std::vector<BamTools::CigarOp> &);
  static const std::vector<clipCoords> getLocalClipCoords(const BamTools::BamAlignment &);
  static const std::vector<std::pair<std::string, int32_t> > countKmersFromText(const std::string &, const std::vector<std::string> &);
  static const int32_t countKmerDepth(const std::vector<std::pair<std::string, int32_t> > &);
  static const int32_t calculateModeKmerDepth(const std::vector<int32_t> &);
  static const std::vector<std::string> kmerize(const std::string &, const int32_t &);
  static const std::vector<BamTools::RefData> populateRefData(const std::string &);
  static const std::string pullRefSequenceFromRegion(const std::pair<int32_t, int32_t> &, const std::string &, const std::vector<BamTools::RefData> &, const int32_t &);
  static const std::string revComp(const std::string);
  static const std::string getChromosomeFromRefID(const int32_t &, const std::vector<BamTools::RefData> &);


 private:
  static const clipCoords isWithinRegion(clipCoords &, const std::pair<int32_t, int32_t> &);
  static const bool anyOverlap(std::vector<int32_t> const &, std::vector<int32_t> const &);
};

#endif // UTIL_H


