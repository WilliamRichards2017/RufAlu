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


class util{

 public:


  static void exec(char const*);
  static std::vector<BamTools::BamAlignment> intersect(const char *, const char *);
  static std::string baseName(std::string);
  static bool checkDoubleStranded(std::vector<polyA>);


  static const std::vector<std::string> Split(const std::string&, const char);
  static const char * getRootDirectory(std::string);


 private:

};

#endif // UTIL_H


