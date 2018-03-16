#ifndef DISCORDANTANDCHIMERICREADS_H
#define DISCORDANTANDCHIMERICREADS_H

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include <vector>

class DACReads{

 public:
  DACReads(std::string);
  ~DACReads();

  std::vector<BamTools::BamAlignment*> * getAllAluCandidateReads();
  std::vector<BamTools::BamAlignment*> * getAllChimericReads();
  std::vector<BamTools::BamAlignment*> * getAllUnmappedReads();
 private:
  std::vector<BamTools::BamAlignment*> *chimericReads_;
  std::vector<BamTools::BamAlignment*> *unmappedReads_;
  void setAllChimericReads(std::string);
  void setAllUnmappedReads(std::string);
};

#endif // DISCORDANTANDCHIMERICREADS_H
