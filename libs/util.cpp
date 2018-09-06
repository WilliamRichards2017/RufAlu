#include <stdexcept>
#include <string>
#include <algorithm>

#include "util.h"
#include "contig.h"
#include "knownAlus.h"

const bool util::isReadLeftBound(const std::vector<BamTools::CigarOp> & cigOps){
  if(cigOps[0].Type == 'S'){
    return true;
  }
  return false;
}

const std::vector<std::string> util::getClipSeqs(const BamTools::BamAlignment & al){
  std::vector<std::string> clipSeqs;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;
  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);

  al.GetSoftClips(clipSizes, readPositions, genomePositions);
  for(int i = 0; i < readPositions.size(); ++i){
    //std::cout << "Clipped seq for read is: " << al.QueryBases.substr(readPositions[i], clipSizes[i]) << std::endl;
    clipSeqs.push_back(al.QueryBases.substr(readPositions[i]+insertionVec[i], clipSizes[i]));
  }
  return clipSeqs;
}

const std::pair<std::string, int32_t> util::getHighestQualityAluHit(const std::vector<std::pair<std::string, int32_t> > & aluHits){
  std::pair<std::string, int32_t> maxQualHit = std::make_pair("", -1);
  for(auto a : aluHits){
    if(a.second > maxQualHit.second){
      maxQualHit = a;
    }
  }
  //std::cout << "highest quality hit for " << maxQualHit.first << " is " << maxQualHit.second << std::endl;
  return maxQualHit;
}

const int32_t util::getLongestTail(const std::vector<polyA> & tails){
  int32_t maxTail = 0;
  for(auto t : tails){
    if(t.getLongestTail() > maxTail){
      maxTail = t.getLongestTail();
    }
  }
  return maxTail;
}

const clipCoords util::isWithinRegion(clipCoords & cc , const std::pair<int32_t, int32_t> & p) {
  //std::cout << "Checking if " << cc.clipStart << " is withing region: " << p.first << "," << p.second << std::endl; 
  if(cc.clipStart >= p.first - 2 and cc.clipStart <= p.second + 2){
    return cc;
  }
  clipCoords nullC = {-1, -1, ltr, 0};
  return nullC;
}

const clipCoords util::intersectPeaksAndClips(const std::vector<std::pair<int32_t, int32_t> > & peakVec, const std::vector<clipCoords> & clipVec){
  for(auto c : clipVec){
    for(auto p : peakVec){
      if(util::isWithinRegion(c, p).clipStart != -1){
	return util::isWithinRegion(c, p);
      }
    }
  }
  clipCoords nullC = {-1, -1, ltr, 0};
  return nullC;
}

const std::vector<int32_t> util::getInsertionVec(const BamTools::BamAlignment & al){
  const std::vector<BamTools::CigarOp> cig = al.CigarData;
  std::vector<int32_t> insertionVec;
  int32_t indel = 0;
  for(auto c : cig){
    if(c.Type =='S'){
      insertionVec.push_back(indel);
      indel = 0;
    }
    else if(c.Type == 'I'){
      //indel += c.Length;
    }
    else if(c.Type == 'D'){
      indel -= c.Length;
    }
  }
  return insertionVec;
}

const std::vector<clipCoords> util::getLocalClipCoords(const BamTools::BamAlignment & al) {
  std::vector<clipCoords> coordsVec = {};
  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;



  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  clipCoords c = {};

  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);

  for(int32_t i = 0; i < readPositions.size(); ++i){

    if(readPositions[i]-clipSizes[i]==0){
      c.clipDir = rtl;
      c.clipStart = readPositions[i] + insertionVec[i];
      c.clipEnd = 0;
    }
    else{
      c.clipDir = ltr;
      c.clipStart = readPositions[i] + insertionVec[i];
      c.clipEnd = c.clipStart + clipSizes[i];
    }
    c.index = i;
    coordsVec.push_back(c);

    //std::cout << "detecting tail for seq: " << al_.QueryBases << std::endl;                                                                                                        
    //std::cout << "with clip coords " << c.clipStart << ", " << c.clipEnd << ", " << c.clipDir << std::endl;                                                                        
  }
  return coordsVec;
}

std::pair<int32_t, int32_t> getWindowPeak(std::string digitWindow){
  int peak = 0;
  for(auto c: digitWindow){
    if(atoi(&c) > peak){
      peak = atoi(&c);
    }
    
  }
  return std::make_pair(0,0);
}


const std::vector<int32_t> util::getPeakVector(const BamTools::BamAlignment & al){
  
  std::vector<int32_t> peakVec;
  //std::cout << "peakVector : ";
  for(auto c : al.Qualities){
    peakVec.push_back(int(c)-33); // ascii conversion
    //std::cout << peakVec.back() << ", ";
  }
  //std::cout << std::endl;
  return peakVec;
  
}

const std::vector<std::pair<int32_t, int32_t> > util::getPeaks(const BamTools::BamAlignment & al) {

  
  auto amp = util::getPeakVector(al);
  int wideStart = -1;                 // The start of any current wide peak
  
  int grad = -1;                      // Sign of gradient (almost)
  //    =  1 for increasing
  //    =  0 for level AND PREVIOUSLY INCREASING (so potential wide peak)
  //    = -1 for decreasing OR level, but previously decreasing
  // A sharp peak is identified by grad=1 -> grad=-1
  // A wide  peak is identified by grad=0 -> grad=-1

  std::vector<std::pair<int32_t, int32_t> > peakCoords;

  for (int i = 0; i < amp.size() - 1; i++) {
    if(amp[i+1] < amp[i]){    
      if(grad == 1){
	//std::cout << "Sharp peak of " << amp[i] << " at i = " << i << '\n';
	peakCoords.push_back(std::make_pair(i,i));
      }
      else if(grad == 0){
	peakCoords.push_back(std::make_pair(wideStart,i));
	//std::cout << "Wide peak of " << amp[i] << " from i = " << wideStart << " to " << i << '\n';
      }
      grad = -1;
    }
    else if(amp[i+1] == amp[i]){   // Check for start of a wide peak
      if(grad == 1){
	wideStart = i;
	grad = 0;
      }
    }
    else{
      grad = 1;
    }
  }
  return peakCoords;
}


const bool util::anyOverlap(std::vector<int32_t> const & a, std::vector<int32_t> const & b){
  return std::find_first_of (a.begin(), a.end(),
			     b.begin(), b.end()) != a.end();
}


//TODO: Refactor to be optimized with new position information rather than left/right
bool util::checkDoubleStranded(std::vector<polyA> t){
  std::vector<int32_t> reverseTailStarts = {};
  std::vector<int32_t> forwardTailStarts = {};
  for(auto it = std::begin(t); it != std::end(t); ++it){
    if(it->isTailReverseStrand()){
      reverseTailStarts.push_back(it->getGlobalClipCoords().clipStart);
    }
    else{
      forwardTailStarts.push_back(it->getGlobalClipCoords().clipStart);
    }
  }
  
  return util::anyOverlap(reverseTailStarts, forwardTailStarts);
}

std::string util::baseName(std::string path){
  return path.substr(path.find_last_of("/\\")+1);
}

void util::exec(char const* cmd) {
  std::cout << "executing command " << cmd << std::endl;
  char buffer[512];
  std::string result = "";
  FILE* pipe = popen(cmd, "r");
  if (!pipe) throw std::runtime_error("popen() failed!");
  try {
    while (!feof(pipe)) {
      if (fgets(buffer, 512, pipe) != NULL)
        result += buffer;
    }
  } catch (...) {
    pclose(pipe);
    throw;
  }
  pclose(pipe);
  return;
}

const std::vector<std::string> util::Split(const std::string& line, const char delim)
{
  std::vector<std::string> tokens;
  std::stringstream lineStream(line);
  std::string token;
  while(getline(lineStream, token, delim)){
    tokens.push_back(token);
  }
  return tokens;
}

const char * util::getRootDirectory(std::string rufAluPath){
  
  std::vector<std::string> tokens = util::Split(rufAluPath, '/');

  std::string rootDir = "";

  while((tokens.front() != "RufAlu") && (!tokens.empty())){
    rootDir+="/";
    rootDir+= tokens.front();
    tokens.erase(tokens.begin());
  }

  if(tokens.front() !="RufAlu"){
    std::cout << "Unable to parse path to RufAlu, exiting now" << std::endl;
    exit (EXIT_FAILURE);
  }
  return rootDir.c_str();

}

void util::printCigar(const std::vector<BamTools::CigarOp> & cig){
  std::cout << std::endl;
  for(auto c : cig){
    std::cout << c.Type << c.Length;
  }
  std::cout << std::endl;
}
