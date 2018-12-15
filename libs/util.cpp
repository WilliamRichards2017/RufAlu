#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unistd.h>


#include "contig.h"
#include "knownAlus.h"
#include "util.h"

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

std::string util::exec(char const* cmd) {
  //std::cout << "executing command " << cmd << std::endl;
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
  return result;
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

const int32_t util::calculateModeKmerDepth(const std::vector<int32_t> & kmerDepth){

  if(kmerDepth.size() == 0){
    return 0;
  }
  int32_t number = kmerDepth.front();
  int32_t mode = number;
  int32_t count = 0;
  int32_t countMode = 0;
  
  for(const auto d : kmerDepth){
    //std::cout << "Number, mode, count, countMode is: " << number << ", " << mode << ", " << count << ", " << countMode << std::endl;
    if(d == number){
      ++count;
    }
    else{
      if (count > countMode){
	countMode = count;
	mode = number;
      }
      count = 1;
      number = d;
    }
  }
  return mode;
}


const int32_t util::countKmerDepth(const std::vector<std::pair<std::string, int32_t> > & kmers){
  std::vector<int32_t> kmerCounts;

  for(const auto & k : kmers){
    if(k.second > 0){
      kmerCounts.push_back(k.second);
    }
  }

  if(kmerCounts.size() == 0){
    return 0;
  }

  return *std::min_element(kmerCounts.begin(), kmerCounts.end());
}

const std::vector<std::string> util::filterKmersFromText(const std::string & textPath, const std::vector<std::string> & kmers){

  std::cout << "Filtering kmers from file: " << textPath << std::endl;
  std::ifstream file(textPath);
  std::string line;

  std::vector<std::string> kmerCounts;
  std::map<std::string, int32_t> kmerMap;

  while(std::getline(file, line)){
    std::istringstream iss(line);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
    if(kmerCount.size() == 2){
      kmerMap.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
    }
    else{
      std::cout << "kmerCount.size() != 2" << std::endl;
    }
  }

  for(auto k : kmers){
    auto it = kmerMap.find(k);
    auto revIt = kmerMap.find(util::revComp(k));
    if(it == kmerMap.end() and revIt == kmerMap.end()){
      std::cout << "Did not find kmer: " << k << " in exclude file" << std::endl;
      kmerCounts.push_back(k);
    }
    else{
      std::cout << "Filtering out reference kmer" << k  << std::endl;
    }
  }
  return kmerCounts;
}


const std::vector<std::pair<std::string, int32_t> > util::countKmersFromText(const std::string & textPath, const std::vector<std::string> & kmers){
  std::ifstream file(textPath);
  std::string line;

  std::vector<std::pair<std::string, int32_t> > kmerCounts;
  std::map<std::string, int32_t> kmerMap;

  while(std::getline(file, line)){
    std::istringstream iss(line);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());
    if(kmerCount.size() == 2){
      kmerMap.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
      //kmerCounts.push_back(std::make_pair(kmerCount[0], atoi(kmerCount[1].c_str())));
      //std::cout << "kmer " << kmerCount[0] << " has count " << kmerCount[1] << std::endl;
    }
  }
  for(auto k : kmers){
    //std::cout << "looking up kmer: " << k << std::endl;
    auto it = kmerMap.find(k);
    auto revIt = kmerMap.find(util::revComp(k));
    if(it != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(it->first, it->second));
      //std::cout << "found kmer in map" << std::endl;
    }
    else if(revIt != kmerMap.end()){
      kmerCounts.push_back(std::make_pair(k, revIt->second));
      //std::cout << "found revComp(kmer) in map" << std::endl;
    }
    else{
      //std::cout << "else statement" << std::endl;
      kmerCounts.push_back(std::make_pair(k, 0));
    }
  }
  //std::cout << "Returning kmerCounts with size of" << kmerCounts.size() << std::endl;
  return kmerCounts;
}

const std::vector<std::pair<std::string, int32_t> > util::countKmersFromJhash(const std::string & jhashPath, const std::vector<std::string> & kmers, const std::string & jellyfishPath){
  //  std::map<std::string, int32_t> ret;

  std::vector<std::pair<std::string, int32_t> > ret;

  for (const auto & kmer : kmers){

    std::string cmd = jellyfishPath + " query " + jhashPath + " " + kmer;
    //std::cout << "executing command: " << cmd << std::endl;
    
    std::string queryOutput = util::exec(cmd.c_str());
    //std::cout << "command output is: " << queryOutput << std::endl;
    std::istringstream iss(queryOutput);
    std::vector<std::string> kmerCount((std::istream_iterator<std::string>(iss)),
				       std::istream_iterator<std::string>());

    if(kmerCount.size() == 2){
      //ret.insert({kmerCount[0], atoi(kmerCount[1].c_str())});
      ret.push_back(std::make_pair(kmerCount[0], atoi(kmerCount[1].c_str())));
    }
  }
  return ret;
}

const std::vector<std::string> util::kmerize(const std::string & sequence, const int32_t & kmerSize){
  int32_t kmercount = 0;
  std::vector<std::string> kmers;

  while(kmercount + kmerSize < sequence.length()){
    std::string kmer = sequence.substr(kmercount, kmerSize);
    kmers.push_back(kmer);
    ++kmercount;
  }
  
  return kmers;
}

const std::vector<BamTools::RefData> util::populateRefData(const std::string & bamPath){
  BamTools::BamReader reader;
  if (!reader.Open(bamPath)){
    std::cout << "Could not open input Bam file" << bamPath << std::endl;
    exit (EXIT_FAILURE);
  }
  return reader.GetReferenceData();
}

const std::string util::pullRefSequenceFromRegion(const std::pair<int32_t, int32_t> & region, const std::string & refPath, const std::vector<BamTools::RefData> & refData, const int32_t & refSize, std::string fastaHackPath){


  std::string cmd = fastaHackPath + " -r " + util::getChromosomeFromRefID(region.first, refData) + ":" + std::to_string(region.second) + ".." + std::to_string(region.second + refSize) + ' ' + refPath;
  std::cout << "Executing command: " << cmd << std::endl;
  std::string out = util::exec(cmd.c_str());

  std::cout << "Returning output: " << out << std::endl;
  return out;
}

const std::string util::revComp (const std::string sequence){
  std::string newString = "";
  //cout << "Start - " << Sequence << "\n";
  for(int i = sequence.size()-1; i>=0; i+= -1) {
    char C = sequence.c_str()[i];
    if (C == 'A')
      {newString += 'T';}
    else if (C == 'C')
      {newString += 'G';}
    else if (C == 'G')
      {newString += 'C';}
    else if (C == 'T')
      {newString += 'A';}
    else if (C == 'N')
      {newString += 'N';}
    else {std::cout << "ERROR IN RevComp - " << C << std::endl;}
  }
  return newString;
}

const std::string util::getChromosomeFromRefID(const int32_t & id, const std::vector<BamTools::RefData> & refData){
  std::string ret = "";
  if(id == -1) {
    ret = "unmapped";
    return ret;
  }
  ret = refData[id].RefName;
  return ret;
}

const bool util::fileExists( const std::string &Filename )
{
  return access( Filename.c_str(), 0 ) == 0;
}


const int32_t util::getLargestClipIndex(const std::vector<int> & clipSizes){
  int32_t index = -1;
  int32_t largestSize = -1;
  for(int32_t i = 0; i < clipSizes.size(); ++i){
    if(clipSizes[i] > largestSize){
      largestSize = clipSizes[i];
      index = i;
    }
  }
  return index;
}

const std::string util::calculateAltSequence(const BamTools::BamAlignment & al){

  int32_t kmerSize = 25;

  std::vector<int> clipSizes;
  std::vector<int> readPositions;
  std::vector<int> genomePositions;

  const std::vector<int32_t> insertionVec = util::getInsertionVec(al);
  al.GetSoftClips(clipSizes, readPositions, genomePositions);

  int32_t index = util::getLargestClipIndex(clipSizes);
  int32_t breakPoint = readPositions[index];


  std::string altSequence = al.QueryBases.substr(breakPoint-kmerSize+1, kmerSize*2-1);

  std::cout << "Alt sequence is: " << altSequence << std::endl;

  return altSequence;
}
