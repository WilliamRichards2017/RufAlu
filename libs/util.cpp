#include <stdexcept>
#include <string>


#include "util.h"

void util::printContigWindow(contigWindow c){
  std::cout << "contig is: " << c.contig.Name << std::endl;
  for (auto it = std::begin(c.window); it != std::end(c.window); ++it){
    std::cout << it->QueryBases << std::endl;
  }
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

bool util::overlap(std::pair<int, int> a, std::vector<std::pair<int, int> > b){
  for (auto it = std::begin(b); it != std::end(b); ++it){
    if((a.first >= it->first && a.first <= it->second) || (a.second >= it->first && a.second <= it->second)){
      return true;
    }
  }
  return false;
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

std::vector<BamTools::BamAlignment> util::intersect(const char * a, const char * b){
  
  std::cout << "Inside intersect" << std::endl;

  std::vector<BamTools::BamAlignment> intersection;
  std::vector<BamTools::BamRegion> coords;

  BamTools::BamReader reader;
  reader.LocateIndex();
  if(!reader.Open(a)){
    std::cout << "Could not open the following input bamfile: " << a << std::endl;
    return intersection;
  }

  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    BamTools::BamRegion region = BamTools::BamRegion(al.RefID, al.Position, al.RefID, al.GetEndPosition()); 
    std::cout << "Pushing back coords  " << al.Position << ", " << al.GetEndPosition() << std::endl;
    //if(region.LeftPosition != 0 || region.RightPosition != 0){
    coords.push_back(region);
    //}
  }
  
  reader.Close();
  
  if(!reader.Open(b)){
    std::cout << "could not open the following input bamfile: " << b << std::endl;
  }

  
  for(auto it = std::begin(coords); it != std::end(coords); ++it){
    //std::cout << "looping through regions" << std::endl;
    std::cout << "setting region to be " << coords.back().LeftPosition << ", " << coords.back().LeftRefID << ", " << coords.back().RightPosition << ", " << coords.back().RightRefID << std::endl;
    reader.SetRegion(coords.back());
    //std::cout << "tryna pop" << std::endl;
    coords.pop_back();
    BamTools::BamAlignment bl;
    while(reader.GetNextAlignment(bl)){
      std::cout << "found overlapping read" << std::endl;
      intersection.push_back(bl);
    }
  }
  
  return intersection;  
}

std::vector<BamTools::BamAlignment> util::intersectBams(const char * a, const char * b){
  std::vector<BamTools::BamAlignment> intersection;

  std::vector<std::pair<int, int> > aCoords;
  BamTools::BamReader reader;
  if (!reader.Open(a)){
    std::cout << "Could not open the following input bamfile: " << a << std::endl;
    return intersection;
  }
  BamTools::BamAlignment al;

  while(reader.GetNextAlignment(al)){
    std::pair<int,int> coords = std::make_pair(al.Position, al.GetEndPosition());
    aCoords.push_back(coords);
    //std::cout << "pushing back coords " << coords.first << ", " << coords.second << std::endl;
  }

  if (!reader.Open(b)){
    std::cout << "Could not open the following input bamfile: " << b << std::endl;
    return intersection;
  }
  BamTools::BamAlignment bl;
  
  while(reader.GetNextAlignment(bl)){
    std::pair<int, int> coords = std::make_pair(bl.Position, bl.GetEndPosition());
    if (util::overlap(coords, aCoords)){
      intersection.push_back(bl);
      //std::cout << "found overlap at coords " << coords.first << ", " << coords.second << std::endl;
    }
  }
  return intersection;
}
  
  const char *util::contigsToFastq(std::vector<fastqRead> *contigs, const char * outFile){

  std::ofstream out;
  out.open(outFile);

  for(auto it = std::begin(*contigs); it != std::end(*contigs); ++it){

    //out << (*it)->id << std::endl;                                                                                                                                                                                                         \
                                                                                                                                                                                                                                              
    out << '@' << (*it).name << std::endl;
    out << (*it).seq << std::endl;
    out << '+' << std::endl;
    out << (*it).qual << std::endl;
    //std::cout << "name: " << (*it).name <<  std::endl;
    //std::cout << "sequence: " << (*it).seq <<  std::endl;
    //std::cout << "quality: " << (*it).qual <<  std::endl;

  }

  return outFile;
}
