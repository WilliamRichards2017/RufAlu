#include <stdexcept>
#include <string>

#include "util.h"
#include "contig.h"
#include "knownAlus.h"

bool util::checkDoubleStranded(std::vector<polyA> t){
  bool f = false;
  bool r = false;
  for(auto it = std::begin(t); it != std::end(t); ++it){
    if(it->isTailReverseStrand()){
      r = true;
    }
    else{
      f = true;
    }
  }
  if(f and r){
    std::cout << "passed check for double stranded tails" << std::endl;
  }
  return (f and r);
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
