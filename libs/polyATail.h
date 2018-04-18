#ifndef POLYA_H
#define POLYA_H

#include <string>

class polyA{
  
 public:
  static std::pair<bool, int> detectPolyATail(std::string);  
  static std::pair<bool, int> detectPolyTTail(std::string);
 private:
};

#endif // POLYA_H
