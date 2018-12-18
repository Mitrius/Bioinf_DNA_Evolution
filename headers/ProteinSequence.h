#include <vector>
#include "base.h"
#include <string>
#include <map>

class ProteinSequence
{
public:
  ProteinSequence();
  ProteinSequence(const std::string &initialSequence);
  std::string baseVec2String();
  base &operator[](int index);
  void pushBack(const base &newBase);
  typename std::vector<base>::iterator begin()
  {
    return inside.begin();
  }
  typename std::vector<base>::iterator end()
  {
    return inside.end();
  }

private:
  std::vector<base> inside;
  std::map<base, char> codings;
  std::map<char, base> deCodings;
};