#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_Utils.h"

// ============================================================================ 
std::string EpetraExt::toString(const int& x)
{ 
  char s[100];
  std::sprintf(s, "%d", x);
  return std::string(s);
}

// ============================================================================ 
std::string EpetraExt::toString(const unsigned int& x)
{
  char s[100];
  std::sprintf(s, "%d", x);
  return std::string(s);
}   
// ============================================================================ 
std::string EpetraExt::toString(const double& x)
{ 
  char s[100];
  std::sprintf(s, "%g", x);
  return std::string(s);
}

