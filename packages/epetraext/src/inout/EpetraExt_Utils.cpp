#include "EpetraExt_ConfigDefs.h"
#include "EpetraExt_Utils.h"

// ============================================================================ 
string EpetraExt::toString(const int& x)
{ 
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

// ============================================================================ 
string EpetraExt::toString(const unsigned int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}   
// ============================================================================ 
string EpetraExt::toString(const double& x)
{ 
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

