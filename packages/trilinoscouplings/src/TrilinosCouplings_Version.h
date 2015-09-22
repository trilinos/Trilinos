#ifndef TRILINOSCOUPLINGS_VERSION_H
#define TRILINOSCOUPLINGS_VERSION_H

#include "TrilinosCouplings_ConfigDefs.h"
#include "Trilinos_version.h"

string TrilinosCouplings_Version() { 
  return("TrilinosCouplings in Trilinos " TRILINOS_VERSION_STRING); 
}

#endif /* TRILINOSCOUPLINGS_VERSION_H */
