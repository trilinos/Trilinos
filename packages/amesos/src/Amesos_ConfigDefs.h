#ifndef AMESOS_CONFIGDEFS
#define AMESOS_CONFIGDEFS

#include "Amesos_config.h"
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_VECTOR
#include <vector>
#else
  Amesos requires STL vector class
#endif

#endif 

