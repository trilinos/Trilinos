#ifndef AMESOS_CONFIGDEFS
#define AMESOS_CONFIGDEFS

#include "Amesos_config.h"
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_VECTOR
#include <vector>
#else
  Amesos requires STL vector class
#endif

  //  Disable Kundert for now (we test the KundertOO interface, 
  //  but support the Epetra_CrsKundertSparse interface
#undef HAVE_AMESOS_KUNDERT
#endif 

