#ifndef AMESOS_CONFIGDEFS
#define AMESOS_CONFIGDEFS

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

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

