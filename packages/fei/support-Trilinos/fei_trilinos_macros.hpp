#ifndef _fei_trilinos_macros_hpp_
#define _fei_trilinos_macros_hpp_

//fei_macros.hpp includes FEI_config.h
#include "fei_macros.hpp"

#undef HAVE_TEUCHOS
#ifndef HAVE_NO_TEUCHOS
#define HAVE_TEUCHOS
#endif

#undef HAVE_EPETRA
#ifndef HAVE_NO_EPETRA
#define HAVE_EPETRA
#endif

#undef HAVE_AZTECOO
#ifndef HAVE_NO_AZTECOO
#define HAVE_AZTECOO
#endif

#undef HAVE_AMESOS
#ifndef HAVE_NO_AMESOS
#define HAVE_AMESOS
#endif

#undef HAVE_ML
#ifndef HAVE_NO_ML
#define HAVE_ML
#endif

#undef HAVE_IFPACK
#ifndef HAVE_NO_IFPACK
#define HAVE_IFPACK
#endif

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef F77_FUNC
#undef F77_FUNC_

#endif

