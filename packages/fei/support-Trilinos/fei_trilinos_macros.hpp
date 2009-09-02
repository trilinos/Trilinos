#ifndef _fei_trilinos_macros_hpp_
#define _fei_trilinos_macros_hpp_

//fei_macros.hpp includes FEI_config.h
#include "fei_macros.hpp"

#ifdef FEI_BYPASS_CONFIG_H

#undef HAVE_FEI_TEUCHOS
#ifndef HAVE_FEI_NO_TEUCHOS
#define HAVE_FEI_TEUCHOS
#endif

#undef HAVE_FEI_EPETRA
#ifndef HAVE_FEI_NO_EPETRA
#define HAVE_FEI_EPETRA
#endif

#undef HAVE_FEI_AZTECOO
#ifndef HAVE_FEI_NO_AZTECOO
#define HAVE_FEI_AZTECOO
#endif

#undef HAVE_FEI_AMESOS
#ifndef HAVE_FEI_NO_AMESOS
#define HAVE_FEI_AMESOS
#endif

#undef HAVE_FEI_ML
#ifndef HAVE_FEI_NO_ML
#define HAVE_FEI_ML
#endif

#undef HAVE_FEI_IFPACK
#ifndef HAVE_FEI_NO_IFPACK
#define HAVE_FEI_IFPACK
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
//FEI_BYPASS_CONFIG_H

#endif
//_fei_trilinos_macros_hpp_

