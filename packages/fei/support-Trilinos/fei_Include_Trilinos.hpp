#ifndef _fei_Include_Trilinos_hpp_
#define _fei_Include_Trilinos_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_trilinos_macros.hpp"

#ifdef HAVE_FEI_EPETRA

#ifndef FEI_SER
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_LinearProblem.h>
#endif

#ifdef HAVE_FEI_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#endif

//some undefs to avoid warnings about "Attempt to redefine..." in
//AztecOO_config.h, which is included downstream of AztecOO.h
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef F77_FUNC
#undef F77_FUNC_

#ifdef HAVE_FEI_AZTECOO
#include <AztecOO.h>
#endif

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#undef F77_FUNC
#undef F77_FUNC_

#ifdef HAVE_FEI_AMESOS

#include <Amesos_config.h>
#include <Amesos.h>

#endif

#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#undef F77_FUNC
#undef F77_FUNC_

#ifdef HAVE_FEI_IFPACK

#include <Ifpack.h>

#endif

#endif

