/*
!
 * \file ml_ifpack_wrap.h
 *
 * \brief Interface to the Trilinos package Ifpack.
 *
 * The ML/Amesos interface allows ML users to apply Ifpack iterative methods
 * as smoothers.
 *
 *
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */


#ifndef ML_IFPACK_EPETRA_WRAP
#define ML_IFPACK_EPETRA_WRAP

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)

/* Revised section for LevelWrap and RefMaxwell */
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
class Epetra_Vector;

namespace ML_Epetra {
  Epetra_Operator* ML_Gen_Smoother_Ifpack_Epetra(const Epetra_Operator *A,const Epetra_Vector* InvDiagonal,Teuchos::ParameterList & List, std::string printMsg,bool verbose);
}
#endif
#endif
