/*!
 * \file ml_ValidateParameters.h
 *
 * \brief Parameter Validation for ML
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#ifndef ML_VALIDATEPARAMETERS_H
#define ML_VALIDATEPARAMETERS_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Teuchos_ParameterList.hpp"

namespace ML_Epetra
{
  //! Builds a list of "valid" parameters for parameter validation for MultiLevelPreconditioner.
  Teuchos::ParameterList * GetValidMLPParameters();

  //! Builds a list of "valid" smoothing parameters for parameter validation for MultiLevelPreconditioner.
  void SetValidSmooParams(Teuchos::ParameterList *PL, Teuchos::Array<std::string> &smootherList);
  //! Builds a list of "valid" aggregation parameters for parameter validation for MultiLevelPreconditioner.
  void SetValidAggrParams(Teuchos::ParameterList *PL);
  //void GetValidSmootherParameters(ParameterList &PL, Array<std::string> smootherList);

  //! Validates the parameters of inList (warning: level-specific parameters
  //! will not be validated) for MultiLevelPreconditioner.
  bool ValidateMLPParameters(const Teuchos::ParameterList &inList, int depth=5);

  //! Builds a list of "valid" parameters for parameter validation for RefMaxwell.
  Teuchos::ParameterList * GetValidRefMaxwellParameters();

  //! Validates the parameters of inList (warning: level-specific parameters
  //! will not be validated) for RefMaxwell.
  bool ValidateRefMaxwellParameters(const Teuchos::ParameterList &inList);
}
#endif
#endif
