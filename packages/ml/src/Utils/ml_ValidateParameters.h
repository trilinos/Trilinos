
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in this data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

/*!
 * \file ml_ValidateParameters.h
 *
 * \brief Parameter Validation for ML
 *
 */
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

#ifndef ML_VALIDATEPARAMETERS_H
#define ML_VALIDATEPARAMETERS_H

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
  bool ValidateMLPParameters(const Teuchos::ParameterList &inList, int depth=0);

  //! Builds a list of "valid" parameters for parameter validation for RefMaxwell.
  Teuchos::ParameterList * GetValidRefMaxwellParameters();
  
  //! Validates the parameters of inList (warning: level-specific parameters
  //! will not be validated) for RefMaxwell.
  bool ValidateRefMaxwellParameters(const Teuchos::ParameterList &inList);    
}
#endif
#endif
