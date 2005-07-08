#ifndef MLAPI_GALLERY_H
#define MLAPI_GALLERY_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

/*!
\file MLAPI_Gallery

\brief Matrix creation functions.

\author Marzio Sala, SNL 9214

\date Last updated on Mar-05.
*/

  class Space;
  class Operator;

  // ====================================================================== 
  //! Creates a matrix using the TRIUTILS gallery.
  // ====================================================================== 
  
  Operator Gallery(const string ProblemType, const Space& MySpace);

  // ====================================================================== 
  //! Creates a 1D shifted Laplacian.
  // ====================================================================== 
  
  Operator GetShiftedLaplacian1D(const int NX, const double Factor = 0.99);

  // ====================================================================== 
  //! Creates a 2D shifted Laplacian.
  // ====================================================================== 
  
  Operator GetShiftedLaplacian2D(const int NX, const int NY, 
                                 const double Factor = 0.99,
                                 const bool RandomScale = false);
  
  // ====================================================================== 
  //! Reads a matrix in MATLAB format.
  // ====================================================================== 
  
  Operator ReadMatrix(const char* FileName);
  
  // ====================================================================== 
  //! Creates a recirculation problem in 2D.
  // ====================================================================== 
  
  Operator GetRecirc2D(const int NX, const int NY, const double conv,
                       const double diff);

  // ====================================================================== 
  //! Populates a list from specified file.
  // ====================================================================== 
  
  Teuchos::ParameterList ReadParameterList(const char* FileName);

}

#endif // HAVE_ML_MLAPI

#endif // MLAPI_GALLERY_H
