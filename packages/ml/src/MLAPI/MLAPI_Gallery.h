#ifndef MLAPI_GALLERY_H
#define MLAPI_GALLERY_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Gallery.h

\brief MLAPI interface to the Galeri package.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_common.h"

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

  Operator Gallery(const std::string ProblemType, const Space& MySpace);

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

#endif // MLAPI_GALLERY_H
