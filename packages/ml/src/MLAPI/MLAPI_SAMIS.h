#ifndef MLAPI_SAMIS_H
#define MLAPI_SAMIS_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
#include <iostream>

namespace MLAPI {

class Operator;

/*!
  \file MLAPI_SAMIS

  \brief Input matrix and kernel vectors from SAMIS format.

  \author Marzio Sala, SNL 9214 and Marian Brezina, UC Boulder.

  \date Last updated on Mar-05.
*/

  // ====================================================================== 
  //! Reads symmetric matrix from SAMIS binary format.
  // ====================================================================== 

  void ReadSAMISMatrix(const char *filen, Operator& A, int& NumPDEEqns);


  // ====================================================================== 
  //! Reads null space vectors from SAMIS binary format.
  // ====================================================================== 

  void ReadSAMISKernel(const char *myKerFileName, MultiVector& A);

} // namespace MLAPI

#endif // HAVE_ML_MLAPI

#endif  // MLAPI_SAMIS_H
