// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef MY_SDM_HELPERS_HPP
#define MY_SDM_HELPERS_HPP

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseHelpers.hpp"

namespace Anasazi{

//! Helper functions for SerialDenseMatrices used in Anasazi regression tests.
template <class ScalarType>
void randomSDM( Teuchos::SerialDenseMatrix<int, ScalarType>& matrix )
{
  Teuchos::randomSyncedMatrix( matrix );
} 

}

#endif // MY_SDM_HELPERS_HPP
