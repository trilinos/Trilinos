// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
#ifndef MY_SDM_HELPERS_HPP
#define MY_SDM_HELPERS_HPP

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi{

//! Helper functions for SerialDenseMatrices used in Anasazi regression tests.
template <class ScalarType>
void randomSDM( Teuchos::SerialDenseMatrix<int, ScalarType>& matrix )
{
  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();

  const int procRank = rank(*comm);

  // Construct a separate serial dense matrix and synchronize it to get around
  // input matrices that are subviews of a larger matrix.
  Teuchos::SerialDenseMatrix<int, ScalarType> newMatrix( matrix.numRows(), matrix.numCols() );
  if (procRank == 0)
    newMatrix.random();
  else
    newMatrix.putScalar( Teuchos::ScalarTraits<ScalarType>::zero() );

  broadcast(*comm, 0, matrix.numRows()*matrix.numCols(), newMatrix.values());

  // Assign the synchronized matrix to the input.
  matrix.assign( newMatrix );
} 

}

#endif // MY_SDM_HELPERS_HPP
