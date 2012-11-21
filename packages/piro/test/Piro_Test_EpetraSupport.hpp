/*
// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER
*/

#ifndef PIRO_TEST_EPETRASUPPORT_HPP
#define PIRO_TEST_EPETRASUPPORT_HPP

#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

namespace Piro {

namespace Test {

Teuchos::Array<double> arrayFromVector(const Epetra_Vector &vec)
{
  const Epetra_BlockMap &vecMap = vec.Map();
  const Epetra_LocalMap localMap(vecMap.NumGlobalElements(), 0, vec.Comm());
  const Epetra_Import localImporter(localMap, vecMap);

  Epetra_Vector localVec(localMap, false);
  localVec.Import(vec, localImporter, Insert);

  return Teuchos::Array<double>(localVec.Values(), localVec.Values() + localVec.MyLength());
}

Teuchos::Array<double> arrayFromVector(const Epetra_MultiVector &mv, int col)
{
  return arrayFromVector(*mv(col));
}

Teuchos::RCP<Epetra_Vector> vectorNew(const Epetra_BlockMap &map)
{
  return Teuchos::rcp(new Epetra_Vector(map));
}

Teuchos::RCP<Epetra_MultiVector> multiVectorNew(const Epetra_BlockMap &map, int vectorCount)
{
  return Teuchos::rcp(new Epetra_MultiVector(map, vectorCount));
}

Teuchos::RCP<Epetra_MultiVector> multiVectorNew(const Epetra_BlockMap &map, const Epetra_BlockMap &colMap)
{
  TEUCHOS_ASSERT(colMap.NumGlobalElements() == colMap.NumMyElements());
  return multiVectorNew(map, colMap.NumGlobalElements());
}

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_EPETRASUPPORT_HPP */
