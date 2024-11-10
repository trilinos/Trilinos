// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TEST_EPETRASUPPORT_HPP
#define PIRO_TEST_EPETRASUPPORT_HPP

#include "Epetra_Comm.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Operator.h"

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

Teuchos::RCP<Epetra_Vector> vectorFromLinOp(const Epetra_Operator &op, int col)
{
  const Teuchos::RCP<Epetra_Vector> result = vectorNew(op.OperatorRangeMap());

  const Teuchos::RCP<Epetra_Vector> rhs = vectorNew(op.OperatorDomainMap());
  const double value = 1.0;
  const int ierr = rhs->ReplaceGlobalValues(1, &value, &col);
  // Some processes might not hold the entry
  TEUCHOS_ASSERT(ierr == 0 || ierr == 1);
  {
    const Epetra_Comm &comm = rhs->Comm();
    int ierrLoc = ierr, ierrSum;
    comm.SumAll(&ierrLoc, &ierrSum, 1);
    // At least one process holds the entry
    TEUCHOS_ASSERT(ierrSum < comm.NumProc());
  }

  op.Apply(*rhs, *result);
  return result;
}

Teuchos::Array<double> arrayFromLinOp(const Epetra_Operator &op, int col)
{
  const Teuchos::RCP<const Epetra_Vector> v = vectorFromLinOp(op, col);
  return arrayFromVector(*v);
}

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_EPETRASUPPORT_HPP */
