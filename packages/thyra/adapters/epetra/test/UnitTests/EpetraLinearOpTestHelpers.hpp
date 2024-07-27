// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "EpetraThyraAdaptersTestHelpers.hpp"


namespace {


RCP<Epetra_CrsMatrix> getEpetraMatrix(int numRows, int numCols, double shift=0.0) 
{

  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const Epetra_Map rowMap(numRows, 0, *comm);
  const Epetra_Map domainMap(numCols, numCols, 0, *comm);
 
  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap, numCols));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j) {
    columnIndices[j] = j;
  }

  const int numLocalRows = rowMap.NumMyElements();

  for (int i = 0; i < numLocalRows; ++i) {
    
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10 + shift;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );

  }

  epetraCrsM->FillComplete(domainMap, rowMap);

  return epetraCrsM;
}


} // namespace

#if defined(Thyra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ThyraEpetraAdapters package is deprecated"
#endif
#endif

