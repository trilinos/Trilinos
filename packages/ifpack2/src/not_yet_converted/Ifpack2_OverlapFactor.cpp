// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Ifpack2_OverlapFactor.hpp"
#include "Tpetra_Comm.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_VbrMatrix.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"

//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Ifpack2_OverlapGraph * OverlapGraph) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(OverlapGraph),
    UserMatrix_(0)
{
}
//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Tpetra_RowMatrix * UserMatrix) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(0),
    UserMatrix_(UserMatrix)
{
}
//==============================================================================
Ifpack2_OverlapFactor::Ifpack2_OverlapFactor(const Ifpack2_OverlapFactor & Source) 
  : Factored_(Source.Factored_),
    Allocated_(Source.Allocated_),
    ValuesInitialized_(Source.ValuesInitialized_),
    OverlapGraph_(Source.OverlapGraph_),
    UserMatrix_(Source.UserMatrix_)
{
}
//==============================================================================
int Ifpack2_OverlapFactor::InitValues(const Tpetra_RowMatrix * UserMatrix) {
  

  if (OverlapGraph_!=0) {

    Tpetra_CrsMatrix * CrsMatrix = dynamic_cast<Tpetra_CrsMatrix *>(UserMatrix);
    if (CrsMatrix!=0) 
  if (!Allocated()) EPETRA_CHK_ERR(-1); //Must be allocated
  if (ValuesInitialized()) EPETRA_CHK_ERR(1); // Values already init'ed, warn caller
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
//==============================================================================
int Ifpack2_OverlapFactor::Factor() {
  
  if (!ValuesInitialized()) EPETRA_CHK_ERR(-1); // Values must be initialized
  if (Factored()) EPETRA_CHK_ERR(1); // Return with a warning that factor already done
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
