/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
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

#include "Ifpack_OverlapFactor.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Ifpack_OverlapGraph * OverlapGraph) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(OverlapGraph),
    UserMatrix_(0)
{
}
//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Epetra_RowMatrix * UserMatrix) 
  : Factored_(false),
    Allocated_(false),
    ValuesInitialized_(false),
    OverlapGraph_(0),
    UserMatrix_(UserMatrix)
{
}
//==============================================================================
Ifpack_OverlapFactor::Ifpack_OverlapFactor(const Ifpack_OverlapFactor & Source) 
  : Factored_(Source.Factored_),
    Allocated_(Source.Allocated_),
    ValuesInitialized_(Source.ValuesInitialized_),
    OverlapGraph_(Source.OverlapGraph_),
    UserMatrix_(Source.UserMatrix_)
{
}
//==============================================================================
int Ifpack_OverlapFactor::InitValues(const Epetra_RowMatrix * UserMatrix) {
  

  if (OverlapGraph_!=0) {

    Epetra_CrsMatrix * CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(UserMatrix);
    if (CrsMatrix!=0) 
  if (!Allocated()) EPETRA_CHK_ERR(-1); //Must be allocated
  if (ValuesInitialized()) EPETRA_CHK_ERR(1); // Values already init'ed, warn caller
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
//==============================================================================
int Ifpack_OverlapFactor::Factor() {
  
  if (!ValuesInitialized()) EPETRA_CHK_ERR(-1); // Values must be initialized
  if (Factored()) EPETRA_CHK_ERR(1); // Return with a warning that factor already done
  
  EPETRA_CHK_ERR(DerivedFactor()); // Call Derived class factorization
  SetValuesInitialized(false);
  SetFactored(true);
  return(0);
}
