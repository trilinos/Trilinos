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

#include "Ifpack_Jacobi.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Ifpack_OverlapGraph * OverlapGraph, bool UseReciprocal, 
			       int NumSteps) 
  : Epetra_Object("Ifpack::Jacobi"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(OverlapGraph),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), OverlapGraph->OverlapGraph().Comm()),
    UseReciprocal_(UseReciprocal),
    NumSteps_(NumSteps),
    DiagValues_(0)
{
}
//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Epetra_RowMatrix * UserMatrix, bool UseReciprocal, 
			       int NumSteps) 
  : Epetra_Object("Ifpack::Jacobi"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(UserMatrix),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), UserMatrix->Comm()),
    UseReciprocal_(UseReciprocal),
    NumSteps_(NumSteps),
    DiagValues_(0)
{
}
//==============================================================================
Ifpack_Jacobi::Ifpack_Jacobi(const Ifpack_Jacobi & Source) 
  : Epetra_Object(Source),
    Epetra_CompObject(Source),
    Ifpack_OverlapFactorObject(Source),
    Ifpack_OverlapSolveObject(Source),
    UseReciprocal_(Source.UseReciprocal_),
    NumSteps_(Source.NumSteps_),
    DiagValues_(Source.DiagValues_)
{
  if (DiagValues_!=0) DiagValues_ = new Epetra_Vector(*DiagValues_);
}
//==============================================================================
int Ifpack_Jacobi::ProcessOverlapMatrix(const Epetra_RowMatrix &A)
{

  if (DiagValues_=0) DiagValues_ = new Epetra_Vector(A.RowMap()); // Allocate if necessary
  EPETRA_CHK_ERR(A.ExtraMyDiagonalCopy(*DiagValues_)); // Get Diagonal Values of A

  // Compute inverse of diagonal if specified
  if (UseReciprocal()) {EPETRA_CHK_ERR(DiagValues().Reciprocal(DiagValues()));}

  if (NumSteps!=1) EPETRA_CHK_ERR(-1); // NumSteps != 1 not supported yet.
  return(0);
}
//==============================================================================
int Ifpack_Jacobi::DerivedFactor()
{

  return(0);
}
