/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "Ifpack_CrsIlut.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

//==============================================================================
Ifpack_CrsIlut::Ifpack_CrsIlut(const Ifpack_OverlapGraph * OverlapGraph, double DropTol, 
			       double FillTol) 
  : Epetra_Object("Ifpack::CrsIlut"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(OverlapGraph),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), OverlapGraph->OverlapGraph().Comm()),
    DropTol_(DropTol),
    FillTol_(FillTol)
{
}
//==============================================================================
Ifpack_CrsIlut::Ifpack_CrsIlut(const Epetra_RowMatrix * UserMatrix, double DropTol, 
			       double FillTol) 
  : Epetra_Object("Ifpack::CrsIlut"),
    Epetra_CompObject(),
    Ifpack_OverlapFactorObject(UserMatrix),
    Ifpack_OverlapSolveObject(Epetra_Object::Label(), UserMatrix->Comm()),
    DropTol_(DropTol),
    FillTol_(FillTol)
{
}
//==============================================================================
Ifpack_CrsIlut::Ifpack_CrsIlut(const Ifpack_CrsIlut & Source) 
  : Epetra_Object(Source),
    Epetra_CompObject(Source),
    Ifpack_OverlapFactorObject(Source),
    Ifpack_OverlapSolveObject(Source),
    DropTol_(Source.DropTol_),
    FillTol_(Source.FillTol_)
{
}
//==============================================================================
int Ifpack_CrsIlut::ProcessOverlapMatrix(const Epetra_RowMatrix &A)
{

  return(0);
}
//==============================================================================
int Ifpack_CrsIlut::DerivedFactor()
{

  return(0);
}
