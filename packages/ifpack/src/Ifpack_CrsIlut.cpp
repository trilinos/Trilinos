/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
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

#ifdef HAVE_IFPACK_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#include <ifp_parameters.h>
#endif

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

#ifdef HAVE_IFPACK_TEUCHOS
//==========================================================================
int Ifpack_CrsIlut::SetParameters(const Teuchos::ParameterList& parameterlist,
				  bool cerr_warning_if_unused)
{
  Ifpack::param_struct params;
  params.double_params[Ifpack::fill_tolerance] = FillTol_;
  params.double_params[Ifpack::drop_tolerance] = DropTol_;

  Ifpack::set_parameters(parameterlist, params, cerr_warning_if_unused);

  FillTol_ = params.double_params[Ifpack::fill_tolerance];
  Droptol_ = params.double_params[Ifpack::drop_tolerance];

  return(0);
}
#endif

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
