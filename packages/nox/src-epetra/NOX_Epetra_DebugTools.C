//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Epetra_DebugTools.H"

#include "Epetra_MapColoring.h" // Needed for FD coloring
#include "Epetra_Time.h" // for performance metrics
#include "Epetra_SerialComm.h" // for performance metrics
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

// Headers needed for FD coloring
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_BlockMapOut.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"

#include "Teuchos_ParameterList.hpp"

//----------------------------------------------------------------------------//
//----------------- static DebugTools methods ---------------------------------//
//----------------------------------------------------------------------------//


NOX::Epetra::DebugTools::NOX_EPETRA_OP_TYPE
NOX::Epetra::DebugTools::get_operator_type( const Epetra_Operator * op )
{

  // Run the gamut of choices to see what we have for op
  // Note that we must order the casts from most to least specialized

  if( NULL == op )
    return NONE;

  const NOX::Epetra::MatrixFree * opMF = dynamic_cast<const NOX::Epetra::MatrixFree *>(op);
  if( opMF )
    return MATRIX_FREE;

  const NOX::Epetra::FiniteDifferenceColoring * opFDC =
    dynamic_cast<const NOX::Epetra::FiniteDifferenceColoring *>(op);
  if( opFDC )
    return FINITE_DIFFERENCE_COLORNG;

  const NOX::Epetra::FiniteDifference * opFD =
    dynamic_cast<const NOX::Epetra::FiniteDifference *>(op);
  if( opFD )
    return FINITE_DIFFERENCE;

  const Epetra_CrsMatrix * crs =
    dynamic_cast<const Epetra_CrsMatrix *>(op);
  if( crs )
    return CRS_MATRIX;

  return NONE;
}

//----------------------------------------------------------------------------

Epetra_CrsMatrix *
NOX::Epetra::DebugTools::compute_matrix_using_operator( const Epetra_Operator * op)
{

  const Epetra_Map & rowMap  =     op->OperatorDomainMap();
  Epetra_CrsMatrix * p_mat   = new Epetra_CrsMatrix(Copy, rowMap, 0);
  Epetra_Vector *    tempVec = new Epetra_Vector(rowMap);
  Epetra_Vector *    tempRes = new Epetra_Vector(rowMap);

  // Show progress on what could be a long operation
  if( 0 == op->Comm().MyPID() )
  {
    std::cout << "****************  CREATING MATRIX FROM OPERATOR ************************ "
	      << std::endl;
    cout << NOX::Utils::fill(72) << endl;
  }
  int totalPerturbations = tempVec->GlobalLength();
  int outFreq = totalPerturbations / 71;
  if( 1 > outFreq )
    outFreq = 1;

  for( int col = 0; col < tempVec->GlobalLength(); ++col )
  {
    tempVec->PutScalar(0.0);
    if( rowMap.MyGID(col) )
      (*tempVec)[col] = 1.0;

    op->Apply(*tempVec, *tempRes);

    // Fill in only the nonzero entries from the apply
    for( int row = 0; row < p_mat->NumMyRows(); ++row)
    {
      if( fabs( (*tempRes)[row] ) > 1.e-12 )
      {
	int ierr = p_mat->InsertGlobalValues( rowMap.GID(row), 1, &(*tempRes)[row], &col );
	if( ierr < 0 )
	{
          std::string msg = //"ERROR (" + ierr + ") : "
	       "NOX::Epetra::DebugTools::compute_matrix_using_operator crsMatrix.ExtractGlobalRowView(...) failed for row : ";// + row;
          throw msg;
	}
      }
    }
    if( (0 == op->Comm().MyPID()) && (0 == (col % outFreq)) )
      cout << "-" << flush;
  }
  if( 0 == op->Comm().MyPID() )
    cout << "*" << endl;

  p_mat->FillComplete();

  delete tempRes;
  delete tempVec;

  return p_mat;
}

#ifdef HAVE_NOX_EPETRAEXT
//----------------------------------------------------------------------------

int
NOX::Epetra::DebugTools::readVector( std::string baseName, const Epetra_Comm & comm, Epetra_Vector*& vec )
{
    int ierr = 0;
    Epetra_Map    * tempMap = NULL;

    std::string mapFilename = baseName + "_map";
    ierr = EpetraExt::MatrixMarketFileToMap( mapFilename.c_str(), comm, tempMap );
    if( ierr != 0 )
      return ierr;
    std::string vecFilename = baseName + "_vec";
    ierr = EpetraExt::MatrixMarketFileToVector( vecFilename.c_str(), *tempMap, vec );

    return ierr;
}

//----------------------------------------------------------------------------

void
NOX::Epetra::DebugTools::writeVector( std::string baseName, const Epetra_Vector & vec, FORMAT_TYPE outFormat, bool writeMap )
{
  if( ASCII == outFormat )
  {
    vec.Print( std::cout );
  }
  else
  {
    std::string fileName = baseName + "_vec";
    EpetraExt::VectorToMatrixMarketFile(fileName.c_str(), vec);
    if( writeMap )
    {
      fileName = baseName + "_map";
      EpetraExt::BlockMapToMatrixMarketFile(fileName.c_str(), vec.Map());
    }
  }
}

//----------------------------------------------------------------------------

void
NOX::Epetra::DebugTools::writeMatrix( std::string baseName, const Epetra_RowMatrix & mat, FORMAT_TYPE outFormat )
{
  if( ASCII == outFormat )
  {
    // We can only Print a CrsMatrix
    const Epetra_CrsMatrix * crsMat = dynamic_cast<const Epetra_CrsMatrix *>(&mat);
    if( crsMat )
      crsMat->Print( std::cout );
    else
    {
      std::string msg = "Could not cast incoming Epetra_RowMatrix to Epetra_CrsMatrix.";
      throw msg;
    }
  }
  else
  {
    std::string fileName = baseName + "_matrix";
    EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), mat);
    fileName = baseName + "_map";
    EpetraExt::BlockMapToMatrixMarketFile(fileName.c_str(), mat.Map());
  }
}

//----------------------------------------------------------------------------

void
NOX::Epetra::DebugTools::readMatrix( std::string baseName, Epetra_CrsMatrix* & p_crsMat )
{

  if( NULL != p_crsMat )
  {
    std::string msg = "Incoming Epetra_CrsMatrix pointer is not NULL.  This could cause a memory leak.";
    throw msg;
  }

  std::string mapFileName = baseName + "_map";

  Epetra_Map * tmpMap = NULL;
  int ierr = EpetraExt::MatrixMarketFileToMap( mapFileName.c_str(), p_crsMat->Comm(), tmpMap );
  if( (0 != ierr) || (NULL == tmpMap) )
  {
    std::string msg = "Could not get Epetra_Map from file."; // \"" + mapFileName + "\"." 
    throw msg;
  }

  std::string matrixFileName = baseName + "_matrix";

  ierr = EpetraExt::MatrixMarketFileToCrsMatrix( matrixFileName.c_str(), *tmpMap, p_crsMat );
  if( (0 != ierr) || (NULL == p_crsMat) )
  {
    std::string msg = "Could not get Epetra_CrsMatrix from file."; // \"" + matrixFileName + "\"." 
    throw msg;
  }

  delete tmpMap;

  return;
}

//----------------------------------------------------------------------------

void
NOX::Epetra::DebugTools::writeOperator( std::string baseName, const Epetra_Operator & op, FORMAT_TYPE outFormat )
{
  NOX_EPETRA_OP_TYPE opTyp = get_operator_type( &op );

  if( NONE == opTyp )
    return;

  // Now get (or compute? - could be a user-option ) the matrix


  switch( opTyp )
  {
    case MATRIX_FREE	  :
      Epetra_CrsMatrix * tmpMatrix;
      tmpMatrix = compute_matrix_using_operator( &op );
      writeMatrix( baseName, *tmpMatrix );
      delete tmpMatrix;
      break;

    case FINITE_DIFFERENCE	   :
    case FINITE_DIFFERENCE_COLORNG :
      writeMatrix( baseName, dynamic_cast<const NOX::Epetra::FiniteDifference &>(op).getUnderlyingMatrix(), outFormat );
      break;

    case CRS_MATRIX	  :
      writeMatrix( baseName, dynamic_cast<const Epetra_CrsMatrix &>(op), outFormat );
      break;

    default			  :
      std::string msg = "Could not get a valid Matrix from incoming Epetra_Operator."; // of type " + opTyp + ".";
      throw msg;
  }

  return;
}
#endif

