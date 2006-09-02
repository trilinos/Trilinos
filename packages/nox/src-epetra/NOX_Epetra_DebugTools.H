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

#ifndef NOX_EPETRA_DEBUGTOOLS_H
#define NOX_EPETRA_DEBUGTOOLS_H

#include "NOX.H"
#include "NOX_Epetra.H"

//----------------------------------------------------------------------------//
//----------------- DebugTools::ObjectOutput ----------------------------------//
//----------------------------------------------------------------------------//
//
// This class allows linear solver tolerance criteria to be used as a NOX
// convergence test.

namespace NOX {
  namespace Epetra {

class DebugTools 
{

 public:

  enum NOX_EPETRA_OP_TYPE { MATRIX_FREE                 ,
                            FINITE_DIFFERENCE           ,
                            FINITE_DIFFERENCE_COLORNG   ,
                            CRS_MATRIX                  ,
                            NONE                           };

  enum FORMAT_TYPE        { ASCII         = 0           ,
                            MATRIX_MARKET = 1              };

  static NOX_EPETRA_OP_TYPE get_operator_type( const Epetra_Operator *                     );
  static NOX_EPETRA_OP_TYPE get_operator_type( const Teuchos::RefCountPtr<Epetra_Operator> );

  // Factory-like method to create an operator - needed RWH ?
  //static Epetra_Operator * create_nonlinear_operator( NOX_EPETRA_OP_TYPE );

  // Method to compute a matrix from an operator acting on unit Cartesian vect// spanning the space of problem dofs
  static Epetra_CrsMatrix * compute_matrix_using_operator( const Epetra_Operator * );

#ifdef HAVE_NOX_EPETRAEXT
  // I/O methods via EpetraExt 
  static int readVector    ( std::string, const Epetra_Comm &, Epetra_Vector*&                    );
  static void writeVector  ( std::string, const Epetra_Vector &    , FORMAT_TYPE outFormat= ASCII, bool writeMap = true );
  static void writeMatrix  ( std::string, const Epetra_RowMatrix & , FORMAT_TYPE outFormat= ASCII );
  static void writeOperator( std::string, const Epetra_Operator &  , FORMAT_TYPE outFormat= ASCII );
  static void readMatrix   ( std::string,       Epetra_CrsMatrix*&                                );
#endif

}; 

} // Epetra namespace
} // NOX namespace

#endif
