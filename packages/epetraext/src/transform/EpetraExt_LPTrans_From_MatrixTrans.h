//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
                                                                                                    
#ifndef EpetraExt_LINEARPROBLEM_MATRIXTRANS_H
#define EpetraExt_LINEARPROBLEM_MATRIXTRANS_H

#include <EpetraExt_Transform.h>

class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_Export;
class Epetra_Import;

namespace EpetraExt {

//! Adaptation of an Epetra_CrsMatrix Transform to a Epetra_LinearProblem Transform
class LinearProblem_MatrixTrans : public SameTypeTransform<Epetra_LinearProblem>
{
  SameTypeTransform<Epetra_CrsMatrix> & matrixTrans_;

  Epetra_Import * Importer_;
  Epetra_Export * Exporter_;

  Epetra_LinearProblem * OldProblem_;
  Epetra_CrsMatrix * OldMatrix_;
  Epetra_MultiVector * OldLHS_;
  Epetra_MultiVector * OldRHS_;
  Epetra_Map * OldRowMap_;

  Epetra_LinearProblem * NewProblem_;
  Epetra_CrsMatrix * NewMatrix_;
  Epetra_MultiVector * NewLHS_;
  Epetra_MultiVector * NewRHS_;

 public:

  //! Destructor
  ~LinearProblem_MatrixTrans();

  //! Constructor
  LinearProblem_MatrixTrans( SameTypeTransform<Epetra_CrsMatrix> & matrix_trans )
  : matrixTrans_(matrix_trans),
    OldProblem_(0),
    OldMatrix_(0),
    OldLHS_(0),
    OldRHS_(0),
    OldRowMap_(0),
    NewProblem_(0),
    NewMatrix_(0),
    NewLHS_(0),
    NewRHS_(0),
    Importer_(0),
    Exporter_(0)
  {}

  //! Transform Operator
  /*! Constructs a Epetra_LinearProblem
   */
  NewTypeRef operator()( OriginalTypeRef orig );

  //! Forward Data Migration
  bool fwd();

  //! Reverse Data Migration
  bool rvs();

};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_MATRIXTRANS_H

