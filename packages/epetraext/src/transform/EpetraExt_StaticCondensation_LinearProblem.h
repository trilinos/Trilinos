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

#ifndef EpetraExt_LINEARPROBLEM_STATICCONDENSATION_H
#define EpetraExt_LINEARPROBLEM_STATICCONDENSATION_H

#include <EpetraExt_Transform.h>

class Epetra_LinearProblem;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_Export;

namespace EpetraExt {

class LinearProblem_StaticCondensation : public SameTypeTransform<Epetra_LinearProblem>
{
  const int degree_;
  const bool verbose_;

  const Epetra_Map * OldRowMap_;
  Epetra_LinearProblem * OldProblem_;
  Epetra_MultiVector * OldRHS_;
  Epetra_MultiVector * OldLHS_;
  const Epetra_CrsGraph * OldGraph_;
  Epetra_CrsMatrix * OldMatrix_;

  Epetra_Export * Exporter_;

  Epetra_Map * NewRowMap_;
  Epetra_Map * NewColMap_;
  Epetra_LinearProblem * NewProblem_;
  Epetra_MultiVector * NewRHS_;
  Epetra_MultiVector * NewLHS_;
  Epetra_CrsGraph * NewGraph_;
  Epetra_CrsMatrix * NewMatrix_;

  Epetra_Map * UMap_;
  Epetra_Map * RMap_;
  Epetra_Map * LMap_;

  Epetra_Export * UExporter_;
  Epetra_Export * RExporter_;
  Epetra_Export * LExporter_;

  Epetra_MultiVector * ULHS_;
  Epetra_MultiVector * RLHS_;
  Epetra_MultiVector * LLHS_;

  Epetra_MultiVector * URHS_;
  Epetra_MultiVector * RRHS_;
  Epetra_MultiVector * LRHS_;

  Epetra_CrsGraph * UUGraph_;
  Epetra_CrsGraph * URGraph_;
  Epetra_CrsGraph * ULGraph_;
  Epetra_CrsGraph * RRGraph_;
  Epetra_CrsGraph * RLGraph_;
  Epetra_CrsGraph * LLGraph_;

  Epetra_CrsMatrix * UUMatrix_;
  Epetra_CrsMatrix * URMatrix_;
  Epetra_CrsMatrix * ULMatrix_;
  Epetra_CrsMatrix * RRMatrix_;
  Epetra_CrsMatrix * RLMatrix_;
  Epetra_CrsMatrix * LLMatrix_;

 public:

  ~LinearProblem_StaticCondensation();

  LinearProblem_StaticCondensation( int degree = 1, bool verbose = false )
  : degree_(degree),
    verbose_(verbose),
    OldRowMap_(0),
    OldProblem_(0),
    OldRHS_(0),
    OldLHS_(0),
    OldGraph_(0),
    OldMatrix_(0),
    Exporter_(0),
    NewRowMap_(0),
    NewColMap_(0),
    NewProblem_(0),
    NewRHS_(0),
    NewLHS_(0),
    NewGraph_(0),
    NewMatrix_(0),
    UMap_(0),
    RMap_(0),
    LMap_(0),
    UExporter_(0),
    RExporter_(0),
    LExporter_(0),
    ULHS_(0),
    RLHS_(0),
    LLHS_(0),
    URHS_(0),
    RRHS_(0),
    LRHS_(0),
    UUGraph_(0),
    URGraph_(0),
    ULGraph_(0),
    RRGraph_(0),
    RLGraph_(0),
    LLGraph_(0),
    UUMatrix_(0),
    URMatrix_(0),
    ULMatrix_(0),
    RRMatrix_(0),
    RLMatrix_(0),
    LLMatrix_(0)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();
};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_STATICCONDENSATION_H
