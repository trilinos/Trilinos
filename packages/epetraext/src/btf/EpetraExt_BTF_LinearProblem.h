// @HEADER
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
// @HEADER

#ifndef EpetraExt_LINEARPROBLEM_BTF_H
#define EpetraExt_LINEARPROBLEM_BTF_H

#include <EpetraExt_Transform.h>

#include <vector>
#include <map>
#include <set>

class Epetra_LinearProblem;
class Epetra_VbrMatrix;
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_CrsGraph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_SerialDenseMatrix;

namespace EpetraExt {

class LinearProblem_BTF : public SameTypeTransform<Epetra_LinearProblem> {

 public:

  ~LinearProblem_BTF();

  LinearProblem_BTF( double thres = 0.0,
                     int verbose = 0 )
  : NewMap_(0),
    NewMatrix_(0),
    NewGraph_(0),
    NewLHS_(0),
    NewRHS_(0),
    NewProblem_(0),
    OrigGraph_(0),
    OrigMatrix_(0),
    OrigRHS_(0),
    OrigLHS_(0),
    OrigProblem_(0),
    threshold_(thres),
    verbose_(verbose),
    changedLP_(false)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

  bool changedLP() { return changedLP_; }

 private:

  void deleteNewObjs_();

  Epetra_BlockMap * NewMap_;
  
  Epetra_LinearProblem * NewProblem_;

  Epetra_VbrMatrix * NewMatrix_;
  Epetra_CrsGraph * NewGraph_;

  Epetra_MultiVector * NewLHS_;
  Epetra_MultiVector * NewRHS_;

  Epetra_Map * OrigRowMap_;
  Epetra_Map * OrigColMap_;
  Epetra_LinearProblem * OrigProblem_;
  Epetra_CrsGraph * OrigGraph_;
  Epetra_CrsMatrix * OrigMatrix_;
  Epetra_MultiVector * OrigLHS_;
  Epetra_MultiVector * OrigRHS_;

  std::vector<int> OldGlobalElements_;

  std::vector< std::set<int> > ZeroElements_;

  std::vector< std::vector<Epetra_SerialDenseMatrix*> > Blocks_;
  std::vector<int> BlockDim_;
  std::vector<int> BlockCnt_;
  std::map<int,int> BlockRowMap_;
  std::map<int,int> SubBlockRowMap_;
  std::map<int,int> BlockColMap_;
  std::map<int,int> SubBlockColMap_;

  std::vector< std::vector<int> > NewBlockRows_;

  const double threshold_;
  const int verbose_;

  bool changedLP_;
};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_BTF_H
