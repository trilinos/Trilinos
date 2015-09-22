//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

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
