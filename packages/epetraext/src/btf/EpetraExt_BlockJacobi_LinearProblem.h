//@HEADER
// ************************************************************************
//
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
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

#ifndef EPETRAEXT_BLOCKJACOBI_LINEARPROBLEM_H
#define EPETRAEXT_BLOCKJACOBI_LINEARPROBLEM_H

#include <EpetraExt_Transform.h>

#include <vector>

class Epetra_LinearProblem;
class Epetra_VbrMatrix;
class Epetra_SerialDenseSVD;
class Epetra_SerialDenseMatrix;

namespace EpetraExt {

class LinearProblem_BlockJacobi : public SameTypeTransform<Epetra_LinearProblem> {

 public:

  ~LinearProblem_BlockJacobi();

  LinearProblem_BlockJacobi( int verbose = 0,
                             int thresholding = 0,
                             double rthresh = 0.0,
                             double athresh = 0.0,
                             bool removeDiag = false )
  : NumBlocks_(0),
    NewProblem_(0),
    NewMatrix_(0),
    thresholding_(thresholding),
    rthresh_(rthresh),
    athresh_(athresh),
    removeDiag_(removeDiag),
    verbose_(verbose)
  {}

  NewTypeRef operator()( OriginalTypeRef orig );

  bool fwd();
  bool rvs();

  void RThresh( double val ) { rthresh_ = val; }

 private:

  int NumBlocks_;

  double rthresh_;
  double athresh_;
  const int thresholding_;

  const bool removeDiag_;

  Epetra_LinearProblem * NewProblem_;
  Epetra_VbrMatrix * NewMatrix_;

  std::vector<Epetra_SerialDenseMatrix**> VbrBlocks_;
  std::vector<int> VbrBlockCnt_;
  std::vector<int> VbrBlockDim_;
  std::vector<int*> VbrBlockIndices_;

  std::vector<Epetra_SerialDenseSVD*> SVDs_;
  std::vector<Epetra_SerialDenseMatrix*> Inverses_;
  std::vector<Epetra_SerialDenseMatrix*> RHSBlocks_;

  const int verbose_;
};

} //namespace EpetraExt

#endif //EPETRAEXT_BLOCKJACOBI_LINEARPROBLEM_H
