
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
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
// ************************************************************************
//@HEADER

#include "Epetra_ConfigDefs.h"

#ifdef HAVE_OSKI
#ifdef HAVE_EPETRA_TEUCHOS
#include "Epetra_OskiPermutation.h"

//=============================================================================


Epetra_OskiPermutation::Epetra_OskiPermutation() 
  : Permutation_(NULL) {}

Epetra_OskiPermutation::Epetra_OskiPermutation(const Epetra_OskiPermutation& Source) 
  : Permutation_(Source.Permutation_) {}

Epetra_OskiPermutation::Epetra_OskiPermutation(bool RowPerm, const Epetra_OskiMatrix& Source) {
  if(RowPerm)
    *this = Source.ViewRowPermutation();
  else
    *this = Source.ViewColumnPermutation();
}

Epetra_OskiPermutation::~Epetra_OskiPermutation() {
  if(Permutation_ != NULL)
    delete this;
}

void Epetra_OskiPermutation::ReplacePermutation (const oski_perm_t& InPerm) {
  Permutation_ = &InPerm;
}

int Epetra_OskiPermutation::PermuteVector(const bool TransA, Epetra_OskiMultiVector& Vector) const {
  int ReturnVal;
  if(TransA)
    ReturnVal = oski_PermuteVecView(*Permutation_, OP_TRANS, Vector.Oski_View());
  else
    ReturnVal = oski_PermuteVecView(*Permutation_, OP_NORMAL, Vector.Oski_View());
  if(ReturnVal)
    std::cerr << "Error in PermuteVector.\n";
  return ReturnVal;
}

#endif
#endif
