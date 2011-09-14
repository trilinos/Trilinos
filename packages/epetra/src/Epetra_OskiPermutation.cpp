
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
