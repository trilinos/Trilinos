
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
#include "Epetra_OskiMultiVector.h"

//=============================================================================

Epetra_OskiMultiVector::Epetra_OskiMultiVector(const Epetra_OskiMultiVector& Source) 
  : Epetra_MultiVector(Source), 
  Epetra_View_(Source.Epetra_View_), 
  Copy_Created_(Source.Copy_Created_) {
    Oski_View_ = oski_CopyVecView(Source.Oski_View_);
}

Epetra_OskiMultiVector::Epetra_OskiMultiVector(const Epetra_MultiVector& Source) 
  : Epetra_MultiVector(Source),
  Epetra_View_(&Source), 
  Copy_Created_(false) {
    double* A;
    double** Aptr;
    int LDA;
    int* LDAptr;
    LDAptr = new int[1];
    Aptr = new double*[1];
    if(Source.ConstantStride() || (Source.NumVectors() == 1)) {
      if(Source.ExtractView(Aptr, LDAptr))
        std::cerr << "Extract view failed\n";
      else
        Oski_View_ = oski_CreateMultiVecView(*Aptr, Source.MyLength(), Source.NumVectors(), LAYOUT_COLMAJ, *LDAptr);
    }
    else {
      Copy_Created_ = true;
      LDA = Source.MyLength();
      A = new double[LDA*Source.NumVectors()];
      if(Source.ExtractCopy(A, LDA))
        std::cerr << "Extract copy failed\n";
      else
        Oski_View_ = oski_CreateMultiVecView(A, Source.MyLength(), Source.NumVectors(), LAYOUT_COLMAJ, LDA);
    }
    delete [] LDAptr;
    delete [] Aptr;
}

Epetra_OskiMultiVector::~Epetra_OskiMultiVector() {
  if(oski_DestroyVecView(Oski_View_))
    std::cerr << "Vector destroy failed\n";
}

bool Epetra_OskiMultiVector::Copy_Created() const {
  return Copy_Created_;
}

oski_vecview_t Epetra_OskiMultiVector::Oski_View() const {
  return Oski_View_;
}

const Epetra_MultiVector* Epetra_OskiMultiVector::Epetra_View() const {
  return Epetra_View_;
}

Epetra_OskiMultiVector& Epetra_OskiMultiVector::operator = (const Epetra_OskiMultiVector& Source) {
  Epetra_View_ = Source.Epetra_View_;
  Oski_View_ = Source.Oski_View_;
  Copy_Created_ = Source.Copy_Created_;
  return(*this);
}

#endif
#endif
