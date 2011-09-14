
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
