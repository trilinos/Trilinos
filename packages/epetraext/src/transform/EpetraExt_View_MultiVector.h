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

#ifndef EpetraExt_MULTIVECTOR_VIEW_H
#define EpetraExt_MULTIVECTOR_VIEW_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_Transform.h>

class Epetra_MultiVector;
class Epetra_BlockMap;

namespace EpetraExt {

//! Generates a sub-block view of a Epetra_MultiVector
class MultiVector_View : public ViewTransform<Epetra_MultiVector> {

  /* const Epetra_BlockMap & OrigMap_; (unused private field) */
  const Epetra_BlockMap & NewMap_;

  const int NumVec_;

 public:

  //! Destructor
  ~MultiVector_View();

  //! Constructor
  MultiVector_View( const Epetra_BlockMap & /* orig_map (unused) */,
                    const Epetra_BlockMap & new_map,
                    const int num_vec = -1 )
    : /* OrigMap_(orig_map), (unused) */
    NewMap_(new_map),
    NumVec_(num_vec)
  {}

  //! Transformation Operator
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_MULTIVECTOR_VIEW_H
