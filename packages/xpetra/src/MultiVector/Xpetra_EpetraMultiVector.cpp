// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_EpetraMultiVector.hpp"

#include "Xpetra_EpetraImport.hpp"
#include "Xpetra_EpetraExport.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_EpetraVector.hpp"

#include "Epetra_SerialComm.h"

namespace Xpetra {

// TODO: move that elsewhere
  template<class GlobalOrdinal>
  const Epetra_MultiVector & toEpetra(const MultiVector<double, int, GlobalOrdinal> & x) {
    XPETRA_DYNAMIC_CAST(const EpetraMultiVectorT<GlobalOrdinal>, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }

  template<class GlobalOrdinal>
  Epetra_MultiVector & toEpetra(MultiVector<double, int, GlobalOrdinal> & x) {
    XPETRA_DYNAMIC_CAST(      EpetraMultiVectorT<GlobalOrdinal>, x, tX, "toEpetra");
    return *tX.getEpetra_MultiVector();
  }
  //

  template<class GlobalOrdinal>
  RCP<MultiVector<double, int, GlobalOrdinal> > toXpetra(RCP<Epetra_MultiVector> vec) {
    if (!vec.is_null())
      return rcp(new EpetraMultiVectorT<GlobalOrdinal>(vec));

    return Teuchos::null;
  }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template RCP<MultiVector<double, int, int> > toXpetra<int>(RCP<Epetra_MultiVector>);
template const Epetra_MultiVector & toEpetra(const MultiVector<double, int, int> &);
template Epetra_MultiVector & toEpetra(MultiVector<double, int, int> &);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template RCP<MultiVector<double, int, long long> > toXpetra<long long>(RCP<Epetra_MultiVector>);
template const Epetra_MultiVector & toEpetra(const MultiVector<double, int, long long> &);
template Epetra_MultiVector & toEpetra(MultiVector<double, int, long long> &);
#endif

} // namespace Xpetra
