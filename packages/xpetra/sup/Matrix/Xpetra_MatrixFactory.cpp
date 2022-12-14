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

#ifndef XPETRA_MATRIXFACTORY_CPP
#define XPETRA_MATRIXFACTORY_CPP

#include "Xpetra_MatrixFactory.hpp"
#if 0
namespace Xpetra {
   template<>
   Teuchos::RCP<Matrix > MatrixFactory2<double,int,int,Xpetra::Matrix<double, int, int>::node_type>::BuildCopy(const Teuchos::RCP<const Matrix > A) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
        RCP<const EpetraCrsMatrixT<GlobalOrdinal,Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal,Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix>     newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal,Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap  (newECrsOp));

      return newOp;
    }
#endif
#endif

    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix>     newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap(newTCrsOp));

      return newOp;
    }

    return Teuchos::null;  // make compiler happy
  }

#ifdef HAVE_XPETRA_INT_LONG_LONG
  template<>
  Teuchos::RCP<Xpetra::Matrix<double,int,long long,typename Xpetra::Matrix<double, int, long long>::node_type> > MatrixFactory2<double,int,long long,typename Xpetra::Matrix<double, int, long long>::node_type>::BuildCopy(const Teuchos::RCP<const Xpetra::Matrix<double,int,long long,typename Matrix<double, int, long long>::node_type> > A) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
        RCP<const EpetraCrsMatrixT<GlobalOrdinal,Node> > oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrixT<GlobalOrdinal,Node> >(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix>     newECrsOp(new EpetraCrsMatrixT<GlobalOrdinal,Node>(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap  (newECrsOp));

      return newOp;
    }
#endif
#endif

    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix>     newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap(newTCrsOp));

      return newOp;
    }

    return Teuchos::null;  // make compiler happy
  }

#endif // HAVE_XPETRA_INT_LONG_LONG

} // namespace Xpetra
#endif // if 0
#endif
