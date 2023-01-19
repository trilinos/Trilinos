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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#include "Xpetra_ConfigDefs.hpp"
#ifdef HAVE_XPETRA_THYRA

#include "Xpetra_BlockedCrsMatrix.hpp"

//#include "Xpetra_ThyraUtils.hpp"


namespace Xpetra {

#ifdef HAVE_XPETRA_EPETRA

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
  // implementation of "toThyra" for full specialization on SC=double, LO=GO=int and NO=EpetraNode
  // We need the specialization in the cpp file due to a circle dependency in the .hpp files for BlockedCrsMatrix
  Teuchos::RCP<Thyra::LinearOpBase<double> >
  ThyraUtils<double, int, int, EpetraNode>::toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, int, EpetraNode> >& mat) {

    int nRows = mat->Rows();
    int nCols = mat->Cols();

    Teuchos::RCP<Xpetra::Matrix<double, int, int, EpetraNode> > Ablock = mat->getInnermostCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrixWrap<double, int, int, EpetraNode> > Ablock_wrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<double, int, int, EpetraNode>>(Ablock);
    TEUCHOS_TEST_FOR_EXCEPT(Ablock_wrap.is_null() == true);

    bool bTpetra = false;
    bool bEpetra = false;
    // Note: Epetra is enabled
#if ((defined(EPETRA_HAVE_OMP)  && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT) && defined(HAVE_TPETRA_INST_DOUBLE)))
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Ablock_wrap->getCrsMatrix());
    if(tpetraMat!=Teuchos::null) bTpetra = true;
#else
    bTpetra = false;
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node> > epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node> >(Ablock_wrap->getCrsMatrix());
    if(epetraMat!=Teuchos::null) bEpetra = true;
#endif

    TEUCHOS_TEST_FOR_EXCEPT(bTpetra == bEpetra); // we only allow Epetra OR Tpetra

    // create new Thyra blocked operator
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> > blockMat =
        Thyra::defaultBlockedLinearOp<Scalar>();

    blockMat->beginBlockFill(nRows,nCols);

    for (int r=0; r<nRows; ++r) {
      for (int c=0; c<nCols; ++c) {
        Teuchos::RCP<Matrix> xpmat = mat->getMatrix(r,c);

        if(xpmat == Teuchos::null) continue; // shortcut for empty blocks

        Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thBlock = Teuchos::null;

        // check whether the subblock is again a blocked operator
        Teuchos::RCP<BlockedCrsMatrix> xpblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(xpmat);
        if(xpblock != Teuchos::null) {
          if(xpblock->Rows() == 1 && xpblock->Cols() == 1) {
            // If it is a single block operator, unwrap it
            Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpblock->getCrsMatrix());
            TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
            thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpwrap->getCrsMatrix());
          } else {
            // recursive call for general blocked operators
            thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpblock);
          }
        } else {
          // check whether it is a CRSMatrix object
          Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpmat);
          TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
          thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpwrap->getCrsMatrix());
        }

        blockMat->setBlock(r,c,thBlock);
      }
    }

    blockMat->endBlockFill();

    return blockMat;
  }
# endif //#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
  // implementation of "toThyra" for full specialization on SC=double, LO=int, GO=long long and NO=EpetraNode
  // We need the specialization in the cpp file due to a circle dependency in the .hpp files for BlockedCrsMatrix
  Teuchos::RCP<Thyra::LinearOpBase<double> >
  ThyraUtils<double, int, long long, EpetraNode>::toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<double, int, long long, EpetraNode> >& mat) {

    int nRows = mat->Rows();
    int nCols = mat->Cols();

    Teuchos::RCP<Xpetra::Matrix<double, int, long long, EpetraNode> > Ablock = mat->getInnermostCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrixWrap<double, int, long long, EpetraNode> > Ablock_wrap = Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<double, int, long long, EpetraNode>>(Ablock);
    TEUCHOS_TEST_FOR_EXCEPT(Ablock_wrap.is_null() == true);

    bool bTpetra = false;
    bool bEpetra = false;
    // Note: Epetra is enabled
#if ((defined(EPETRA_HAVE_OMP)  && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)) || \
     (!defined(EPETRA_HAVE_OMP) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG) && defined(HAVE_TPETRA_INST_DOUBLE)))
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Ablock_wrap->getCrsMatrix());
    if(tpetraMat!=Teuchos::null) bTpetra = true;
#else
    bTpetra = false;
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node> > epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node> >(Ablock_wrap->getCrsMatrix());
    if(epetraMat!=Teuchos::null) bEpetra = true;
#endif

    TEUCHOS_TEST_FOR_EXCEPT(bTpetra == bEpetra); // we only allow Epetra OR Tpetra

    // create new Thyra blocked operator
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> > blockMat =
        Thyra::defaultBlockedLinearOp<Scalar>();

    blockMat->beginBlockFill(nRows,nCols);

    for (int r=0; r<nRows; ++r) {
      for (int c=0; c<nCols; ++c) {
        Teuchos::RCP<Matrix> xpmat = mat->getMatrix(r,c);

        if(xpmat == Teuchos::null) continue; // shortcut for empty blocks

        Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thBlock = Teuchos::null;

        // check whether the subblock is again a blocked operator
        Teuchos::RCP<BlockedCrsMatrix> xpblock = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(xpmat);
        if(xpblock != Teuchos::null) {
          if(xpblock->Rows() == 1 && xpblock->Cols() == 1) {
            // If it is a single block operator, unwrap it
            Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpblock->getCrsMatrix());
            TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
            thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpwrap->getCrsMatrix());
          } else {
            // recursive call for general blocked operators
            thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpblock);
          }
        } else {
          // check whether it is a CRSMatrix object
          Teuchos::RCP<CrsMatrixWrap> xpwrap = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xpmat);
          TEUCHOS_TEST_FOR_EXCEPT(xpwrap.is_null() == true);
          thBlock = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpwrap->getCrsMatrix());
        }

        blockMat->setBlock(r,c,thBlock);
      }
    }

    blockMat->endBlockFill();

    return blockMat;
  }
#endif // #ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES

#endif

} // namespace Xpetra

#endif // HAVE_XPETRA_THYRA
