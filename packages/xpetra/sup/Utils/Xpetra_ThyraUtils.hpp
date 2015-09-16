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
#ifndef XPETRA_THYRAUTILS_HPP
#define XPETRA_THYRAUTILS_HPP

#ifdef HAVE_XPETRA_THYRA

#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Tpetra_ConfigDefs.hpp"
#endif

#ifdef HAVE_XPETRA_EPETRA
#include "Epetra_CombineMode.h"
#endif

#include "Xpetra_Map.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"
#include "Xpetra_MapExtractor.hpp"

#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>

//#include <Thyra_LinearOpBase.hpp>
#ifdef HAVE_XPETRA_TPETRA
//#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_SpmdVectorBase.hpp>
#include <Thyra_get_Epetra_Operator.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Xpetra_EpetraMap.hpp>
#endif

namespace Xpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class BlockedCrsMatrix;

template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = KokkosClassic::DefaultNode::DefaultNodeType>
class ThyraUtils {

private:
#undef XPETRA_THYRAUTILS_SHORT
#include "Xpetra_UseShortNames.hpp"

public:

  static Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> >
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, std::vector<size_t>& stridingInfo, LocalOrdinal stridedBlockId = -1, GlobalOrdinal offset = 0) {

    Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map = toXpetra(vectorSpace);

    if(stridedBlockId == -1) {
      TEUCHOS_TEST_FOR_EXCEPT(map->getNodeNumElements() % stridingInfo.size() != 0);
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPT(map->getNodeNumElements() % stridingInfo[stridedBlockId] != 0);
    }

    Teuchos::RCP<Xpetra::StridedMap<LocalOrdinal,GlobalOrdinal,Node> > ret = Xpetra::StridedMapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(map, stridingInfo, stridedBlockId, offset);
    return ret;
  }

  static Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  toXpetra(const Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >& vectorSpace, const Teuchos::RCP<const Teuchos::Comm<int> >& comm) {

    // check whether we have a Tpetra based Thyra operator
    bool bIsTpetra = false;
  #ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<const Thyra::TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetra_vsc = Teuchos::rcp_dynamic_cast<const Thyra::TpetraVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(vectorSpace);
    bIsTpetra = Teuchos::is_null(tpetra_vsc) ? false : true;
  #endif

    // check whether we have an Epetra based Thyra operator
    bool bIsEpetra = !bIsTpetra; // note: this is a little bit fragile!

  #ifdef HAVE_XPETRA_TPETRA
    if(bIsTpetra) {
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node> TOE;
      //Teuchos::RCP<const Thyra::TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgVec = vectorSpace->createMember();
      Teuchos::RCP<Thyra::VectorBase<Scalar> > rgVec = Thyra::createMember<Scalar>(vectorSpace, std::string("label"));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgVec));
      Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rgTpetraVec = TOE::getTpetraVector(rgVec);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraVec));
      Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rgTpetraMap = rgTpetraVec->getMap();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgTpetraMap));

      Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node > > rgXpetraMap = Xpetra::toXpetraNonConst(rgTpetraMap);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgXpetraMap));
      return rgXpetraMap;
    }
  #endif

  #ifdef HAVE_XPETRA_EPETRA
    if(bIsEpetra) {
      /*const Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<double> > spmd_vs =
        Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<double> >(vectorSpace, true);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(spmd_vs));*/
      //Teuchos::RCP<const Epetra_Comm> epComm = Thyra::get_Epetra_Comm(*spmd_vs->getComm());
      //Teuchos::RCP<const Epetra_Comm> epComm = Thyra::get_Epetra_Comm(*comm);
      Teuchos::RCP<const Epetra_Comm> epComm = Xpetra::toEpetra(comm);

      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epComm));
      //Teuchos::RCP<const Epetra_Comm> epComm = Thyra::get_Epetra_Comm(const Teuchos::Comm<Ordinal>& comm);
      Teuchos::RCP<const Epetra_Map> rgEpetraMap = Thyra::get_Epetra_Map(*vectorSpace, epComm);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgEpetraMap));

      Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node > > rgXpetraMap = Teuchos::rcp(new Xpetra::EpetraMap(rgEpetraMap));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgXpetraMap));
      return rgXpetraMap;
    }
  #endif
    return Teuchos::null;
  }

  static bool isTpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > & op){
    // check whether we have a Tpetra based Thyra operator
    bool bIsTpetra = false;
#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tpetra_op = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(op);
    bIsTpetra = Teuchos::is_null(tpetra_op) ? false : true;
#endif
    return bIsTpetra;
  }

  static bool isEpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > & op){
    // check whether we have an Epetra based Thyra operator
    bool bIsEpetra = false;
#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Thyra::EpetraLinearOp> epetra_op = Teuchos::rcp_dynamic_cast<const Thyra::EpetraLinearOp>(op);
    bIsEpetra = Teuchos::is_null(epetra_op) ? false : true;
#endif
    return bIsEpetra;
  }

  static Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  toXpetra(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >& op) {

#ifdef HAVE_XPETRA_TPETRA
    if(isTpetra(op)) {
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node> TOE;
      Teuchos::RCP<const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraOp = TOE::getConstTpetraOperator(op);
      // we should also add support for the const versions!
      //getConstTpetraOperator(const RCP<const LinearOpBase<Scalar> > &op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
      Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraRowMat = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraRowMat));
      Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraCrsMat = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraRowMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraCrsMat));
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraNcnstCrsMat = Teuchos::rcp_const_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraCrsMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraNcnstCrsMat));

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xTpetraCrsMat =
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(TpetraNcnstCrsMat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));
      return xTpetraCrsMat;
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if(isEpetra(op)) {
      Teuchos::RCP<const Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator( *op );
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<const Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(epetra_op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_rowmat));
      Teuchos::RCP<const Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(epetra_rowmat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_crsmat));
      Teuchos::RCP<Epetra_CrsMatrix> epetra_ncnstcrsmat = Teuchos::rcp_const_cast<Epetra_CrsMatrix>(epetra_crsmat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_ncnstcrsmat));

      Teuchos::RCP<Xpetra::EpetraCrsMatrix > xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epetra_ncnstcrsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));
      return xEpetraCrsMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  toXpetra(const Teuchos::RCP<Thyra::LinearOpBase<Scalar> >& op) {

#ifdef HAVE_XPETRA_TPETRA
    if(isTpetra(op)) {
      typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node> TOE;
      Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraOp = TOE::getTpetraOperator(op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraOp));
      Teuchos::RCP<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraRowMat = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraRowMat));
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > TpetraCrsMat = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(TpetraRowMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TpetraCrsMat));

      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xTpetraCrsMat =
          Teuchos::rcp(new Xpetra::TpetraCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(TpetraCrsMat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpetraCrsMat));
      return xTpetraCrsMat;
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    if(isEpetra(op)) {
      Teuchos::RCP<Epetra_Operator> epetra_op = Thyra::get_Epetra_Operator( *op );
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_op));
      Teuchos::RCP<Epetra_RowMatrix> epetra_rowmat = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(epetra_op);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_rowmat));
      Teuchos::RCP<Epetra_CrsMatrix> epetra_crsmat = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(epetra_rowmat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epetra_crsmat));

      Teuchos::RCP<Xpetra::EpetraCrsMatrix > xEpetraCrsMat =
          Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epetra_crsmat));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpetraCrsMat));
      return xEpetraCrsMat;
    }
#endif
    return Teuchos::null;
  }

  static Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
  toThyra(const Teuchos::RCP<const Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(mat);
    if(tpetraMat!=Teuchos::null) {
      Teuchos::RCP<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xTpCrsMat = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(mat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpCrsMat));
      Teuchos::RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpCrsMat = xTpCrsMat->getTpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpRowMat   = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpCrsMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpRowMat));
      Teuchos::RCP<const Tpetra::Operator <Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpOperator = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpRowMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpOperator));

      thyraOp = Thyra::createConstLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<const Xpetra::EpetraCrsMatrix> epetraMat = Teuchos::rcp_dynamic_cast<const Xpetra::EpetraCrsMatrix>(mat);
    if(epetraMat!=Teuchos::null) {
      Teuchos::RCP<Xpetra::EpetraCrsMatrix > xEpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrix >(mat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpCrsMat));
      Teuchos::RCP<const Epetra_CrsMatrix> epCrsMat = xEpCrsMat->getEpetra_CrsMatrix();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<const Thyra::EpetraLinearOp> thyraEpOp = Thyra::epetraLinearOp(epCrsMat,"op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
  toThyra(const Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& mat) {
    // create a Thyra operator from Xpetra::CrsMatrix
    Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thyraOp = Teuchos::null;

#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(mat);
    if(tpetraMat!=Teuchos::null) {
      Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xTpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(mat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xTpCrsMat));
      Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpCrsMat = xTpCrsMat->getTpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));

      Teuchos::RCP<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpRowMat   = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpCrsMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpRowMat));
      Teuchos::RCP<Tpetra::Operator <Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(tpRowMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpOperator));

      thyraOp = Thyra::createLinearOp(tpOperator);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraOp));
    }
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrix> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrix>(mat);
    if(epetraMat!=Teuchos::null) {
      Teuchos::RCP<Xpetra::EpetraCrsMatrix > xEpCrsMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrix >(mat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xEpCrsMat));
      Teuchos::RCP<Epetra_CrsMatrix> epCrsMat = xEpCrsMat->getEpetra_CrsMatrixNonConst();
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(epCrsMat));

      Teuchos::RCP<Thyra::EpetraLinearOp> thyraEpOp = Thyra::nonconstEpetraLinearOp(epCrsMat,"op");
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraEpOp));
      thyraOp = thyraEpOp;
    }
#endif
    return thyraOp;
  }

  static Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
  toThyra(const Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& mat) {

    int nRows = mat->Rows();
    int nCols = mat->Cols();

    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ablock = mat->getMatrix(0,0);

    bool bTpetra = false;
    bool bEpetra = false;
#ifdef HAVE_XPETRA_TPETRA
    Teuchos::RCP<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > tpetraMat = Teuchos::rcp_dynamic_cast<Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(Ablock);
    if(tpetraMat!=Teuchos::null) bTpetra = true;
#endif

#ifdef HAVE_XPETRA_EPETRA
    Teuchos::RCP<Xpetra::EpetraCrsMatrix> epetraMat = Teuchos::rcp_dynamic_cast<Xpetra::EpetraCrsMatrix>(Ablock);
    if(epetraMat!=Teuchos::null) bEpetra = true;
 #endif

    TEUCHOS_TEST_FOR_EXCEPT(bTpetra == bEpetra); // we only allow Epetra OR Tpetra

    // create new Thyra blocked operator
    Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<Scalar> > blockMat =
        Thyra::defaultBlockedLinearOp<Scalar>();

    blockMat->beginBlockFill(nRows,nCols);

    for (int r=0; r<nRows; ++r) {
      for (int c=0; c<nCols; ++c) {
        Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thBlock =
          Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(mat->getMatrix(r,c));
        std::stringstream label; label << "A" << r << c;
        blockMat->setBlock(r,c,thBlock);
      }
    }

    blockMat->endBlockFill();

    return blockMat;
  }

}; // end Utils class

} // end namespace Xpetra

#define XPETRA_THYRAUTILS_SHORT
#endif // HAVE_XPETRA_THYRA

#endif // XPETRA_THYRAUTILS_HPP
