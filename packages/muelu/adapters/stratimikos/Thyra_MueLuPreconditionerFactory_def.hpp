// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
#ifndef THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP
#define THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP

#include "Thyra_MueLuPreconditionerFactory_decl.hpp"

#ifdef HAVE_MUELU_STRATIMIKOS

namespace Thyra {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;


  // Constructors/initializers/accessors

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MueLuPreconditionerFactory() :
      paramList_(rcp(new ParameterList()))
  {}

  // Overridden from PreconditionerFactoryBase

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
    const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

#ifdef HAVE_MUELU_TPETRA
    if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isTpetra(fwdOp)) return true;
#endif

    if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isBlockedOperator(fwdOp)) return true;

    return false;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<PreconditionerBase<Scalar> > MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const {
    return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  initializePrec(const RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc, PreconditionerBase<Scalar>* prec, const ESupportSolveUse supportSolveUse) const {
    using Teuchos::rcp_dynamic_cast;

    // we are using typedefs here, since we are using objects from different packages (Xpetra, Thyra,...)
    typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                     XpMap;
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>      XpOp;
    typedef Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>       XpThyUtils;
    typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XpCrsMat;
    typedef Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpBlockedCrsMat;
    typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>           XpMat;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      XpMultVec;
    typedef Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType,LocalOrdinal,GlobalOrdinal,Node>      XpMultVecDouble;
    typedef Thyra::LinearOpBase<Scalar>                                      ThyLinOpBase;
#ifdef HAVE_MUELU_TPETRA
    typedef MueLu::TpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> MueTpOp;
    typedef Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TpOp;
    typedef Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> ThyTpLinOp;
#endif

    // Check precondition
    TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
    TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
    TEUCHOS_ASSERT(prec);

    // Create a copy, as we may remove some things from the list
    ParameterList paramList = *paramList_;

    // Retrieve wrapped concrete Xpetra matrix from FwdOp
    const RCP<const ThyLinOpBase> fwdOp = fwdOpSrc->getOp();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

    // Check whether it is Epetra/Tpetra
    bool bIsEpetra  = XpThyUtils::isEpetra(fwdOp);
    bool bIsTpetra  = XpThyUtils::isTpetra(fwdOp);
    bool bIsBlocked = XpThyUtils::isBlockedOperator(fwdOp);
    TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == true  && bIsTpetra == true));
    TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == bIsTpetra) && bIsBlocked == false);
    TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra != bIsTpetra) && bIsBlocked == true);

    RCP<XpMat> A = Teuchos::null;
    if(bIsBlocked) {
      Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > ThyBlockedOp =
          Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase<Scalar> >(fwdOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(ThyBlockedOp));

      TEUCHOS_TEST_FOR_EXCEPT(ThyBlockedOp->blockExists(0,0)==false);

      Teuchos::RCP<const LinearOpBase<Scalar> > b00 = ThyBlockedOp->getBlock(0,0);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(b00));

      RCP<const XpCrsMat > xpetraFwdCrsMat00 = XpThyUtils::toXpetra(b00);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMat00));

      // MueLu needs a non-const object as input
      RCP<XpCrsMat> xpetraFwdCrsMatNonConst00 = Teuchos::rcp_const_cast<XpCrsMat>(xpetraFwdCrsMat00);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMatNonConst00));

      // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
      RCP<XpMat> A00 = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xpetraFwdCrsMatNonConst00));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A00));

      RCP<const XpMap> rowmap00 = A00->getRowMap();
      RCP< const Teuchos::Comm< int > > comm = rowmap00->getComm();

      // create a Xpetra::BlockedCrsMatrix which derives from Xpetra::Matrix that MueLu can work with
      RCP<XpBlockedCrsMat> bMat = Teuchos::rcp(new XpBlockedCrsMat(ThyBlockedOp, comm));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(bMat));

      // save blocked matrix
      A = bMat;
    } else {
      RCP<const XpCrsMat > xpetraFwdCrsMat = XpThyUtils::toXpetra(fwdOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMat));

      // MueLu needs a non-const object as input
      RCP<XpCrsMat> xpetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<XpCrsMat>(xpetraFwdCrsMat);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(xpetraFwdCrsMatNonConst));

      // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
      A = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xpetraFwdCrsMatNonConst));
    }
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A));

    // Retrieve concrete preconditioner object
    const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

    // extract preconditioner operator
    RCP<ThyLinOpBase> thyra_precOp = Teuchos::null;
    thyra_precOp = rcp_dynamic_cast<Thyra::LinearOpBase<Scalar> >(defaultPrec->getNonconstUnspecifiedPrecOp(), true);

    // Variable for multigrid hierarchy: either build a new one or reuse the existing hierarchy
    RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > H = Teuchos::null;

    // make a decision whether to (re)build the multigrid preconditioner or reuse the old one
    // rebuild preconditioner if startingOver == true
    // reuse preconditioner if startingOver == false
    const bool startingOver = (thyra_precOp.is_null() || !paramList.isParameter("reuse: type") || paramList.get<std::string>("reuse: type") == "none");

    if (startingOver == true) {
      // extract coordinates from parameter list
      Teuchos::RCP<XpMultVecDouble> coordinates = Teuchos::null;
      coordinates = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ExtractCoordinatesFromParameterList(paramList);

      // TODO check for Xpetra or Thyra vectors?
      RCP<XpMultVec> nullspace = Teuchos::null;
#ifdef HAVE_MUELU_TPETRA
      if (bIsTpetra) {
        typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tMV;
        RCP<tMV> tpetra_nullspace = Teuchos::null;
        if (paramList.isType<Teuchos::RCP<tMV> >("Nullspace")) {
          tpetra_nullspace = paramList.get<RCP<tMV> >("Nullspace");
          paramList.remove("Nullspace");
          nullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tpetra_nullspace);
          TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(nullspace));
        }
      }
#endif
      // build a new MueLu hierarchy
      H = MueLu::CreateXpetraPreconditioner(A, paramList, coordinates, nullspace);

    } else {
      // reuse old MueLu hierarchy stored in MueLu Tpetra/Epetra operator and put in new matrix

      // get old MueLu hierarchy
#if defined(HAVE_MUELU_TPETRA)
      if (bIsTpetra) {

        RCP<ThyTpLinOp> tpetr_precOp = rcp_dynamic_cast<ThyTpLinOp>(thyra_precOp);
        RCP<MueTpOp>    muelu_precOp = rcp_dynamic_cast<MueTpOp>(tpetr_precOp->getTpetraOperator(),true);

        H = muelu_precOp->GetHierarchy();
      }
#endif
      // TODO add the blocked matrix case here...

      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetNumLevels(), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no levels in it");
      TEUCHOS_TEST_FOR_EXCEPTION(!H->GetLevel(0)->IsAvailable("A"), MueLu::Exceptions::RuntimeError,
                                 "Thyra::MueLuPreconditionerFactory: Hierarchy has no fine level operator");
      RCP<MueLu::Level> level0 = H->GetLevel(0);
      RCP<XpOp>    O0 = level0->Get<RCP<XpOp> >("A");
      RCP<XpMat>   A0 = rcp_dynamic_cast<XpMat>(O0);

      if (!A0.is_null()) {
        // If a user provided a "number of equations" argument in a parameter list
        // during the initial setup, we must honor that settings and reuse it for
        // all consequent setups.
        A->SetFixedBlockSize(A0->GetFixedBlockSize());
      }

      // set new matrix
      level0->Set("A", A);

      H->SetupRe();
    }

    // wrap hierarchy H in thyraPrecOp
    RCP<ThyLinOpBase > thyraPrecOp = Teuchos::null;
#if defined(HAVE_MUELU_TPETRA)
    if (bIsTpetra) {
      RCP<MueTpOp> muelu_tpetraOp = rcp(new MueTpOp(H));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(muelu_tpetraOp));
      RCP<TpOp> tpOp = Teuchos::rcp_dynamic_cast<TpOp>(muelu_tpetraOp);
      thyraPrecOp = Thyra::createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpOp);
    }
#endif

    if(bIsBlocked) {
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::nonnull(thyraPrecOp));

      typedef MueLu::XpetraOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node>    MueXpOp;
      //typedef Thyra::XpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>    ThyXpLinOp; // unused
      const RCP<MueXpOp> muelu_xpetraOp = rcp(new MueXpOp(H));

      RCP<const VectorSpaceBase<Scalar> > thyraRangeSpace  = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(muelu_xpetraOp->getRangeMap());
      RCP<const VectorSpaceBase<Scalar> > thyraDomainSpace = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(muelu_xpetraOp->getDomainMap());

      RCP <Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpOp = Teuchos::rcp_dynamic_cast<Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(muelu_xpetraOp);
      thyraPrecOp = Thyra::xpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(thyraRangeSpace, thyraDomainSpace,xpOp);
    }

    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraPrecOp));

    defaultPrec->initializeUnspecified(thyraPrecOp);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  uninitializePrec(PreconditionerBase<Scalar>* prec, RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
    TEUCHOS_ASSERT(prec);

    // Retrieve concrete preconditioner object
    const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

    if (fwdOp) {
      // TODO: Implement properly instead of returning default value
      *fwdOp = Teuchos::null;
    }

    if (supportSolveUse) {
      // TODO: Implement properly instead of returning default value
      *supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED;
    }

    defaultPrec->uninitialize();
  }


  // Overridden from ParameterListAcceptor
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameterList(RCP<ParameterList> const& paramList) {
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
    paramList_ = paramList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstParameterList() {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unsetParameterList() {
    RCP<ParameterList> savedParamList = paramList_;
    paramList_ = Teuchos::null;
    return savedParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getParameterList() const {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getValidParameters() const {
    static RCP<const ParameterList> validPL;

    if (Teuchos::is_null(validPL))
      validPL = rcp(new ParameterList());

    return validPL;
  }

  // Public functions overridden from Teuchos::Describable
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MueLuPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
    return "Thyra::MueLuPreconditionerFactory";
  }
} // namespace Thyra

#endif // HAVE_MUELU_STRATIMIKOS

#endif // ifdef THYRA_MUELU_PRECONDITIONER_FACTORY_DEF_HPP
