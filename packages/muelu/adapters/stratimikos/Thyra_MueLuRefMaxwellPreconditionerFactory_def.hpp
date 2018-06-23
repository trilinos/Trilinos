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
#ifndef THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DEF_HPP
#define THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DEF_HPP

#include "Thyra_MueLuRefMaxwellPreconditionerFactory_decl.hpp"

#ifdef HAVE_MUELU_STRATIMIKOS

namespace Thyra {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;


  // Constructors/initializers/accessors

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MueLuRefMaxwellPreconditionerFactory() :
      paramList_(rcp(new ParameterList()))
  {}

  // Overridden from PreconditionerFactoryBase

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
    const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

#ifdef HAVE_MUELU_TPETRA
    if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isTpetra(fwdOp)) return true;
#endif

    if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isBlockedOperator(fwdOp)) return true;

    return false;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<PreconditionerBase<Scalar> > MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const {
    return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
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
    typedef Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node>      XpMultVecDouble;
    typedef Thyra::LinearOpBase<Scalar>                                      ThyLinOpBase;
    typedef Thyra::DiagonalLinearOpBase<Scalar>                              ThyDiagLinOpBase;
#ifdef HAVE_MUELU_TPETRA
    typedef Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> ThyTpLinOp;
#endif
    Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec")));

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

    // Variable for RefMaxwell preconditioner: either build a new one or reuse the existing preconditioner
    RCP<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner = Teuchos::null;

    // make a decision whether to (re)build the multigrid preconditioner or reuse the old one
    // rebuild preconditioner if startingOver == true
    // reuse preconditioner if startingOver == false
    const bool startingOver = (thyra_precOp.is_null() || !paramList.isParameter("reuse: type") || paramList.get<std::string>("reuse: type") == "none");

    if (startingOver == true) {
      // extract coordinates from parameter list
      Teuchos::RCP<XpMultVecDouble> coordinates = Teuchos::null;
      {
        Teuchos::TimeMonitor tM_coords(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec get coords")));
        coordinates = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ExtractCoordinatesFromParameterList(paramList);
        paramList.set<RCP<XpMultVecDouble> >("Coordinates", coordinates);
      }

      // TODO check for Xpetra or Thyra vectors?
#ifdef HAVE_MUELU_TPETRA
      if (bIsTpetra) {
        typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>      tV;
        typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tMV;
        typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TpCrsMat;
        Teuchos::TimeMonitor tMwrap(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap objects")));
        if (paramList.isType<Teuchos::RCP<tMV> >("Nullspace")) {
          Teuchos::TimeMonitor tM_nullspace(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap nullspace")));
          RCP<tMV> tpetra_nullspace = paramList.get<RCP<tMV> >("Nullspace");
          paramList.remove("Nullspace");
          RCP<XpMultVec> nullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tpetra_nullspace);
          paramList.set<RCP<XpMultVec> >("Nullspace", nullspace);
          TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(nullspace));
        }

        if (paramList.isParameter("M1")) {
          if (paramList.isType<Teuchos::RCP<TpCrsMat> >("M1")) {
            Teuchos::TimeMonitor tM_M1(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap M1")));
            RCP<TpCrsMat> tM1 = paramList.get<RCP<TpCrsMat> >("M1");
            paramList.remove("M1");
            RCP<XpCrsMat> xM1 = rcp_dynamic_cast<XpCrsMat>(tM1, true);
            paramList.set<RCP<XpCrsMat> >("M1", xM1);
          } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M1")) {
            Teuchos::TimeMonitor tM_M1(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap M1")));
            RCP<const ThyLinOpBase> thyM1 = paramList.get<RCP<const ThyLinOpBase> >("M1");
            paramList.remove("M1");
            RCP<const XpCrsMat> crsM1 = XpThyUtils::toXpetra(thyM1);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsM1));
            // MueLu needs a non-const object as input
            RCP<XpCrsMat> crsM1NonConst = Teuchos::rcp_const_cast<XpCrsMat>(crsM1);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsM1NonConst));
            // wrap as an Xpetra::Matrix that MueLu can work with
            RCP<XpMat> M1 = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsM1NonConst));
            paramList.set<RCP<XpMat> >("M1", M1);
          }
        }  else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Need to specify matrix M1.");

        if (paramList.isParameter("D0")) {
          if (paramList.isType<Teuchos::RCP<TpCrsMat> >("D0")) {
            Teuchos::TimeMonitor tM_D0(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap D0")));
            RCP<TpCrsMat> tD0 = paramList.get<RCP<TpCrsMat> >("D0");
            paramList.remove("D0");
            RCP<XpCrsMat> xD0 = rcp_dynamic_cast<XpCrsMat>(tD0, true);
            paramList.set<RCP<XpCrsMat> >("D0", xD0);
          } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("D0")) {
            Teuchos::TimeMonitor tM_D0(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap D0")));
            RCP<const ThyLinOpBase> thyD0 = paramList.get<RCP<const ThyLinOpBase> >("D0");
            paramList.remove("D0");
            RCP<const XpCrsMat> crsD0 = XpThyUtils::toXpetra(thyD0);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsD0));
            // MueLu needs a non-const object as input
            RCP<XpCrsMat> crsD0NonConst = Teuchos::rcp_const_cast<XpCrsMat>(crsD0);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsD0NonConst));
            // wrap as an Xpetra::Matrix that MueLu can work with
            RCP<XpMat> D0 = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsD0NonConst));
            paramList.set<RCP<XpMat> >("D0", D0);
          }
        } else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Need to specify matrix D0.");

        if (paramList.isParameter("M0inv")) {
          if (paramList.isType<Teuchos::RCP<TpCrsMat> >("M0inv")) {
            Teuchos::TimeMonitor tM_M0inv(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap M0inv")));
            RCP<TpCrsMat> tM0inv = paramList.get<RCP<TpCrsMat> >("M0inv");
            paramList.remove("M0inv");
            RCP<XpCrsMat> xM0inv = rcp_dynamic_cast<XpCrsMat>(tM0inv, true);
            paramList.set<RCP<XpCrsMat> >("M0inv", xM0inv);
          } else if (paramList.isType<Teuchos::RCP<const ThyDiagLinOpBase> >("M0inv")) {
            Teuchos::TimeMonitor tM_M0inv(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap M0inv")));
            RCP<const ThyDiagLinOpBase> thyM0inv = paramList.get<RCP<const ThyDiagLinOpBase> >("M0inv");
            paramList.remove("M0inv");
            RCP<const Thyra::VectorBase<Scalar> > diag = thyM0inv->getDiag();
            RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraVector(diag);
            RCP<XpMat> M0inv = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Xpetra::toXpetra(tDiag));
            paramList.set<RCP<XpMat> >("M0inv", M0inv);
          } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M0inv")) {
            Teuchos::TimeMonitor tM_M0inv(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec wrap M0inv")));
            RCP<const ThyLinOpBase> thyM0inv = paramList.get<RCP<const ThyLinOpBase> >("M0inv");
            paramList.remove("M0inv");
            RCP<const XpCrsMat> crsM0inv = XpThyUtils::toXpetra(thyM0inv);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsM0inv));
            // MueLu needs a non-const object as input
            RCP<XpCrsMat> crsM0invNonConst = Teuchos::rcp_const_cast<XpCrsMat>(crsM0inv);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(crsM0invNonConst));
            // wrap as an Xpetra::Matrix that MueLu can work with
            RCP<XpMat> M0inv = rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(crsM0invNonConst));
            paramList.set<RCP<XpMat> >("M0inv", M0inv);
          }
        } else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Need to specify matrix M0inv.");
        
      }
#endif

      {
        // build a new MueLu RefMaxwell preconditioner
        Teuchos::TimeMonitor tMbuild(*Teuchos::TimeMonitor::getNewTimer(std::string("ThyraMueLuRefMaxwell::initializePrec build prec")));
        preconditioner = rcp(new MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, paramList, true));
      }

    } else {
      // reuse old MueLu preconditioner stored in MueLu Xpetra operator and put in new matrix

      // get old MueLu preconditioner
#if defined(HAVE_MUELU_TPETRA)
      if (bIsTpetra) {

        RCP<ThyTpLinOp> tpetr_precOp = rcp_dynamic_cast<ThyTpLinOp>(thyra_precOp);
        // RCP<XpOp>    muelu_precOp = rcp_dynamic_cast<XpOp>(tpetr_precOp->getTpetraOperator(),true);
        preconditioner = rcp_dynamic_cast<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpetr_precOp->getTpetraOperator(),true);
      }
#endif
      // TODO add the blocked matrix case here...

    }

    // wrap preconditioner in thyraPrecOp
    RCP<ThyLinOpBase > thyraPrecOp = Teuchos::null;
    RCP<const VectorSpaceBase<Scalar> > thyraRangeSpace  = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(preconditioner->getRangeMap());
    RCP<const VectorSpaceBase<Scalar> > thyraDomainSpace = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(preconditioner->getDomainMap());

    RCP<XpOp> xpOp = Teuchos::rcp_dynamic_cast<XpOp>(preconditioner);
    thyraPrecOp = Thyra::xpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(thyraRangeSpace, thyraDomainSpace,xpOp);

    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraPrecOp));

    defaultPrec->initializeUnspecified(thyraPrecOp);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
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
  void MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameterList(RCP<ParameterList> const& paramList) {
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
    paramList_ = paramList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstParameterList() {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unsetParameterList() {
    RCP<ParameterList> savedParamList = paramList_;
    paramList_ = Teuchos::null;
    return savedParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getParameterList() const {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getValidParameters() const {
    static RCP<const ParameterList> validPL;

    if (Teuchos::is_null(validPL))
      validPL = rcp(new ParameterList());

    return validPL;
  }

  // Public functions overridden from Teuchos::Describable
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
    return "Thyra::MueLuRefMaxwellPreconditionerFactory";
  }
} // namespace Thyra

#endif // HAVE_MUELU_STRATIMIKOS

#endif // ifdef THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DEF_HPP
