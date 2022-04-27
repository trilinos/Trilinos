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
#include <list>

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

namespace Thyra {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_const_cast;

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

    return false;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<PreconditionerBase<Scalar> > MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const {
    return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuRefMaxwellPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  initializePrec(const RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc, PreconditionerBase<Scalar>* prec, const ESupportSolveUse supportSolveUse) const {

    // we are using typedefs here, since we are using objects from different packages (Xpetra, Thyra,...)
    typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>      XpOp;
    typedef Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>       XpThyUtils;
    typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>           XpMat;
    typedef Thyra::LinearOpBase<Scalar>                                      ThyLinOpBase;
    typedef Thyra::XpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyXpOp;
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
    typedef Xpetra::TpetraHalfPrecisionOperator<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpHalfPrecOp;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>           XpMV;
    typedef typename XpHalfPrecOp::HalfScalar                                     HalfScalar;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType                 Magnitude;
    typedef typename Teuchos::ScalarTraits<Magnitude>::halfPrecision              HalfMagnitude;
    typedef Xpetra::MultiVector<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>       XphMV;
    typedef Xpetra::MultiVector<Magnitude,LocalOrdinal,GlobalOrdinal,Node >       XpmMV;
    typedef Xpetra::MultiVector<HalfMagnitude,LocalOrdinal,GlobalOrdinal,Node >   XphmMV;
    typedef Xpetra::Matrix<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>            XphMat;
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
    TEUCHOS_TEST_FOR_EXCEPT((bIsEpetra == true  && bIsTpetra == true));

    // wrap the forward operator as an Xpetra::Matrix that MueLu can work with
    // MueLu needs a non-const object as input
    RCP<XpMat> A = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(fwdOp));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(A));

    // Retrieve concrete preconditioner object
    const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec = Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

    // extract preconditioner operator
    RCP<ThyLinOpBase> thyra_precOp = Teuchos::null;
    thyra_precOp = rcp_dynamic_cast<Thyra::LinearOpBase<Scalar> >(defaultPrec->getNonconstUnspecifiedPrecOp(), true);

    // make a decision whether to (re)build the multigrid preconditioner or reuse the old one
    // rebuild preconditioner if startingOver == true
    // reuse preconditioner if startingOver == false
    const bool startingOver = (thyra_precOp.is_null() || !paramList.isParameter("refmaxwell: enable reuse") || !paramList.get<bool>("refmaxwell: enable reuse"));
    const bool useHalfPrecision = paramList.get<bool>("half precision", false) && bIsTpetra;

    RCP<XpOp> xpPrecOp;
    if (startingOver == true) {

      // Convert to Xpetra
      std::list<std::string> convertXpetra = {"Coordinates", "Nullspace", "M1", "Ms", "D0", "M0inv"};
      for (auto it = convertXpetra.begin(); it != convertXpetra.end(); ++it)
        replaceWithXpetra<Scalar,LocalOrdinal,GlobalOrdinal,Node>(paramList,*it);

      paramList.set<bool>("refmaxwell: use as preconditioner", true);
      if (useHalfPrecision) {
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)

        // convert to half precision
        RCP<XphMat> halfA = Xpetra::convertToHalfPrecision(A);
        if (paramList.isType<RCP<XpmMV> >("Coordinates")) {
          RCP<XpmMV> coords = paramList.get<RCP<XpmMV> >("Coordinates");
          paramList.remove("Coordinates");
          RCP<XphmMV> halfCoords = Xpetra::convertToHalfPrecision(coords);
          paramList.set("Coordinates",halfCoords);
        }
        if (paramList.isType<RCP<XpMV> >("Nullspace")) {
          RCP<XpMV> nullspace = paramList.get<RCP<XpMV> >("Nullspace");
          paramList.remove("Nullspace");
          RCP<XphMV> halfNullspace = Xpetra::convertToHalfPrecision(nullspace);
          paramList.set("Nullspace",halfNullspace);
        }
        std::list<std::string> convertMat = {"M1", "Ms", "D0", "M0inv"};
        for (auto it = convertMat.begin(); it != convertMat.end(); ++it) {
          if (paramList.isType<RCP<XpMat> >(*it)) {
            RCP<XpMat> M = paramList.get<RCP<XpMat> >(*it);
            paramList.remove(*it);
            RCP<XphMat> halfM = Xpetra::convertToHalfPrecision(M);
            paramList.set(*it,halfM);
          }
        }

        // build a new half-precision MueLu RefMaxwell preconditioner
        RCP<MueLu::RefMaxwell<HalfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfPrec = rcp(new MueLu::RefMaxwell<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>(halfA, paramList, true));
        xpPrecOp = rcp(new XpHalfPrecOp(halfPrec));
#else
        TEUCHOS_TEST_FOR_EXCEPT(true);
#endif
      } else
        {
          // build a new MueLu RefMaxwell preconditioner
          RCP<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner = rcp(new MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, paramList, true));
          xpPrecOp = rcp_dynamic_cast<XpOp>(preconditioner);
        }
    } else {
      // reuse old MueLu preconditioner stored in MueLu Xpetra operator and put in new matrix

      RCP<ThyXpOp> thyXpOp = rcp_dynamic_cast<ThyXpOp>(thyra_precOp, true);
      RCP<XpOp>    xpOp    = thyXpOp->getXpetraOperator();
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
      RCP<XpHalfPrecOp> xpHalfPrecOp = rcp_dynamic_cast<XpHalfPrecOp>(xpOp);
      if (!xpHalfPrecOp.is_null()) {
        RCP<MueLu::RefMaxwell<HalfScalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner = rcp_dynamic_cast<MueLu::RefMaxwell<HalfScalar,LocalOrdinal,GlobalOrdinal,Node>>(xpHalfPrecOp->GetHalfPrecisionOperator(), true);
        RCP<XphMat> halfA = Xpetra::convertToHalfPrecision(A);
        preconditioner->resetMatrix(halfA);
        xpPrecOp = rcp_dynamic_cast<XpOp>(preconditioner);
      } else
#endif
        {
          RCP<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> > preconditioner = rcp_dynamic_cast<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(xpOp, true);
          preconditioner->resetMatrix(A);
          xpPrecOp = rcp_dynamic_cast<XpOp>(preconditioner);
        }
    }

    // wrap preconditioner in thyraPrecOp
    RCP<const VectorSpaceBase<Scalar> > thyraRangeSpace  = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpPrecOp->getRangeMap());
    RCP<const VectorSpaceBase<Scalar> > thyraDomainSpace = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpPrecOp->getDomainMap());

    RCP<ThyLinOpBase > thyraPrecOp = Thyra::xpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(thyraRangeSpace, thyraDomainSpace, xpPrecOp);
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
