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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP
#define THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP

#include <MueLu_ConfigDefs.hpp>

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

// Stratimikos needs Thyra, so we don't need special guards for Thyra here
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
#include <Thyra_MueLuPreconditionerFactory.hpp>
#ifdef HAVE_MUELU_TPETRA
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#endif
#ifdef HAVE_MUELU_EPETRA
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_ThyraUtils.hpp>

#include <MueLu_XpetraOperator_decl.hpp> // todo fix me
#include <MueLu_RefMaxwell.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#include <Xpetra_TpetraOperator.hpp>
#include <Xpetra_TpetraHalfPrecisionOperator.hpp>
#endif
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#include <Xpetra_EpetraOperator.hpp>
#endif

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Kokkos_DefaultNode.hpp"

#include <list>

namespace Thyra {

  /** @brief Concrete preconditioner factory subclass for Thyra based on MueLu.
      @ingroup MueLuAdapters
      Add support for MueLu preconditioners in Thyra. This class provides an interface both
      for Epetra and Tpetra.

      The general implementation only handles Tpetra. For Epetra there is a specialization
      on SC=double, LO=int, GO=int and NO=EpetraNode.
  */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class MueLuRefMaxwellPreconditionerFactory : public PreconditionerFactoryBase<Scalar> {
  public:

    /** @name Constructors/initializers/accessors */
    //@{

    /** \brief . */
    MueLuRefMaxwellPreconditionerFactory();
    //@}

    /** @name Overridden from PreconditionerFactoryBase */
    //@{

    /** \brief . */
    bool isCompatible(const LinearOpSourceBase<Scalar>& fwdOp) const;
    /** \brief . */
    Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const;
    /** \brief . */
    void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOp,
                        PreconditionerBase<Scalar>* prec,
                        const ESupportSolveUse supportSolveUse
                       ) const;
    /** \brief . */
    void uninitializePrec(PreconditionerBase<Scalar>* prec,
                          Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp,
                          ESupportSolveUse* supportSolveUse
                         ) const;

    //@}

    /** @name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    /** \brief . */
    void                                          setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          unsetParameterList();
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          getNonconstParameterList();
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getParameterList() const;
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getValidParameters() const;
    //@}

    /** \name Public functions overridden from Describable. */
    //@{

    /** \brief . */
    std::string description() const;

    // ToDo: Add an override of describe(...) to give more detail!

    //@}

  private:

    Teuchos::RCP<Teuchos::ParameterList> paramList_;

  };

#ifdef HAVE_MUELU_EPETRA
  /** @brief Concrete preconditioner factory subclass for Thyra based on MueLu.
      @ingroup MueLuAdapters
      Add support for MueLu preconditioners in Thyra. This class provides an interface both
      for Epetra and Tpetra.

      Specialization for Epetra
  */
  template <>
  class MueLuRefMaxwellPreconditionerFactory<double,int,int,Xpetra::EpetraNode> : public PreconditionerFactoryBase<double> {
  public:
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Xpetra::EpetraNode Node;

    /** @name Constructors/initializers/accessors */
    //@{

    /** \brief . */
    MueLuRefMaxwellPreconditionerFactory()  : paramList_(rcp(new ParameterList())) { }

    //@}

    /** @name Overridden from PreconditionerFactoryBase */
    //@{

    /** \brief . */
    bool isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
      const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

#ifdef HAVE_MUELU_TPETRA
      if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isTpetra(fwdOp)) return true;
#endif

#ifdef HAVE_MUELU_EPETRA
      if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isEpetra(fwdOp)) return true;
#endif

      return false;
    }

    /** \brief . */
    Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const {
      return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
    }

    /** \brief . */
    void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc,
                        PreconditionerBase<Scalar>* prec,
                        const ESupportSolveUse /* supportSolveUse */
                       ) const {
      using Teuchos::rcp_dynamic_cast;

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

    /** \brief . */
    void uninitializePrec(PreconditionerBase<Scalar>* prec,
                          Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp,
                          ESupportSolveUse* supportSolveUse
                         ) const {
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

    //@}

    /** @name Overridden from Teuchos::ParameterListAcceptor */
    //@{

    /** \brief . */
    void                                          setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
      paramList_ = paramList;
    }
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          unsetParameterList() {
      RCP<ParameterList> savedParamList = paramList_;
      paramList_ = Teuchos::null;
      return savedParamList;
    }
    /** \brief . */
    Teuchos::RCP<Teuchos::ParameterList>          getNonconstParameterList() { return paramList_; }
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getParameterList() const {   return paramList_; }
    /** \brief . */
    Teuchos::RCP<const Teuchos::ParameterList>    getValidParameters() const {
      static RCP<const ParameterList> validPL;

      if (Teuchos::is_null(validPL))
        validPL = rcp(new ParameterList());

      return validPL;
    }
    //@}

    /** \name Public functions overridden from Describable. */
    //@{

    /** \brief . */
    std::string description() const { return "Thyra::MueLuRefMaxwellPreconditionerFactory"; }

    // ToDo: Add an override of describe(...) to give more detail!

    //@}

  private:
    Teuchos::RCP<Teuchos::ParameterList> paramList_;
  }; // end specialization for Epetra

#endif // HAVE_MUELU_EPETRA

} // namespace Thyra

#endif // #ifdef HAVE_MUELU_STRATIMIKOS

#endif // THYRA_MUELU_REFMAXWELL_PRECONDITIONER_FACTORY_DECL_HPP
