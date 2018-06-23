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

#ifdef HAVE_MUELU_STRATIMIKOS

// Stratimikos needs Thyra, so we don't need special guards for Thyra here
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_XpetraLinearOp.hpp"
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

#include <MueLu_Hierarchy.hpp>
#include <MueLu_HierarchyManager.hpp>
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_MasterList.hpp>
#include <MueLu_XpetraOperator_decl.hpp> // todo fix me
#include <MueLu_RefMaxwell.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <MueLu_TpetraOperator.hpp>
#endif
#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#include <Xpetra_EpetraOperator.hpp>
#endif

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Kokkos_DefaultNode.hpp"


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

      if (Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isBlockedOperator(fwdOp)) return true;

      return false;
    }

    /** \brief . */
    Teuchos::RCP<PreconditionerBase<Scalar> > createPrec() const {
      return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
    }

    /** \brief . */
    void initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc,
                        PreconditionerBase<Scalar>* prec,
                        const ESupportSolveUse supportSolveUse
                       ) const {
      using Teuchos::rcp_dynamic_cast;

      // we are using typedefs here, since we are using objects from different packages (Xpetra, Thyra,...)
      typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>                     XpMap;
      typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>      XpOp;
      typedef Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>       XpThyUtils;
      typedef Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>    XpCrsMatWrap;
      typedef Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>        XpCrsMat;
      typedef Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> XpBlockedCrsMat;
      typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>           XpMat;
      typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>      XpMultVec;
      typedef Xpetra::MultiVector<double,LocalOrdinal,GlobalOrdinal,Node>      XpMultVecDouble;
      typedef Thyra::LinearOpBase<Scalar>                                      ThyLinOpBase;
      typedef Thyra::DiagonalLinearOpBase<Scalar>                              ThyDiagLinOpBase;
#ifdef HAVE_MUELU_TPETRA
      // TAW 1/26/2016: We deal with Tpetra objects
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
      typedef Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> ThyTpLinOp;
#endif
#endif
#if defined(HAVE_MUELU_EPETRA)
      typedef Thyra::EpetraLinearOp                                         ThyEpLinOp;
      typedef Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>                  XpEpCrsMat;
#endif

      //std::cout << "-======---------------------------------" << std::endl;
      //std::cout << *paramList_ << std::endl;
      //std::cout << "-======---------------------------------" << std::endl;

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
        coordinates = MueLu::Utilities<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ExtractCoordinatesFromParameterList(paramList);
        paramList.set<RCP<XpMultVecDouble> >("Coordinates", coordinates);

        // TODO check for Xpetra or Thyra vectors?
#ifdef HAVE_MUELU_TPETRA
        if (bIsTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
          typedef Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>      tV;
          typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tMV;
          typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>      TpCrsMat;
          if (paramList.isType<Teuchos::RCP<tMV> >("Nullspace")) {
            RCP<tMV> tpetra_nullspace = paramList.get<RCP<tMV> >("Nullspace");
            paramList.remove("Nullspace");
            RCP<XpMultVec> nullspace = MueLu::TpetraMultiVector_To_XpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(tpetra_nullspace);
            paramList.set<RCP<XpMultVec> >("Nullspace", nullspace);
            TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(nullspace));
          }

          if (paramList.isParameter("M1")) {           
            if (paramList.isType<Teuchos::RCP<TpCrsMat> >("M1")) {
              RCP<TpCrsMat> tM1 = paramList.get<RCP<TpCrsMat> >("M1");
              paramList.remove("M1");
              RCP<XpCrsMat> xM1 = rcp_dynamic_cast<XpCrsMat>(tM1, true);
              paramList.set<RCP<XpCrsMat> >("M1", xM1);
            } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M1")) {
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
          } else
            TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Need to specify matrix M1.");

          if (paramList.isParameter("D0")) {
            if(paramList.isType<Teuchos::RCP<TpCrsMat> >("D0")) {
              RCP<TpCrsMat> tD0 = paramList.get<RCP<TpCrsMat> >("D0");
              paramList.remove("D0");
              RCP<XpCrsMat> xD0 = rcp_dynamic_cast<XpCrsMat>(tD0, true);
              paramList.set<RCP<XpCrsMat> >("D0", xD0);
            } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("D0")) {
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
              RCP<TpCrsMat> tM0inv = paramList.get<RCP<TpCrsMat> >("M0inv");
              paramList.remove("M0inv");
              RCP<XpCrsMat> xM0inv = rcp_dynamic_cast<XpCrsMat>(tM0inv, true);
              paramList.set<RCP<XpCrsMat> >("M0inv", xM0inv);
            } else if (paramList.isType<Teuchos::RCP<const ThyDiagLinOpBase> >("M0inv")) {
              RCP<const ThyDiagLinOpBase> thyM0inv = paramList.get<RCP<const ThyDiagLinOpBase> >("M0inv");
              paramList.remove("M0inv");
              RCP<const Thyra::VectorBase<Scalar> > diag = thyM0inv->getDiag();
              RCP<const Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ttDiag = rcp_dynamic_cast<const Thyra::TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(diag);
              RCP<const tV> tDiag = Thyra::TpetraOperatorVectorExtraction<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getConstTpetraVector(diag);
              RCP<XpMat> M0inv = Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Xpetra::toXpetra(tDiag));
              paramList.set<RCP<XpMat> >("M0inv", M0inv);
            } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M0inv")) {
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
#else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                                     "Thyra::MueLuRefMaxwellPreconditionerFactory: Tpetra does not support GO=int and or EpetraNode.");
#endif
        }
#endif
#ifdef HAVE_MUELU_EPETRA
        if (bIsEpetra) {
          if (paramList.isType<RCP<Epetra_MultiVector> >("Nullspace")) {
            RCP<Epetra_MultiVector> epetra_nullspace = Teuchos::null;
            epetra_nullspace = paramList.get<RCP<Epetra_MultiVector> >("Nullspace");
            paramList.remove("Nullspace");
            RCP<Xpetra::EpetraMultiVectorT<int,Node> > xpEpNullspace = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int,Node>(epetra_nullspace));
            RCP<Xpetra::MultiVector<double,int,int,Node> > xpEpNullspaceMult = rcp_dynamic_cast<Xpetra::MultiVector<double,int,int,Node> >(xpEpNullspace, true);
            RCP<XpMultVec> nullspace = rcp_dynamic_cast<XpMultVec>(xpEpNullspaceMult, true);
            paramList.set<RCP<XpMultVec> >("Nullspace", nullspace);
          }

          if (paramList.isParameter("M1")) {
            if (paramList.isType<Teuchos::RCP<Epetra_CrsMatrix> >("M1")) {
              RCP<Epetra_CrsMatrix> eM1 = paramList.get<RCP<Epetra_CrsMatrix> >("M1");
              paramList.remove("M1");
              RCP<XpEpCrsMat> xeM1 = Teuchos::rcp(new XpEpCrsMat(eM1));
              RCP<XpCrsMat> xCrsM1 = rcp_dynamic_cast<XpCrsMat>(xeM1, true);
              RCP<XpCrsMatWrap> xwM1 = Teuchos::rcp(new XpCrsMatWrap(xCrsM1));
              RCP<XpMat> xM1 = rcp_dynamic_cast<XpMat>(xwM1);
              paramList.set<RCP<XpMat> >("M1", xM1);
            }
            else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M1")) {
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
          } else
            TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "Need to specify matrix M1.");

          if (paramList.isParameter("D0")) {
            if (paramList.isType<Teuchos::RCP<Epetra_CrsMatrix> >("D0")) {
              RCP<Epetra_CrsMatrix> eD0 = paramList.get<RCP<Epetra_CrsMatrix> >("D0");
              paramList.remove("D0");
              RCP<XpEpCrsMat> xeD0 = Teuchos::rcp(new XpEpCrsMat(eD0));
              RCP<XpCrsMat> xCrsD0 = rcp_dynamic_cast<XpCrsMat>(xeD0, true);
              RCP<XpCrsMatWrap> xwD0 = Teuchos::rcp(new XpCrsMatWrap(xCrsD0));
              RCP<XpMat> xD0 = rcp_dynamic_cast<XpMat>(xwD0);
              paramList.set<RCP<XpMat> >("D0", xD0);
            }
            else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("D0")) {
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
            if (paramList.isType<Teuchos::RCP<Epetra_CrsMatrix> >("M0inv")) {
              RCP<Epetra_CrsMatrix> eM0inv = paramList.get<RCP<Epetra_CrsMatrix> >("M0inv");
              paramList.remove("M0inv");
              RCP<XpEpCrsMat> xeM0inv = Teuchos::rcp(new XpEpCrsMat(eM0inv));
              RCP<XpCrsMat> xCrsM0inv = rcp_dynamic_cast<XpCrsMat>(xeM0inv, true);
              RCP<XpCrsMatWrap> xwM0inv = Teuchos::rcp(new XpCrsMatWrap(xCrsM0inv));
              RCP<XpMat> xM0inv = rcp_dynamic_cast<XpMat>(xwM0inv);
              paramList.set<RCP<XpMat> >("M0inv", xM0inv);
            }
            else if (paramList.isType<Teuchos::RCP<const ThyDiagLinOpBase> >("M0inv")) {
              RCP<const ThyDiagLinOpBase> thyM0inv = paramList.get<RCP<const ThyDiagLinOpBase> >("M0inv");
              paramList.remove("M0inv");

              RCP<const Teuchos::Comm<int> > comm = A->getDomainMap()->getComm();
              RCP<const Epetra_Map> map = Thyra::get_Epetra_Map(*(thyM0inv->range()), Xpetra::toEpetra(comm));
              // RCP<XpMap> map = XpThyUtils::toXpetra(thyM0inv->range(), comm);
              RCP<const Thyra::VectorBase<double> > diag = thyM0inv->getDiag();
              RCP<const Epetra_Vector> eDiag = Thyra::get_Epetra_Vector(*map, diag);
              RCP<Epetra_Vector> nceDiag = Teuchos::rcp_const_cast<Epetra_Vector>(eDiag);
              RCP<Xpetra::EpetraVectorT<int,Node> > xpEpDiag = Teuchos::rcp(new Xpetra::EpetraVectorT<int,Node>(nceDiag));
              RCP<const Xpetra::Vector<double,int,int,Node> > xpDiag = rcp_dynamic_cast<const Xpetra::Vector<double,int,int,Node> >(xpEpDiag, true);
              RCP<XpMat> M0inv = Xpetra::MatrixFactory<double,int,int,Node>::Build(xpDiag);
              paramList.set<RCP<XpMat> >("M0inv", M0inv);
            } else if (paramList.isType<Teuchos::RCP<const ThyLinOpBase> >("M0inv")) {
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
        // build a new MueLu RefMaxwell preconditioner
        preconditioner = rcp(new MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A, paramList, true));

      } else {
        // reuse old MueLu preconditioner stored in MueLu Xpetra operator and put in new matrix

        // get old MueLu preconditioner
#if defined(HAVE_MUELU_TPETRA)
        if (bIsTpetra) {
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
          RCP<ThyTpLinOp> tpetr_precOp = rcp_dynamic_cast<ThyTpLinOp>(thyra_precOp);
          preconditioner = rcp_dynamic_cast<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tpetr_precOp->getTpetraOperator(),true);
#else
          TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError,
                                     "Thyra::MueLuRefMaxwellPreconditionerFactory: Tpetra does not support GO=int and or EpetraNode.");
#endif
        }
#endif
#if defined(HAVE_MUELU_EPETRA)// && defined(HAVE_MUELU_SERIAL)
        if (bIsEpetra) {
          RCP<ThyEpLinOp> epetr_precOp = rcp_dynamic_cast<ThyEpLinOp>(thyra_precOp);
          preconditioner = rcp_dynamic_cast<MueLu::RefMaxwell<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(epetr_precOp->epetra_op(),true);
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
