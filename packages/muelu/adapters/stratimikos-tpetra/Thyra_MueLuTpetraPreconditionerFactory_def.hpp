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
#ifndef THYRA_MUELU_TPETRA_PRECONDITIONER_FACTORY_DEF_HPP
#define THYRA_MUELU_TPETRA_PRECONDITIONER_FACTORY_DEF_HPP

#include "Thyra_MueLuTpetraPreconditionerFactory_decl.hpp"

#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

#include "TpetraCore_config.h"
#include "Tpetra_CrsMatrix.hpp"

#include "Teuchos_Ptr.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include <string>

namespace Thyra {

  using Teuchos::RCP;
  using Teuchos::ParameterList;


  // Constructors/initializers/accessors

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::MueLuTpetraPreconditionerFactory()
  {}

  // Overridden from PreconditionerFactoryBase

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  bool MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isCompatible(const LinearOpSourceBase<Scalar>& fwdOpSrc) const {
    const RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc.getOp();

    typedef Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyraTpetraLinOp;
    const RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);

    typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinOp;
    const RCP<const TpetraLinOp> tpetraFwdOp = Teuchos::nonnull(thyraTpetraFwdOp) ? thyraTpetraFwdOp->getConstTpetraOperator() : Teuchos::null;

    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMat;
    const RCP<const TpetraCrsMat> tpetraFwdCrsMat = Teuchos::rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);

    return Teuchos::nonnull(tpetraFwdCrsMat);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<PreconditionerBase<Scalar> > MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::createPrec() const {
    return Teuchos::rcp(new DefaultPreconditioner<Scalar>);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  initializePrec(const Teuchos::RCP<const LinearOpSourceBase<Scalar> >& fwdOpSrc, PreconditionerBase<Scalar>* prec, const ESupportSolveUse supportSolveUse) const {
    // Check precondition

    TEUCHOS_ASSERT(Teuchos::nonnull(fwdOpSrc));
    TEUCHOS_ASSERT(this->isCompatible(*fwdOpSrc));
    TEUCHOS_ASSERT(prec);

    Teuchos::Time totalTimer(""), timer("");
    totalTimer.start(true);

    const RCP<Teuchos::FancyOStream> out = this->getOStream();
    const Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
    Teuchos::OSTab tab(out);
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
      *out << "\nEntering Thyra::MueLuTpetraPreconditionerFactory::initializePrec(...) ...\n";
    }

    // Retrieve wrapped concrete Tpetra matrix from FwdOp

    const Teuchos::RCP<const LinearOpBase<Scalar> > fwdOp = fwdOpSrc->getOp();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(fwdOp));

    typedef Thyra::TpetraLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node> ThyraTpetraLinOp;
    const Teuchos::RCP<const ThyraTpetraLinOp> thyraTpetraFwdOp = Teuchos::rcp_dynamic_cast<const ThyraTpetraLinOp>(fwdOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thyraTpetraFwdOp));

    typedef Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraLinOp;
    const Teuchos::RCP<const TpetraLinOp> tpetraFwdOp = thyraTpetraFwdOp->getConstTpetraOperator();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdOp));

    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraCrsMat;
    const Teuchos::RCP<const TpetraCrsMat> tpetraFwdCrsMat = Teuchos::rcp_dynamic_cast<const TpetraCrsMat>(tpetraFwdOp);
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpetraFwdCrsMat));

    // Retrieve concrete preconditioner object

    const Teuchos::Ptr<DefaultPreconditioner<Scalar> > defaultPrec =
        Teuchos::ptr(dynamic_cast<DefaultPreconditioner<Scalar> *>(prec));
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(defaultPrec));

    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW))
      *out << "\nCreating a new MueLu::TpetraOperator object...\n";
    timer.start(true);

    // Workaround since MueLu interface does not accept const matrix as input
    const Teuchos::RCP<TpetraCrsMat> tpetraFwdCrsMatNonConst = Teuchos::rcp_const_cast<TpetraCrsMat>(tpetraFwdCrsMat);

    // Create and compute the initial preconditioner

    // Create a copy, as we may remove some things from the list
    Teuchos::ParameterList paramList = *paramList_;

    //Tpetra does not instantiate on Scalar=float by default, so we must check for this
    //FIXME This will still break if LO != int or GO != int
# if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_MUELU_INST_FLOAT_INT_INT)
    typedef Tpetra::MultiVector<float, LocalOrdinal, GlobalOrdinal, Node> fMV;
    Teuchos::RCP<fMV> floatCoords;
# endif
    typedef Tpetra::MultiVector<double, LocalOrdinal, GlobalOrdinal, Node> dMV;
    Teuchos::RCP<dMV> doubleCoords;
    if (paramList.isType<Teuchos::RCP<dMV> >("Coordinates")) {
      doubleCoords = paramList.get<Teuchos::RCP<dMV> >("Coordinates");
      paramList.remove("Coordinates");
    }
# if !defined(HAVE_TPETRA_EXPLICIT_INSTANTIATION) || defined(HAVE_MUELU_INST_FLOAT_INT_INT)
    else if (paramList.isType<Teuchos::RCP<fMV> >("Coordinates")) {
      floatCoords = paramList.get<Teuchos::RCP<fMV> >("Coordinates");
      paramList.remove("Coordinates");
      doubleCoords = Teuchos::rcp(new dMV(floatCoords->getMap(),floatCoords->getNumVectors()));
      deep_copy(*doubleCoords,*floatCoords);
    }
# endif

    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MV;
    Teuchos::RCP<MV> null_space;
    if (paramList.isType<Teuchos::RCP<MV> >("Nullspace")) {
      null_space = paramList.get<Teuchos::RCP<MV> >("Nullspace");
      paramList.remove("Nullspace");
    }

    typedef MueLu::TpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> MueLuOperator;
    const Teuchos::RCP<MueLuOperator> mueluPrecOp = MueLu::CreateTpetraPreconditioner(tpetraFwdCrsMatNonConst, paramList, doubleCoords, null_space);

    timer.stop();
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
      Teuchos::OSTab(out).o() << "> Creation time = " << timer.totalElapsedTime() << " sec\n";
    }

    const Teuchos::RCP<LinearOpBase<Scalar> > thyraPrecOp = Thyra::createLinearOp(Teuchos::RCP<TpetraLinOp>(mueluPrecOp));
    defaultPrec->initializeUnspecified(thyraPrecOp);

    totalTimer.stop();
    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_LOW)) {
      *out << "\nTotal time in Thyra::MueLuTpetraPreconditionerFactory::initializePrec(...) = " << totalTimer.totalElapsedTime() << " sec\n";
    }

    if (Teuchos::nonnull(out) && Teuchos::includesVerbLevel(verbLevel, Teuchos::VERB_MEDIUM)) {
      *out << "\nLeaving Thyra::MueLuTpetraPreconditionerFactory::initializePrec(...) ...\n";
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  uninitializePrec(PreconditionerBase<Scalar>* prec, Teuchos::RCP<const LinearOpSourceBase<Scalar> >* fwdOp, ESupportSolveUse* supportSolveUse) const {
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
  void MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setParameterList(Teuchos::RCP<ParameterList> const& paramList) {
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(paramList));
    paramList_ = paramList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getNonconstParameterList() {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<ParameterList> MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::unsetParameterList() {
    Teuchos::RCP<ParameterList> savedParamList = paramList_;
    paramList_ = Teuchos::null;
    return savedParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getParameterList() const {
    return paramList_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getValidParameters() const {
    static RCP<const ParameterList> validPL;

    if (Teuchos::is_null(validPL))
      validPL = Teuchos::rcp(new ParameterList());

    return validPL;
  }


  // Public functions overridden from Teuchos::Describable
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string MueLuTpetraPreconditionerFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const {
    return "Thyra::MueLuTpetraPreconditionerFactory";
  }

} // namespace Thyra

#endif // ifdef THYRA_MUELU_TPETRA_PRECONDITIONER_FACTORY_DEF_HPP
