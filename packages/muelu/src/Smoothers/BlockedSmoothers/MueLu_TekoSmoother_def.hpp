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
#ifndef MUELU_TEKOSMOOTHER_DEF_HPP_
#define MUELU_TEKOSMOOTHER_DEF_HPP_

#ifdef HAVE_MUELU_TEKO

#include "Teuchos_ScalarTraits.hpp"

#include <Thyra_DefaultProductMultiVector.hpp>

#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ThyraUtils.hpp>

#include "MueLu_TekoSmoother.hpp"
#include "MueLu_MergedBlockedMatrixFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TekoSmoother()
    : type_("Teko smoother"), A_(Teuchos::null), bA_(Teuchos::null), bThyOp_(Teuchos::null), tekoParams_(Teuchos::null), inverseOp_(Teuchos::null)
  {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A");
    validParamList->set< std::string >           ("Inverse Type", "", "Name of parameter list within 'Teko parameters' containing the Teko smoother parameters.");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    FactoryMonitor m(*this, "Setup TekoSmoother", currentLevel);
    if (this->IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::TekoSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_  = Factory::Get< RCP<Matrix> >(currentLevel, "A"); // A needed for extracting map extractors
    bA_ = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: input matrix A is not of type BlockedCrsMatrix.");

    bThyOp_ = bA_->getThyraOperator();
    TEUCHOS_TEST_FOR_EXCEPTION(bThyOp_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: Could not extract thyra operator from BlockedCrsMatrix.");

    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > thyOp = Teuchos::rcp_dynamic_cast<const Thyra::LinearOpBase<Scalar> >(bThyOp_);
    TEUCHOS_TEST_FOR_EXCEPTION(thyOp.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: Downcast of Thyra::BlockedLinearOpBase to Teko::LinearOp failed.");

    // parameter list contains TekoSmoother parameters but does not handle the Teko parameters itself!
    const ParameterList& pL = Factory::GetParameterList();
    std::string smootherType = pL.get<std::string>("Inverse Type");
    TEUCHOS_TEST_FOR_EXCEPTION(smootherType.empty(), Exceptions::RuntimeError,
        "MueLu::TekoSmoother::Build: You must provide a 'Smoother Type' name that is defined in the 'Teko parameters' sublist.");
    type_ = smootherType;

    TEUCHOS_TEST_FOR_EXCEPTION(tekoParams_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: No Teko parameters have been set.");

    Teuchos::RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*tekoParams_);
    Teuchos::RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory(smootherType);

    inverseOp_ = Teko::buildInverse(*inverse, thyOp);
    TEUCHOS_TEST_FOR_EXCEPTION(inverseOp_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: Failed to build Teko inverse operator. Probably a problem with the Teko parameters.");

    this->IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(this->IsSetup() == false, Exceptions::RuntimeError,
                               "MueLu::TekoSmoother::Apply(): Setup() has not been called");

    Teuchos::RCP<const Teuchos::Comm<int> > comm = X.getMap()->getComm();

    Teuchos::RCP<const MapExtractor> rgMapExtractor = bA_->getRangeMapExtractor();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgMapExtractor));

    // copy initial solution vector X to Ptr<Thyra::MultiVectorBase> YY

    // create a Thyra RHS vector
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > thyB = Thyra::createMembers(Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(bThyOp_->productRange()),Teuchos::as<int>(B.getNumVectors()));
    Teuchos::RCP<Thyra::ProductMultiVectorBase<Scalar> > thyProdB =
        Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(thyB);
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdB.is_null(), Exceptions::BadCast,
        "MueLu::TekoSmoother::Apply: Failed to cast range space to product range space.");

    // copy RHS vector B to Thyra::MultiVectorBase thyProdB
    Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateThyra(Teuchos::rcpFromRef(B), rgMapExtractor, thyProdB);

    // create a Thyra SOL vector
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > thyX = Thyra::createMembers(Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(bThyOp_->productDomain()),Teuchos::as<int>(X.getNumVectors()));
    Teuchos::RCP<Thyra::ProductMultiVectorBase<Scalar> > thyProdX =
        Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(thyX);
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdX.is_null(), Exceptions::BadCast,
        "MueLu::TekoSmoother::Apply: Failed to cast domain space to product domain space.");

    // copy RHS vector X to Thyra::MultiVectorBase thyProdX
    Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::updateThyra(Teuchos::rcpFromRef(X), rgMapExtractor, thyProdX);

    inverseOp_->apply(
      Thyra::NOTRANS,
      *thyB, //const MultiVectorBase<Scalar> &X,
      thyX.ptr(), //const Ptr<MultiVectorBase<Scalar> > &Y,
      1.0,
      0.0);

    // copy back content of Ptr<Thyra::MultiVectorBase> thyX into X
    Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > XX =
          Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toXpetra(thyX, comm);

    X.update(Teuchos::ScalarTraits<Scalar>::one(), *XX, Teuchos::ScalarTraits<Scalar>::zero());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    return rcp (new TekoSmoother (*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void TekoSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }

} // namespace MueLu

#endif // HAVE_MUELU_TEKO

#endif /* MUELU_TEKOSMOOTHER_DEF_HPP_ */
