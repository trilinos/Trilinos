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

#ifndef MUELU_TEKOSMOOTHER_DECL_HPP_
#define MUELU_TEKOSMOOTHER_DECL_HPP_

#ifdef HAVE_MUELU_TEKO

#include "Teko_Utilities.hpp"

#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_MapExtractor_fwd.hpp>

#include "MueLu_TekoSmoother_fwd.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

/*!
  @class TekoSmoother
  @brief Interface to block smoothers in Teko package

  This is a dummy implementation for arbitrary scalars (different than SC=double).
  We use Teko and the Teko::LinearOp which is declared as
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > with ST=double. See e.g. Teko_ConfigDefs.hpp.
*/

template <class Scalar        = SmootherPrototype<>::scalar_type,
          class LocalOrdinal  = typename SmootherPrototype<Scalar>::local_ordinal_type,
          class GlobalOrdinal = typename SmootherPrototype<Scalar, LocalOrdinal>::global_ordinal_type,
          class Node          = typename SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal>::node_type>
class TekoSmoother : public SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_TEKOSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor
   */
  TekoSmoother()
    : type_("Teko smoother") {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::TekoSmoother: Teko can only be used with SC=double. For more information refer to the doxygen documentation of TekoSmoother.");
  };

  //! Destructor
  virtual ~TekoSmoother() {}
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    return validParamList;
  }

  void DeclareInput(Level &currentLevel) const {}

  void SetTekoParameters(RCP<ParameterList> tekoParams){};
  //@}

  //! @name Setup and Apply methods.
  //@{

  /*! @brief Setup routine
   */
  void Setup(Level &currentLevel) {}

  /*! @brief Apply the Teko smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero This option has no effect.
  */
  void Apply(MultiVector &X, const MultiVector &B, bool InitialGuessIsZero = false) const {}
  //@}

  RCP<SmootherPrototype> Copy() const { return Teuchos::null; }

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const;

  //@}

 private:
  //! smoother type
  std::string type_;
};  // class TekoSmoother

/*!
  @class TekoSmoother
  @brief Interface to block smoothers in Teko package

  This is the specialization for SC=double and LO=int, which actually implements the Teko smoother interface.
  We use Teko and the Teko::LinearOp which is declared as
  Teuchos::RCP<const Thyra::LinearOpBase<ST> > with ST=double. See e.g. Teko_ConfigDefs.hpp.

  In Teko SC=double and LO=int are hard-coded. The node is Tpetra::Map<>::node_type.
  The global ordinal GO is chosen as a conservative choice from what Tpetra supports
  (see Teko_ConfigDefs.hpp for more details)

  Note, that the underlying Thyra package is only templated on the Scalar (not the node type).
*/
template <class GlobalOrdinal,
          class Node>
class TekoSmoother<double, int, GlobalOrdinal, Node> : public SmootherPrototype<double, int, GlobalOrdinal, Node> {
  typedef int LocalOrdinal;
  typedef double Scalar;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;

#undef MUELU_TEKOSMOOTHER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors / destructors
  //@{

  /*! @brief Constructor
   */
  TekoSmoother()
    : type_("Teko smoother")
    , A_(Teuchos::null)
    , bA_(Teuchos::null)
    , bThyOp_(Teuchos::null)
    , tekoParams_(Teuchos::null)
    , inverseOp_(Teuchos::null){};

  //! Destructor
  virtual ~TekoSmoother() {}
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set<RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A");
    validParamList->set<std::string>("Inverse Type", "", "Name of parameter list within 'Teko parameters' containing the Teko smoother parameters.");

    return validParamList;
  }

  void DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  void SetTekoParameters(RCP<ParameterList> tekoParams) { tekoParams_ = tekoParams; };
  //@}

  //! @name Setup and Apply methods.
  //@{

  /*! @brief Setup routine
   */
  void Setup(Level &currentLevel) {
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    FactoryMonitor m(*this, "Setup TekoSmoother", currentLevel);
    if (this->IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::TekoSmoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_  = Factory::Get<RCP<Matrix> >(currentLevel, "A");  // A needed for extracting map extractors
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
    const ParameterList &pL  = Factory::GetParameterList();
    std::string smootherType = pL.get<std::string>("Inverse Type");
    TEUCHOS_TEST_FOR_EXCEPTION(smootherType.empty(), Exceptions::RuntimeError,
                               "MueLu::TekoSmoother::Build: You must provide a 'Smoother Type' name that is defined in the 'Teko parameters' sublist.");
    type_ = smootherType;

    TEUCHOS_TEST_FOR_EXCEPTION(tekoParams_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: No Teko parameters have been set.");

    Teuchos::RCP<Teko::InverseLibrary> invLib  = Teko::InverseLibrary::buildFromParameterList(*tekoParams_);
    Teuchos::RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory(smootherType);

    inverseOp_ = Teko::buildInverse(*inverse, thyOp);
    TEUCHOS_TEST_FOR_EXCEPTION(inverseOp_.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Build: Failed to build Teko inverse operator. Probably a problem with the Teko parameters.");

    this->IsSetup(true);
  }

  /*! @brief Apply the Teko smoother.

  @param X initial guess
  @param B right-hand side
  @param InitialGuessIsZero This option has no effect.
  */
  void Apply(MultiVector &X, const MultiVector &B, bool /* InitialGuessIsZero */ = false) const {
    TEUCHOS_TEST_FOR_EXCEPTION(this->IsSetup() == false, Exceptions::RuntimeError,
                               "MueLu::TekoSmoother::Apply(): Setup() has not been called");

    Teuchos::RCP<const Teuchos::Comm<int> > comm = X.getMap()->getComm();

    Teuchos::RCP<const MapExtractor> rgMapExtractor = bA_->getRangeMapExtractor();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(rgMapExtractor));

    // copy initial solution vector X to Ptr<Thyra::MultiVectorBase> YY

    // create a Thyra RHS vector
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > thyB = Thyra::createMembers(Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(bThyOp_->productRange()), Teuchos::as<int>(B.getNumVectors()));
    Teuchos::RCP<Thyra::ProductMultiVectorBase<Scalar> > thyProdB =
        Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(thyB);
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdB.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Apply: Failed to cast range space to product range space.");

    // copy RHS vector B to Thyra::MultiVectorBase thyProdB
    Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::updateThyra(Teuchos::rcpFromRef(B), rgMapExtractor, thyProdB);

    // create a Thyra SOL vector
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > thyX = Thyra::createMembers(Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<Scalar> >(bThyOp_->productDomain()), Teuchos::as<int>(X.getNumVectors()));
    Teuchos::RCP<Thyra::ProductMultiVectorBase<Scalar> > thyProdX =
        Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<Scalar> >(thyX);
    TEUCHOS_TEST_FOR_EXCEPTION(thyProdX.is_null(), Exceptions::BadCast,
                               "MueLu::TekoSmoother::Apply: Failed to cast domain space to product domain space.");

    // copy RHS vector X to Thyra::MultiVectorBase thyProdX
    Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::updateThyra(Teuchos::rcpFromRef(X), rgMapExtractor, thyProdX);

    inverseOp_->apply(
        Thyra::NOTRANS,
        *thyB,       // const MultiVectorBase<Scalar> &X,
        thyX.ptr(),  // const Ptr<MultiVectorBase<Scalar> > &Y,
        1.0,
        0.0);

    // copy back content of Ptr<Thyra::MultiVectorBase> thyX into X
    Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > XX =
        Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(thyX, comm);

    X.update(Teuchos::ScalarTraits<Scalar>::one(), *XX, Teuchos::ScalarTraits<Scalar>::zero());
  }
  //@}

  RCP<SmootherPrototype> Copy() const { return Teuchos::rcp(new MueLu::TekoSmoother<double, int, GlobalOrdinal, Node>(*this)); }

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  //! Print the object with some verbosity level to an FancyOStream object.
  // using MueLu::Describable::describe; // overloading, not hiding
  void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }

  //! Get a rough estimate of cost per iteration
  size_t getNodeSmootherComplexity() const {
    size_t cplx = 0;
    return cplx;
  }

  //@}

 private:
  //! smoother type
  std::string type_;

  //! block operator
  RCP<Matrix> A_;  // < ! internal blocked operator "A" generated by AFact_
  RCP<BlockedCrsMatrix> bA_;
  RCP<const Thyra::BlockedLinearOpBase<Scalar> > bThyOp_;

  //! Teko parameters
  RCP<ParameterList> tekoParams_;  // < ! parameter list containing Teko parameters. These parameters are not administrated by the factory and not validated.

  Teko::LinearOp inverseOp_;  // < ! Teko inverse operator
};                            // class TekoSmoother (specialization on SC=double)
}  // namespace MueLu

#define MUELU_TEKOSMOOTHER_SHORT

#endif  // HAVE_MUELU_TEKO

#endif /* MUELU_TEKOSMOOTHER_DECL_HPP_ */
