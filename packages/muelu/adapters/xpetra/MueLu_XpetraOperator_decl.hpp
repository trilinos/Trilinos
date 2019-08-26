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
#ifndef MUELU_XPETRAOPERATOR_DECL_HPP
#define MUELU_XPETRAOPERATOR_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"

namespace MueLu {

/*!  @brief Wraps an existing MueLu::Hierarchy as a Xpetra::Operator.
*/
  template <class Scalar = DefaultScalar,
            class LocalOrdinal = DefaultLocalOrdinal,
            class GlobalOrdinal = DefaultGlobalOrdinal,
            class Node = DefaultNode>
  class XpetraOperator : public Xpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  protected:
    XpetraOperator() { }
  public:

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    XpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H) : Hierarchy_(H) { }

    //! Destructor.
    virtual ~XpetraOperator() { }

    //@}

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const {
      typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

      RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
      return A->getDomainMap();
    }

    //! Returns the Tpetra::Map object associated with the range of this operator.
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const {
      typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;

      RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");
      return A->getRangeMap();
    }

    //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
    /*!
      \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(const Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                         Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                         Teuchos::ETransp /* mode */ = Teuchos::NO_TRANS,
                                         Scalar /* alpha */ = Teuchos::ScalarTraits<Scalar>::one(),
                                         Scalar /* beta */  = Teuchos::ScalarTraits<Scalar>::one()) const{
      try {
#ifdef HAVE_MUELU_DEBUG
        typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
        RCP<Matrix> A = Hierarchy_->GetLevel(0)->template Get<RCP<Matrix> >("A");

        // X is supposed to live in the range map of the operator (const rhs = B)
        RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xop =
            Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A->getRangeMap(),X.getNumVectors());
        RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Yop =
            Xpetra::MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A->getDomainMap(),Y.getNumVectors());
        TEUCHOS_TEST_FOR_EXCEPTION(A->getRangeMap()->isSameAs(*(Xop->getMap())) == false, std::logic_error,
                                   "MueLu::XpetraOperator::apply: map of X is incompatible with range map of A");
        TEUCHOS_TEST_FOR_EXCEPTION(A->getDomainMap()->isSameAs(*(Yop->getMap())) == false, std::logic_error,
                                   "MueLu::XpetraOperator::apply: map of Y is incompatible with domain map of A");
#endif

        Y.putScalar(Teuchos::ScalarTraits<Scalar>::zero());
        Hierarchy_->Iterate(X, Y, 1, true);
      } catch (std::exception& e) {
        //FIXME add message and rethrow
        std::cerr << "Caught an exception in MueLu::XpetraOperator::apply():" << std::endl
            << e.what() << std::endl;
      }
    }

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const { return false; }

#ifdef HAVE_MUELU_DEPRECATED_CODE
    template <class NewNode>
    Teuchos::RCP< XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, NewNode> >
    MUELU_DEPRECATED
    clone(const RCP<NewNode>& new_node) const {
      return Teuchos::rcp (new XpetraOperator<Scalar, LocalOrdinal, GlobalOrdinal, NewNode> (Hierarchy_->template clone<NewNode> (new_node)));
    }
#endif

    //! @name MueLu specific
    //@{

    //! Direct access to the underlying MueLu::Hierarchy.
    RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GetHierarchy() const { return Hierarchy_; }

    //@}

  private:
    RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
  };

} // namespace

#endif // MUELU_XPETRAOPERATOR_DECL_HPP
