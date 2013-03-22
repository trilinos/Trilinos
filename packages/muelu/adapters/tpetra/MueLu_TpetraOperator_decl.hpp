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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_TPETRAOPERATOR_DECL_HPP
#define MUELU_TPETRAOPERATOR_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_TPETRA
#include <Tpetra_Operator.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"
//TODO: Kokkos headers

/*! @class TpetraOperator
    Wraps an existing MueLu::Hierarchy as a Tpetra::Operator.
*/

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal,
            class Node = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps = typename Kokkos::DefaultKernels<Scalar, LocalOrdinal, Node>::SparseOps >
  class TpetraOperator
    : public Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Constructor
    TpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & H) : Hierarchy_(H) { }

    //! Destructor.
    virtual ~TpetraOperator() { }

    //@}

    //! Returns the Tpetra::Map object associated with the domain of this operator.
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

    //! Returns the Tpetra::Map object associated with the range of this operator.
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

    //! Returns in Y the result of a Tpetra::Operator applied to a Tpetra::MultiVector X.
    /*!
      \param[in] X - Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
      \param[out] Y -Tpetra::MultiVector of dimension NumVectors containing result.

    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                         Teuchos::ETransp mode = Teuchos::NO_TRANS,
                                         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                                         Scalar beta = Teuchos::ScalarTraits<Scalar>::one()) const;

    //! Indicates whether this operator supports applying the adjoint operator.
    bool hasTransposeApply() const;

  private:

    RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Hierarchy_;

  };

} // namespace

#endif //ifdef HAVE_MUELU_TPETRA

#endif // MUELU_TPETRAOPERATOR_DECL_HPP
