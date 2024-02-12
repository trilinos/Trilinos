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

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy_decl.hpp"

namespace MueLu {

/*!  @brief Wraps an existing MueLu::Hierarchy as a Xpetra::Operator.
 */
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class XpetraOperator : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 protected:
  XpetraOperator() {}

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  XpetraOperator(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H) {}

  //! Destructor.
  virtual ~XpetraOperator() {}

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > getRangeMap() const;

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  /*!
    \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar /* alpha */    = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar /* beta */     = Teuchos::ScalarTraits<Scalar>::one()) const;

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const { return false; }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const;

  //! @name MueLu specific
  //@{

  //! Direct access to the underlying MueLu::Hierarchy.
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > GetHierarchy() const { return Hierarchy_; }

  //@}

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
};

}  // namespace MueLu

#endif  // MUELU_XPETRAOPERATOR_DECL_HPP
