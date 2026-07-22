// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_MUELU_ADAPTER_HPP
#define BELOS_MUELU_ADAPTER_HPP

// TAW: 3/4/2016: we use the Xpetra macros
//                These are available and Xpetra is prerequisite for MueLu

#ifdef HAVE_MUELU_AMGX
#include "MueLu_AMGXOperator.hpp"
#endif

#include <BelosOperatorT.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Hierarchy.hpp"

namespace Belos {
using Teuchos::RCP;
using Teuchos::rcpFromRef;

//
//! @name MueLu Adapter Exceptions
//@{

/** \brief MueLuOpFailure is thrown when a return value from an MueLu
 * call on an Xpetra::Operator or MueLu::Hierarchy is non-zero.
 */
class MueLuOpFailure : public BelosError {
 public:
  MueLuOpFailure(const std::string& what_arg)
    : BelosError(what_arg) {}
};

/*! @class MueLuOp
 *
 * @brief MueLuOp derives from Belos::OperatorT and administrates a MueLu::Hierarchy. It implements the apply
 *        call which represents the effect of the multigrid preconditioner on a given vector.
 *        Note, in contrast to Belos::XpetraOp this operator has the multigrid hierarchy.
 *
 *        The Belos::OperatorT class is a generalization of the Belos::Operator<> class, which
 *        deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator<> interface does).
 *
 *        This is the general implementation for Tpetra only.
 */
template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MueLuOp : public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >,
                public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  //! @name Constructor/Destructor
  //@{

  //! Default constructor
  MueLuOp(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H) {}
#ifdef HAVE_MUELU_AMGX
  MueLuOp(const RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
    : AMGX_(A) {}
#endif
  //! Destructor.
  virtual ~MueLuOp() {}
  //@}

  //! @name Operator application method
  //@{

  /*! \brief This routine takes the Xpetra::MultiVector \c x and applies the operator
    to it resulting in the Xpetra::MultiVector \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    // This does not matter for Hierarchy, but matters for AMGX
    y.putScalar(0.0);

#ifdef HAVE_MUELU_AMGX
    if (!AMGX_.is_null()) {
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX = Xpetra::toTpetra(x);
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY = Xpetra::toTpetra(y);

      AMGX_->apply(tX, tY);
    }
#endif
    if (!Hierarchy_.is_null())
      Hierarchy_->Iterate(x, y, 1, true);
  }
  //@}

  // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
  /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
    to it resulting in the Tpetra::MultiVector \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    // FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate(), but it matters for AMGX
    y.putScalar(0.0);

#ifdef HAVE_MUELU_AMGX
    if (!AMGX_.is_null())
      AMGX_->apply(x, y);
#endif

    if (!Hierarchy_.is_null()) {
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& temp_x = const_cast<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&>(x);

      const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
      Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY(rcpFromRef(y));
      Hierarchy_->Iterate(tX, tY, 1, true);
    }
  }

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
#ifdef HAVE_MUELU_AMGX
  RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AMGX_;
#endif
};

}  // namespace Belos

#endif  // BELOS_MUELU_ADAPTER_HPP
