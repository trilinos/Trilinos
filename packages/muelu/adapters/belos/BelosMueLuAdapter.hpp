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
#ifndef BELOS_MUELU_ADAPTER_HPP
#define BELOS_MUELU_ADAPTER_HPP

// TAW: 3/4/2016: we use the Xpetra macros
//                These are available and Xpetra is prerequisite for MueLu
#ifdef HAVE_XPETRA_EPETRA
#include <Epetra_config.h>
#include <BelosOperator.hpp>
#endif

//#ifdef HAVE_XPETRA_TPETRA
//#include "TpetraCore_config.h"
//#endif

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
class MueLuOp : public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#ifdef HAVE_XPETRA_TPETRA
  ,
                public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#endif
{
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

#ifdef HAVE_XPETRA_TPETRA
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
#endif

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
#ifdef HAVE_MUELU_AMGX
  RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AMGX_;
#endif
};

#ifdef HAVE_XPETRA_EPETRA
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
/*! @class MueLuOp
 *
 * @brief MueLuOp derives from Belos::OperatorT and administrates a MueLu::Hierarchy. It implements the apply
 *        call which represents the effect of the multigrid preconditioner on a given vector.
 *        Note, in contrast to Belos::XpetraOp this operator has the multigrid hierarchy.
 *
 *        The Belos::OperatorT class is a generalization of the Belos::Operator<> class, which
 *        deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator<> interface does).
 *
 *        This is the specialization for <double,int,int,Xpetra::EpetraNode>
 */
template <>
class MueLuOp<double, int, int, Xpetra::EpetraNode> : public OperatorT<Xpetra::MultiVector<double, int, int, Xpetra::EpetraNode> >
#ifdef HAVE_XPETRA_TPETRA
// check whether Tpetra is instantiated on double,int,int,EpetraNode
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
  ,
                                                      public OperatorT<Tpetra::MultiVector<double, int, int, Xpetra::EpetraNode> >
#endif
#endif
#ifdef HAVE_XPETRA_EPETRA
  ,
                                                      public OperatorT<Epetra_MultiVector>,
                                                      public Belos::Operator<double>
#endif
{
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef Xpetra::EpetraNode Node;

 public:
  MueLuOp(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H) {}
#ifdef HAVE_MUELU_AMGX
  MueLuOp(const RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
    : AMGX_(A) {}
#endif
  virtual ~MueLuOp() {}

  void Apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    // FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate(), but it matters for AMGX
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

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT))))
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

      tY.putScalar(0.0);

      Hierarchy_->Iterate(tX, tY, 1, true);
    }
  }
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
  // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
  /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
    to it resulting in the Tpetra::MultiVector \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    Epetra_MultiVector& temp_x = const_cast<Epetra_MultiVector&>(x);

    const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
    Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> tY(rcpFromRef(y));

    // FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
    tY.putScalar(0.0);

    Hierarchy_->Iterate(tX, tY, 1, true);
  }

  /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
    to it resulting in the Belos::MultiVec \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Belos::MultiVec<double>& x, Belos::MultiVec<double>& y, ETrans trans = NOTRANS) const {
    const Epetra_MultiVector* vec_x = dynamic_cast<const Epetra_MultiVector*>(&x);
    Epetra_MultiVector* vec_y       = dynamic_cast<Epetra_MultiVector*>(&y);

    TEUCHOS_TEST_FOR_EXCEPTION(vec_x == NULL || vec_y == NULL, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");

    Apply(*vec_x, *vec_y, trans);
  }
#endif

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
#ifdef HAVE_MUELU_AMGX
  RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AMGX_;
#endif
};
#endif  // !EPETRA_NO_32BIT_GLOBAL_INDICES
#endif  // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_EPETRA
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
/*! @class MueLuOp
 *
 * @brief MueLuOp derives from Belos::OperatorT and administrates a MueLu::Hierarchy. It implements the apply
 *        call which represents the effect of the multigrid preconditioner on a given vector.
 *        Note, in contrast to Belos::XpetraOp this operator has the multigrid hierarchy.
 *
 *        The Belos::OperatorT class is a generalization of the Belos::Operator<> class, which
 *        deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator<> interface does).
 *
 *        This is the specialization for <double,int,long long,Xpetra::EpetraNode>
 */
template <>
class MueLuOp<double, int, long long, Xpetra::EpetraNode> : public OperatorT<Xpetra::MultiVector<double, int, long long, Xpetra::EpetraNode> >
#ifdef HAVE_XPETRA_TPETRA
// check whether Tpetra is instantiated on double,int,int,EpetraNode
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
  ,
                                                            public OperatorT<Tpetra::MultiVector<double, int, long long, Xpetra::EpetraNode> >
#endif
#endif
#ifdef HAVE_XPETRA_EPETRA
  ,
                                                            public OperatorT<Epetra_MultiVector>,
                                                            public Belos::Operator<double>
#endif
{
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef Xpetra::EpetraNode Node;

 public:
  MueLuOp(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& H)
    : Hierarchy_(H) {}
#ifdef HAVE_MUELU_AMGX
  MueLuOp(const RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A)
    : AMGX_(A) {}
#endif
  virtual ~MueLuOp() {}

  void Apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    // FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate(), but it matters for AMGX
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

#ifdef HAVE_XPETRA_TPETRA
#if ((defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
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

      tY.putScalar(0.0);

      Hierarchy_->Iterate(tX, tY, 1, true);
    }
  }
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
  // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
  /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
    to it resulting in the Tpetra::MultiVector \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans = NOTRANS) const {
    TEUCHOS_TEST_FOR_EXCEPTION(trans != NOTRANS, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

    Epetra_MultiVector& temp_x = const_cast<Epetra_MultiVector&>(x);

    const Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
    Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node> tY(rcpFromRef(y));

    // FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
    tY.putScalar(0.0);

    Hierarchy_->Iterate(tX, tY, 1, true);
  }

  /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
    to it resulting in the Belos::MultiVec \c y, which is returned.
    \note It is expected that any problem with applying this operator to \c x will be
    indicated by an std::exception being thrown.
  */
  void Apply(const Belos::MultiVec<double>& x, Belos::MultiVec<double>& y, ETrans trans = NOTRANS) const {
    const Epetra_MultiVector* vec_x = dynamic_cast<const Epetra_MultiVector*>(&x);
    Epetra_MultiVector* vec_y       = dynamic_cast<Epetra_MultiVector*>(&y);

    TEUCHOS_TEST_FOR_EXCEPTION(vec_x == NULL || vec_y == NULL, MueLuOpFailure,
                               "Belos::MueLuOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");

    Apply(*vec_x, *vec_y, trans);
  }
#endif

 private:
  RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Hierarchy_;
#ifdef HAVE_MUELU_AMGX
  RCP<MueLu::AMGXOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > AMGX_;
#endif
};
#endif  // !EPETRA_NO_64BIT_GLOBAL_INDICES
#endif  // HAVE_XPETRA_EPETRA
}  // namespace Belos

#endif  // BELOS_MUELU_ADAPTER_HPP
