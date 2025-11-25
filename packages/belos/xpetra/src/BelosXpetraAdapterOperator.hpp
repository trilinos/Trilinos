// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_XPETRA_ADAPTER_OPERATOR_HPP
#define BELOS_XPETRA_ADAPTER_OPERATOR_HPP

//Note: using MACRO HAVE_XPETRA_ instead of HAVE_MUELU_ because this file will eventually be moved to Xpetra

#include "Xpetra_ConfigDefs.hpp"

#include <BelosOperatorT.hpp>

namespace Belos {
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;

  //
  //! @name MueLu Adapter Exceptions
  //@{

  /** \brief XpetraOpFailure is thrown when a return value from an MueLu
   * call on an Xpetra::Operator or MueLu::Hierarchy is non-zero.
   */
  class XpetraOpFailure : public BelosError {public:
    XpetraOpFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //@}

  //! @name Belos operator for Xpetra
  //@{

  /*! @class XpetraOp
   *
   * @brief Implementation of the Belos::XpetraOp. It derives from the Belos::OperatorT templated on
   *        the Xpetra::MultiVector and the Tpetra::MultiVector (if Teptra is enabled)
   *        Note, in contrast to Belos::MueLuOp this operator administrates an Xpetra::Matrix<> object
   *        and implements the effect of a vector applied to the stored matrix.
   *
   *        The Belos::OperatorT class is a generalization of the Belos::Operator<> class, which
   *        deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator<> interface does).
   *
   *        This is the general implementation for Tpetra only.
   */
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class XpetraOp :
    public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    , public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Default constructor
    XpetraOp (const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & Op) : Op_(Op) {}

    //! Destructor.
    virtual ~XpetraOp() {};
    //@}

    //! @name Operator application method
    //@{

    /*! \brief This routine takes the Xpetra::MultiVector \c x and applies the operator
      to it resulting in the Xpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::XpetraOp::Apply, transpose mode != NOTRANS not supported.");

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
      to it resulting in the Tpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::MueLuTpetraOp::Apply, transpose mode != NOTRANS not supported.");


      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & temp_x = const_cast<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &>(x);

      const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
      Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Op_->apply(tX,tY);
    }

    RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const { return Op_; }

  private:

    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op_;
  };

  //@}

} // namespace Belos

#endif // BELOS_XPETRA_ADAPTER_OPERATOR_HPP
