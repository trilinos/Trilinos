// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSXPETRADIAGOPERATOR_HPP
#define BELOSXPETRADIAGOPERATOR_HPP


//Note: using MACRO HAVE_XPETRA_ instead of HAVE_MUELU_ because this file will eventually be moved to Xpetra

#include "Xpetra_ConfigDefs.hpp"

#include <BelosOperatorT.hpp>

namespace Belos {
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;

  //! @name Belos diagonal operator for Xpetra
  //@{

  /*! @class XpetraDiagOp
   *
   * @brief Implementation of the Belos::XpetraDiagOp. It derives from the Belos::OperatorT templated on
   *        the Xpetra::MultiVector and the Tpetra::MultiVector (if Tpetra is enabled)
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
  class XpetraDiagOp :
    public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    , public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Default constructor
    XpetraDiagOp (const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & diag) : diag_(diag) {}

    //! Destructor.
    virtual ~XpetraDiagOp() {};
    //@}

    //! @name Operator application method
    //@{

    /*! \brief This routine takes the Xpetra::MultiVector \c x and applies the operator
      to it resulting in the Xpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      y.elementWiseMultiply(Teuchos::ScalarTraits<Scalar>::one(), *diag_, x, Teuchos::ScalarTraits<Scalar>::zero());
    }

    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
      to it resulting in the Tpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & temp_x = const_cast<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &>(x);

      const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
      Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY(rcpFromRef(y));

      tY.elementWiseMultiply(Teuchos::ScalarTraits<Scalar>::one(), *diag_, tX, Teuchos::ScalarTraits<Scalar>::zero());
    }

    RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getDiag() const { return diag_; }

  private:

    RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node >> diag_;
  };

  //@}

} // namespace Belos

#endif
