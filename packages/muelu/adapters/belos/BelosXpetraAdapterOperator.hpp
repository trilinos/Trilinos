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
#ifndef BELOS_XPETRA_ADAPTER_OPERATOR_HPP
#define BELOS_XPETRA_ADAPTER_OPERATOR_HPP

//Note: using MACRO HAVE_XPETRA_ instead of HAVE_MUELU_ because this file will eventually be moved to Xpetra

#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Epetra_config.h>
#include <BelosOperator.hpp>
#endif

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
#ifdef HAVE_XPETRA_TPETRA
    , public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#endif
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

#ifdef HAVE_XPETRA_TPETRA
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
#endif

    RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const { return Op_; }

  private:

    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op_;
  };

#ifdef HAVE_XPETRA_EPETRA
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  /*! @class XpetraOp
   *
   * @brief Implementation of the Belos::XpetraOp. It derives from the Belos::OperatorT templated on
   *        the Xpetra::MultiVector and/or the Tpetra::MultiVector (if Teptra is enabled) and/or the
   *        EpetraMultiVector (if Epetra is enabled)
   *
   *        This is the specialization for <double,int,int,Xpetra::EpetraNode>
   */
  template <>
  class XpetraOp<double, int, int, Xpetra::EpetraNode>
    :
    public OperatorT<Xpetra::MultiVector<double, int, int, Xpetra::EpetraNode> >
#ifdef HAVE_XPETRA_TPETRA
#if !((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
    , public OperatorT<Tpetra::MultiVector<double, int, int, Xpetra::EpetraNode> >
#endif
#endif
#ifdef HAVE_XPETRA_EPETRA
    , public OperatorT<Epetra_MultiVector>
#endif
  {
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Xpetra::EpetraNode Node;

  public:

    XpetraOp(const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & Op) : Op_(Op) {}

    virtual ~XpetraOp() {};

    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::XpetraOp::Apply, transpose mode != NOTRANS not supported.");

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

#ifdef HAVE_XPETRA_TPETRA
#if !((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))
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
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Epetra_MultiVector \c x and applies the operator
      to it resulting in the Epetra_MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::MueLuTpetraOp::Apply, transpose mode != NOTRANS not supported.");

      Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(x);

      const Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node> tX(rcpFromRef(temp_x));
      Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node>       tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Op_->apply(tX,tY);
    }
#endif

    RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const { return Op_; }

  private:

    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op_;
  };
#endif // !EPETRA_NO_32BIT_GLOBAL_INDICES
#endif // HAVE_XPETRA_EPETRA

#ifdef HAVE_XPETRA_EPETRA
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  /*! @class XpetraOp
   *
   * @brief Implementation of the Belos::XpetraOp. It derives from the Belos::OperatorT templated on
   *        the Xpetra::MultiVector and/or the Tpetra::MultiVector (if Teptra is enabled) and/or the
   *        EpetraMultiVector (if Epetra64 is enabled). Please be aware that Epetra64 is not supported
   *        by MueLu. Don't expect it to work or produce reasonable results.
   *
   *        This is the specialization for <double,int,long long,Xpetra::EpetraNode>
   */
  template <>
  class XpetraOp<double, int, long long, Xpetra::EpetraNode>
    :
    public OperatorT<Xpetra::MultiVector<double, int, long long, Xpetra::EpetraNode> >
#ifdef HAVE_XPETRA_TPETRA
#if !((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
    , public OperatorT<Tpetra::MultiVector<double, int, long long, Xpetra::EpetraNode> >
#endif
#endif
#ifdef HAVE_XPETRA_EPETRA
    , public OperatorT<Epetra_MultiVector>
#endif
  {
    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef Xpetra::EpetraNode Node;

  public:

    XpetraOp(const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > & Op) : Op_(Op) {}

    virtual ~XpetraOp() {};

    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::XpetraOp::Apply, transpose mode != NOTRANS not supported.");

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

#ifdef HAVE_XPETRA_TPETRA
#if !((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
     (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))
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
#endif
#endif

#ifdef HAVE_XPETRA_EPETRA
    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Epetra_MultiVector \c x and applies the operator
      to it resulting in the Epetra_MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure,
                         "Belos::MueLuTpetraOp::Apply, transpose mode != NOTRANS not supported.");

      Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(x);

      const Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node> tX(rcpFromRef(temp_x));
      Xpetra::EpetraMultiVectorT<GlobalOrdinal,Node>       tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Op_->apply(tX,tY);
    }
#endif

    RCP<const Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const { return Op_; }

  private:

    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Op_;
  };
#endif // !EPETRA_NO_64BIT_GLOBAL_INDICES
#endif // HAVE_XPETRA_EPETRA

  //@}

} // namespace Belos

#endif // BELOS_XPETRA_ADAPTER_OPERATOR_HPP
