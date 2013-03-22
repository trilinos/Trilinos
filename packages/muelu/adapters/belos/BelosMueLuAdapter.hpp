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
#ifndef BELOS_MUELU_ADAPTER_HPP
#define BELOS_MUELU_ADAPTER_HPP

#ifdef HAVE_MUELU_EPETRA
#include <BelosOperator.hpp>
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
  class MueLuOpFailure : public BelosError {public:
    MueLuOpFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  //TODO: doc
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = Kokkos::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class MueLuOp :
    public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#ifdef HAVE_MUELU_TPETRA
    , public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#endif
  {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Default constructor
    MueLuOp(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & H) : Hierarchy_(H) {}

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
    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {

      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure,
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Hierarchy_->Iterate( x, 1, y , true);

    }

#ifdef HAVE_MUELU_TPETRA
    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
      to it resulting in the Tpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {

      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure,
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & temp_x = const_cast<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &>(x);

      const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
      Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Hierarchy_->Iterate( tX, 1, tY , true);

    }
#endif

  private:

    RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Hierarchy_;
  };

  template <>
  class MueLuOp<double, int, int> :
    public OperatorT<Xpetra::MultiVector<double, int, int> >
#ifdef HAVE_MUELU_TPETRA
    , public OperatorT<Tpetra::MultiVector<double, int, int> >
#endif
#ifdef HAVE_MUELU_EPETRA
    , public OperatorT<Epetra_MultiVector>
    , public Belos::Operator<double>
#endif
  {

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

  public:

    MueLuOp(const RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & H) : Hierarchy_(H) {}

    virtual ~MueLuOp() {}

    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {

      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure,
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Hierarchy_->Iterate( x, 1, y , true);

    }

#ifdef HAVE_MUELU_TPETRA
    void Apply ( const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {

      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure,
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

      Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & temp_x = const_cast<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &>(x);

      const Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tX(rcpFromRef(temp_x));
      Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Hierarchy_->Iterate( tX, 1, tY , true);

    }
#endif

#ifdef HAVE_MUELU_EPETRA
    // TO SKIP THE TRAIT IMPLEMENTATION OF XPETRA::MULTIVECTOR
    /*! \brief This routine takes the Tpetra::MultiVector \c x and applies the operator
      to it resulting in the Tpetra::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans=NOTRANS ) const {

      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure,
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners.");

      Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(x);

      const Xpetra::EpetraMultiVector tX(rcpFromRef(temp_x));
      Xpetra::EpetraMultiVector       tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Hierarchy_->Iterate( tX, 1, tY , true);

    }

    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Belos::MultiVec<double>& x, Belos::MultiVec<double>& y, ETrans trans=NOTRANS ) const {
      const Epetra_MultiVector* vec_x = dynamic_cast<const Epetra_MultiVector*>(&x);
      Epetra_MultiVector*       vec_y = dynamic_cast<Epetra_MultiVector*>(&y);

      TEUCHOS_TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, MueLuOpFailure,
                                  "Belos::MueLuOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");

      Apply(*vec_x, *vec_y, trans);
    }
#endif

  private:

    RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Hierarchy_;
  };

} // namespace Belos

#endif // BELOS_MUELU_ADAPTER_HPP
