#ifndef BELOS_XPETRA_ADAPTER_OPERATOR_HPP
#define BELOS_XPETRA_ADAPTER_OPERATOR_HPP

//Note: using MACRO HAVE_XPETRA_ instead of HAVE_MUELU_ because this file will eventually be moved to Xpetra

#ifdef HAVE_XPETRA_EPETRAEXT
#include <BelosOperator.hpp>
#endif

#include <BelosOperatorT.hpp>

#include "Xpetra_ConfigDefs.hpp"

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

  // TODO: doc
  // TODO: should be it named XpetraOp (if Xpetra used by other packages) ?
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > 
  class MueLuOp : 
    public OperatorT<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#ifdef HAVE_XPETRA_TPETRA
    , public OperatorT<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
#endif
  {  

  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    MueLuOp(const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & Op) : Op_(Op) {}
    
    //! Destructor.
    virtual ~MueLuOp() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Cthulu::MultiVector \c x and applies the operator
      to it resulting in the Cthulu::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure, 
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported."); 

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

#ifdef HAVE_XPETRA_TPETRA
    // TO SKIP THE TRAIT IMPLEMENTATION OF CTHULU::MULTIVECTOR
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

  private:
  
    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op_;
  };

  template <> 
  class MueLuOp<double, int, int>   
    : 
    public OperatorT<Xpetra::MultiVector<double, int, int> >
#ifdef HAVE_XPETRA_TPETRA
    , public OperatorT<Tpetra::MultiVector<double, int, int> >
#endif
#ifdef HAVE_XPETRA_EPETRA
    , public OperatorT<Epetra_MultiVector>
#endif
  {  

    typedef double Scalar;
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;
    typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;

  public:

    MueLuOp(const RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & Op) : Op_(Op) {}
    
    virtual ~MueLuOp() {};

    void Apply ( const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure, 
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported."); 

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

#ifdef HAVE_XPETRA_TPETRA
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

#ifdef HAVE_XPETRA_EPETRAEXT
    // TO SKIP THE TRAIT IMPLEMENTATION OF CTHULU::MULTIVECTOR
    /*! \brief This routine takes the Epetra_MultiVector \c x and applies the operator
      to it resulting in the Epetra_MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans=NOTRANS ) const {
      TEUCHOS_TEST_FOR_EXCEPTION(trans!=NOTRANS, XpetraOpFailure, 
                         "Belos::MueLuTpetraOp::Apply, transpose mode != NOTRANS not supported."); 

      Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(x);

      const Xpetra::EpetraMultiVector tX(rcpFromRef(temp_x));
      Xpetra::EpetraMultiVector tY(rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      tY.putScalar(0.0);

      Op_->apply(tX,tY);
    }
#endif

  private:
  
    RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op_;
  };

} // namespace Belos

#endif // BELOS_XPETRA_ADAPTER_OPERATOR_HPP
