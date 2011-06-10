#ifndef BELOS_MUELU_ADAPTER_HPP
#define BELOS_MUELU_ADAPTER_HPP

#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"

#ifdef HAVE_CTHULHU_EPETRA
#include "BelosEpetraAdapter.hpp"
#endif

#include "MueLu_Hierarchy.hpp"

#include "BelosMueLuAdapterMultiVector.hpp" // this defines the MultiVecTraits for Cthulhu::MultiVector

// Here is some documentation about the adapter interface of Belos (not MueLu specific, 06/2011):
// Belos uses the traits techniques for its adapters. Traits OP and MV must be implemented for your operator and multivector classes.
// What is somehow confusing is that Belos also provides an interface Belos::Operator and Belos::MultiVec.
// Internally, Belos only use the traits, not the interface. But traits specialization for Belos::Operator and Belos::MultiVec are provided, so you can either:
// - implements directly the traits Belos::OperatorTraits and Belos::MultiVecTraits
// - implements the interface Belos::Operator and Belos::MultiVec
//
// Here are the adapaters provided by Belos:
// Traits are implemented for three couples of <MV,OP>, and Belos works out of the box with the following combinaison of <MV,OP>:
// - MV=Belos::MultiVec<...>, OP=Belos::Operator<...>
// - MV=Epetra_MultiVector, OP=Epetra_Operator
// - MV=Tpetra_MultiVector<...>, MV=Tpetra_MultiVector<...>
//
// In addition, wrappers around Epetra_MultiVector and Epetra_Operator are provided to turn these Epetra objects to Belos::MultiVec and Belos::Operator and use MV=Belos::MultiVec<...>, OP=Belos::Operator<...> with Epetra.
// The wrappers are the classes Belos::EpetraMultiVec and Belos::EpetraOp 
// So when using Epetra, you have the choice to:
// - use Belos::LinearProblem<double, Belos::Operator<double>, Belos::MultiVec<double>>
// - or Belos::LinearProblem<double, Epetra_Operator<double>, Epetra_MultiVec<double>>
//
// If you use Epetra, you have to be carreful with the meaning for Epetra_Operator::Apply:
// For instance, Ifpack smoothers implement the Epetra_Operator interface but to apply the preconditionner, you have to use ApplyInverse() instead of Apply()!
// To swap the method Apply() and ApplyInverse() of an Epetra_Operator, you can use the class Belos::EpetraPrecOp. This class can be used with both OP=Belos::Operator<...> and OP=Epetra_Operator.
//
// Belos files:
// - src/BelosMultiVecTraits.hpp src/BelosOperatorTraits.hpp : Traits used internally by Belos
// - tpetra/src/BelosTpetraAdapter.*                         : Specialization of the Traits for Tpetra
// - epetra/src/BelosEpetraAdapter.*                         : Specialization of the Traits for Epetra + Implementation of the interface Belos::MultiVec and Belos::Operator for Epetra + Implementation of Belos::EpetraPrecOp
// - src/BelosMultiVec.hpp src/BelosOperator.hpp             : Definition of the interfaces Belos::MultiVec and Belos::Operator + Specialization of the Traits for these interfaces
// - epetra/src/BelosEpetraOperator.*                        : Epetra adapter that wrap Belos objects into an Epetra_Operator (do not have anything to do with the current discussion)

//

//

namespace Belos { 
  // TODO: Should this file be moved to Belos ?
  // TODO: The relation between Belos and MueLu is: Belos uses MueLu as a Preconditionner. So it makes more sense to me.

  // Here are a list of the Belos adapters for MueLu. To use Belos::LinearProblem<ScalarType,MV,OP> with:
  // A - MV=Belos::MultiVec<ScalarType> and OP=Belos::Operator<ScalarType>, turns your MueLu::Hierarchy into a Belos::MueLuEpetraPrecOp
  // B - MV=Epetra_MultiVector          and OP=Epetra_Operator            , turns your MueLu::Hierarchy into a Belos::MueLuEpetraPrecOp (TODO: not available yet)
  // C - MV=Tpetra::MultiVector<...>    and OP=Tpetra_Operator<...>       , turns your MueLu::Hierarchy into a Belos::MueLuTpetraPrecOp (TODO: not available yet)
  // D - MV=Cthulhu::MultiVector<...>   and OP=Cthulhu::Operator<...>     , turns your MueLu::Hierarchy into a Belos::MueLuCthulhuPrecOp => TODO: description not exactly right
  // TODO: I can also quickly implements couples Tpetra::MultiVector/Cthulhu::Operator and Epetra_MultiVector/Cthulhu::Operator

  // -----------------------------------------------------------------------------------------------------------------------------------
  //  A: MV=Belos::MultiVec<ScalarType> and OP=Belos::Operator<ScalarType>
  // -----------------------------------------------------------------------------------------------------------------------------------

  // Turns a MueLu::Hierarchy<ScalarType,...> object to a Belos::Operator<ScalarType>.
  // This allows to use MueLu as a preconditionner for a Belos::LinearProblem
  // with ScalarType=, MV=Belos::MultiVec<ScalarType> and OP=Belos::Operator<ScalarType>
  //
  // Note: this adapter is implemented only for Epetra (and ScalarType=double), because the interface Belos::Operator and Belos::MultiVec is only implemented for Epetra in Belos.
  // For Tpetra, you can use directly the adapter provided for Belos::LinearProblem where OP=Tpetra::Operator<...> or OP=Cthulhu::Operator<...>
  //
  // TODO: This adapter may also be used for OP=Epetra_Operator if we implements the interface Epetra_Operator by adding inheritence to public virtual Epetra_Operator and a bunch of methods
  // (see also the class Belos::EpetraPrecOp in package/belos/epetra/src/ for guidance)
  class MueLuEpetraPrecOp : public Belos::Operator<double> { 
    
    typedef MueLu::Hierarchy<double,int, int, Kokkos::DefaultNode::DefaultNodeType, Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> Hierarchy; // TODO: remove template parameters of Hierarchy
    typedef Cthulhu::MultiVector<double, int, int> MultiVector;

  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    MueLuEpetraPrecOp(const Teuchos::RCP<Hierarchy> & H) : Hierarchy_(H) {}
    
    //! Destructor.
    virtual ~MueLuEpetraPrecOp() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Epetra_MultiVector \c x and applies the operator
      to it resulting in the Epetra_MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    // Note: throw EpetraOpFailure exceptions as Belos::EpetraOp
    void Apply ( const Epetra_MultiVector& x, Epetra_MultiVector& y, ETrans trans=NOTRANS ) const {

      TEST_FOR_EXCEPTION(trans!=NOTRANS, EpetraOpFailure, 
                         "Belos::MueLuEpetraPrecOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners."); 

      Epetra_MultiVector & temp_x = const_cast<Epetra_MultiVector &>(x);

      const Cthulhu::EpetraMultiVector eX(Teuchos::rcpFromRef(temp_x));
      Cthulhu::EpetraMultiVector eY(Teuchos::rcpFromRef(y));

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      eY.putScalar(0.0);

      Hierarchy_->Iterate( eX, 1, eY , true);
      
    }

    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    void Apply ( const Belos::MultiVec<double>& x, Belos::MultiVec<double>& y, ETrans trans=NOTRANS ) const {
      const Epetra_MultiVector* vec_x = dynamic_cast<const Epetra_MultiVector*>(&x);
      Epetra_MultiVector*       vec_y = dynamic_cast<Epetra_MultiVector*>(&y);

      TEST_FOR_EXCEPTION( vec_x==NULL || vec_y==NULL, EpetraOpFailure, 
                          "Belos::MueLuEpetraPrecOp::Apply, x and/or y cannot be dynamic cast to an Epetra_MultiVector.");

      Apply(*vec_x, *vec_y, trans);
    }
  
  private:
  
    Teuchos::RCP<Hierarchy> Hierarchy_;
  };

  // -----------------------------------------------------------------------------------------------------------------------------------
  //  D: MV=Cthulhu::MultiVector<...>   and OP=Cthulhu::Operator<...>
  // -----------------------------------------------------------------------------------------------------------------------------------

  // JG Notes about the implementation of Belos adapters for Cthulhu objects:
  // To use Belos with Cthulhu, we need here a common ancestor between classes Cthulhu::Operator and MueLu::Hierarchy
  // and then implement the Belos::OperatorTraits for this new class hierarchy ancestor.
  //
  // There is several way to do that:
  // 1) Cthulhu::Operator can be the common ancestor:
  //  - 1a) MueLu::Hierarchy implements directly Cthulhu::Operator
  //  - 1b) MueLu::Hierarchy is wrapped to an object that implements this interface
  // 2) Creates a new common interface and:
  // - 2a) Both MueLu::Hierarchy and Cthulhu::Operator inherit from it directly.
  // - 2b) Wrap both MueLu::Hierarchy and Cthulhu::Operator to respect the new interface.
  //
  // PB of 1): Right now, Cthulhu::Operator is way to complicated to be the common interface.
  //           At some point, Cthulhu::Operator should correspond to the Tpetra::Operator interface and the old one should be renamed
  //           1a) is the approach of Ifpack in some sense: Ifpack Preconditionner implements Epetra_Operator but there is problem 
  //           with Apply vs ApplyInverse and so it ends up with an approach more similar to 1b).
  //
  // PB of 2b): If the new interface is only for Belos, we should not introduce inheritence dependency to it
  //            Approach 2b) is very close to 1a) if the Cthulhu::Operator is changed to be === to Tpetra::Operator
  //
  // Righ now, I adopt the approach 2b). The common interface is OperatorT. 
  // This approach is very similar with what is done with Belos::Operator but here, we don't need to wrap Cthulhu::MultiVector
  // I think that the OperatorT interface should be provided by Belos and replace the Belos::Operator

  // TODO: remplace OperatorT<Scalar,LocalOrdinal...> by OperatorT<MV> for more genericity if this class moves to Belos (it can replace Belos::Operator)

  //
  // Base class for the Belos OP template type (as Belos::Operator<>) but this one deals with Cthulhu::MultiVector
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > 
  class OperatorT {
 
    typedef Cthulhu::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
    
  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    OperatorT() {};
    
    //! Destructor.
    virtual ~OperatorT() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    virtual void Apply ( const MultiVector & x, MultiVector & y, ETrans trans=NOTRANS ) const = 0;
  };
  
  /// \brief Specialization of OperatorTraits for OperatorT.
  ///
  /// This is a partial template specialization of
  /// Belos::OperatorTraits class using the Belos::OperatorT 
  /// abstract interface. Any class that inherits
  /// from Belos::OperatorT will be accepted by the Belos templated
  /// solvers, due to this specialization of Belos::OperatorTraits.
  template <class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class OperatorTraits<ScalarType, Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node>, OperatorT<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > 
  {

    typedef Cthulhu::MultiVector<ScalarType, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

  public:
    //! Specialization of Apply() for OperatorT.
    static void Apply (const OperatorT<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& Op, 
                       const MultiVector& x, 
                       MultiVector& y, ETrans trans=NOTRANS) { 
      Op.Apply (x, y, trans); 
    }
  };

  // 
  //! @name MueLu Adapter Exceptions
  //@{

  /** \brief MueLuOpFailure is thrown when a return value from an MueLu
   * call on an Cthulhu::Operator or MueLu::Hierarchy is non-zero.
   */
  class MueLuOpFailure : public BelosError {public:
    MueLuOpFailure(const std::string& what_arg) : BelosError(what_arg)
    {}};

  // TODO: doc
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > 
  class MueLuOp : public OperatorT<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {  //TODO: should be it named CthulhuOp (if Cthulhu used by other packages) ?
    
    typedef Cthulhu::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
    typedef Cthulhu::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    MueLuOp(const Teuchos::RCP<Operator> & Op) : Op_(Op) {}
    
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
    void Apply ( const MultiVector& x, MultiVector& y, ETrans trans=NOTRANS ) const {
      TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure, 
                         "Belos::MueLuOp::Apply, transpose mode != NOTRANS not supported."); 

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Op_->apply(x,y);
    }

  private:
  
    Teuchos::RCP<Operator> Op_;
  };

  //TODO: doc
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps > 
  class MueLuPrecOp : public OperatorT<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> { 
    
    typedef MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Hierarchy;
    typedef Cthulhu::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    MueLuPrecOp(const Teuchos::RCP<Hierarchy> & H) : Hierarchy_(H) {}
    
    //! Destructor.
    virtual ~MueLuPrecOp() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Cthulu::MultiVector \c x and applies the operator
      to it resulting in the Cthulu::MultiVector \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    // Note: throw EpetraOpFailure exceptions as Belos::EpetraOp
    void Apply ( const MultiVector& x, MultiVector& y, ETrans trans=NOTRANS ) const {

      TEST_FOR_EXCEPTION(trans!=NOTRANS, MueLuOpFailure, 
                         "Belos::MueLuPrecOp::Apply, transpose mode != NOTRANS not supported by MueLu preconditionners."); 

      //FIXME InitialGuessIsZero currently does nothing in MueLu::Hierarchy.Iterate().
      y.putScalar(0.0);

      Hierarchy_->Iterate( x, 1, y , true);
      
    }
  
  private:
  
    Teuchos::RCP<Hierarchy> Hierarchy_;
  };

} // namespace Belos

#endif // BELOS_MUELU_ADAPTER_HPP
