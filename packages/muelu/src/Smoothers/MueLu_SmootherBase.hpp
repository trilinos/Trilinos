#ifndef MUELU_SMOOTHERBASE_HPP
#define MUELU_SMOOTHERBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

  /*!
    @class SmootherBase
    @brief Base class for smoothers

    This has the signature for the required Apply function and contains data that is generic across all
    smoothers.
  */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherBase : public Teuchos::VerboseObject<SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > {

#include "MueLu_UseShortNames.hpp"

  private:
    std::string type_;

  public:
    //@{ Constructors/Destructors.
    SmootherBase() {}

    virtual ~SmootherBase() {}
    //@}

    //! @name Apply methods.
    //@{

    //! Apply smoother.
    //virtual void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero) = 0;
    virtual void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) = 0;

    //@}

    //! @name Set/Get methods.
    //@{

    //! Prototype of the method to setup the number of iteration of the smoother.
    virtual void SetNIts(LO const &Nits) { //TODO: should be removed from the interface
       TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "SetNIts() not implemented for this smoother");
    };

    //! Get the smoother type.
    std::string GetType() const {
      return type_;
    }

    /*! @brief Set the smoother type.

    This method must be called by constructors of derived classes.
    */
    void SetType(std::string type) {
      type_ = type;
    }

    //! @name Utilities.
    //@{

    //! @brief Print information about the smoother
    virtual void Print(std::string prefix) const = 0;

    //@}

  }; //class SmootherBase

} //namespace MueLu

#define MUELU_SMOOTHERBASE_SHORT

#endif //ifndef MUELU_SMOOTHERBASE_HPP
