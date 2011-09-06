#ifndef MUELU_SMOOTHERBASE_HPP
#define MUELU_SMOOTHERBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

  /*!
    @class SmootherBase
    @brief Base class for smoothers

    This has the signature for the required Apply function and contains data that is generic across all
    smoothers.
  */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherBase : public BaseClass {

#include "MueLu_UseShortNames.hpp"

  public:
    //@{ Constructors/Destructors.
    SmootherBase() {}

    virtual ~SmootherBase() {}
    //@}

    //! @name Apply methods.
    //@{

    //! Apply smoother.
    virtual void Apply(MultiVector &x, MultiVector const &rhs, bool const &InitialGuessIsZero) = 0;

    //@}

  }; //class SmootherBase

} //namespace MueLu

#define MUELU_SMOOTHERBASE_SHORT

#endif //ifndef MUELU_SMOOTHERBASE_HPP


// SmootherBase = Interface used by Hierarchy.Iterate(). Minimal condition to be used by a smoother
// SmootherPrototype = minimal interface used by the generic SmootherFactory.
// Note that one can implements and use his own SmootherFactory. In this case, SmootherBase is enough.
// AdvSmootherPrototype = for more complex case of reusing setup between presmoother and postsmoother
