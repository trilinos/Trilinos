#ifndef MUELU_SMOOTHERFACTORYBASE_HPP
#define MUELU_SMOOTHERFACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {
  class Level;

  /*!
    @class Smoother factory base class.
    @brief Base class for smoother factories.

    Has only one real capability, which is to record (for the Hierarchy class)
    whether a smoother should be built on the coarsest level.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class SmootherFactoryBase : public SingleLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"

  private:
    bool DoCoarsest_;

  public:
    //@{ Constructors/Destructors.
    SmootherFactoryBase() {}

    virtual ~SmootherFactoryBase() {}
    //@}

    //! @name Build methods.
    //@{

    //! Build pre-smoother and/or post-smoother
    virtual bool Build(Level & currentLevel) const = 0;
  
    virtual bool BuildSmoother(Level & currentLevel, PreOrPost const &pop = BOTH) const = 0;
    //@}

  }; //class SmootherFactoryBase

} //namespace MueLu

#define MUELU_SMOOTHERFACTORYBASE_SHORT

#endif //ifndef MUELU_SMOOTHERFACTORYBASE_HPP

//TODO: remove this interface?
