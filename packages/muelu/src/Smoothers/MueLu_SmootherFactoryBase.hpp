#ifndef MUELU_SMOOTHERFACTORYBASE_HPP
#define MUELU_SMOOTHERFACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"

#include <iostream>

//#include "MueLu_Smoother.hpp"
#include "MueLu_SmootherPrototype.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

/*!
  @class Smoother factory base class.
  @brief Base class for smoother factories.

  Has only one real capability, which is to record (for the Hierarchy class)
  whether a smoother should be built on the coarsest level.
*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class SmootherFactoryBase {

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
    //virtual bool Build(Level &level) = 0;
    virtual void Build(RCP<Level> level, RCP<SmootherPrototype> &preSmoo, RCP<SmootherPrototype> &postSmoo) const = 0;

    //@}

    //! @name Set/Get methods.
    //@{

    //! Set whether a smoother should be built on the coarsest level.
    void SetCoarsest(bool ToF) {
      DoCoarsest_ = ToF;
    }

    //! If true, a smoother should be built on the coarsest level.
    bool DoCoarsest() {
      return DoCoarsest_;
    }
    //@}

}; //class SmootherFactoryBase

} //namespace MueLu

#define MUELU_SMOOTHERFACTORYBASE_SHORT

#endif //ifndef MUELU_SMOOTHERFACTORYBASE_HPP
