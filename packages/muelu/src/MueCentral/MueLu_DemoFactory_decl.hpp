/*
 * MueLu_DemoFactory_decl.hpp
 *
 *  Created on: 06.12.2011
 *      Author: tobias
 */

#ifndef MUELU_DEMOFACTORY_DECL_HPP_
#define MUELU_DEMOFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"  // TODO: include base factory (in this case DemoFactory derives from SingleLevelFactoryBase)
#include "MueLu_DemoFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class DemoFactory class.
    @brief empty factory for demonstration

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ThresholdAFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_DEMOFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    DemoFactory();

    //! Destructor.
    virtual ~DemoFactory();
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const;

    //@}

  private:
    // TODO add member variables

  }; // class DemoFactory

} // namespace MueLu

#define MUELU_DEMOFACTORY_SHORT


#endif /* MUELU_DEMOFACTORY_DECL_HPP_ */
