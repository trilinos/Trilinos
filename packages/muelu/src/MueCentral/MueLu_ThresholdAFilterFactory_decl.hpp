#ifndef MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
#define MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP

/*
 * MueLu_ThresholdAFilterFactory.hpp
 *
 *  Created on: 14.10.2011
 *      Author: tobias
 */

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_ThresholdAFilterFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class ThresholdAFilterFactory class.
    @brief Factory for building a thresholded operator.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ThresholdAFilterFactory : public SingleLevelFactoryBase {
#undef MUELU_THRESHOLDAFILTERFACTORY_SHORT
    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ThresholdAFilterFactory(const std::string& ename, const FactoryBase* fac, const Scalar threshold);

    //! Destructor.
    virtual ~ThresholdAFilterFactory();
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
    std::string        varName_;   ///< name of input and output variable
    const FactoryBase* factory_;   ///< generating factory of input variable
    const Scalar       threshold_; ///< threshold parameter


  }; // class ThresholdAFilterFactory

} // namespace MueLu

#define MUELU_THRESHOLDAFILTERFACTORY_SHORT
#endif // MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
