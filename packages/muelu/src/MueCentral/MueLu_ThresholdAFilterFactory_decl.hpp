#ifndef MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
#define MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

/*
 * MueLu_ThresholdAFilterFactory.hpp
 *
 *  Created on: 14.10.2011
 *      Author: tobias
 */

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"



namespace MueLu {

  /*!
    @class ThresholdAFilterFactory class.
    @brief Factory for building a thresholded operator.

  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ThresholdAFilterFactory : public SingleLevelFactoryBase {

    #include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    ThresholdAFilterFactory(const std::string& ename, const FactoryBase* fac, const Scalar threshold)
      : varName_(ename), factory_(fac), threshold_(threshold)
    ;

    //! Destructor.
    virtual ~ThresholdAFilterFactory() ;
    //@}

    //! Input
    //@{

    void DeclareInput(Level &currentLevel) const ;

    //@}

    //@{
    //! @name Build methods.

    //! Build an object with this factory.
    void Build(Level & currentLevel) const ;

    //@}

  private:
    std::string         varName_;   ///< name of input and output variable
    const FactoryBase*  factory_;   ///< generating factory of input variable
    const Scalar        threshold_; ///< threshold parameter


  }; // class ThresholdAFilterFactory

} //#ifndef MUELU_THRESHOLDAFILTERFACTORY_HPP_namespace MueLu

#define MUELU_THRESHOLDAFILTERFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_THRESHOLDAFILTERFACTORY_DECL_HPP
