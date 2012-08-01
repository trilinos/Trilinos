/*
 * MueLu_VariableTransferFactory_decl.hpp
 *
 *  Created on: Jul 30, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_
#define MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_AggStatTransferFactory_fwd.hpp"

namespace MueLu {

  /*!
    @class AggStatTransferFactory class.
    @brief Simple class for transferring aggregation status information to next coarser level

    This is to be used in conjunction with Muelu::RAPFactory::AddTransferFactory().
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AggStatTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_AGGSTATTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Constructor.

       @param varName The name of the quantity to be restricted.
       @param genFact The generating factory
    */
    AggStatTransferFactory(std::string const & varName, RCP<const FactoryBase> const &genFact=Teuchos::null);

    //! Destructor.
    virtual ~AggStatTransferFactory();

    //@}

    //! @name Input
    //@{

    /*! @brief Specifies the data that this class needs, and the factories that generate that data.

        If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
        will fall back to the settings in FactoryManager.
    */
    void DeclareInput(Level &finelevel, Level &coarseLevel) const;

    //@}

    //! @name Build methods.
    //@{

    //! Build an object with this factory.
    void Build(Level & fineLevel, Level &coarseLevel) const;

    //@}

  private:
    //! name of variable to be transfered.
    std::string varName_;
    //! factory that generates the variable data
    RCP<const FactoryBase> genFact_;

  }; // class MultiVectorTransferFactory

} // namespace MueLu

#define MUELU_AGGSTATTRANSFERFACTORY_SHORT
#endif /* MUELU_AGGSTATTRANSFERFACTORY_DECL_HPP_ */
