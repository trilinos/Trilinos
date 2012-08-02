/*
 * MueLu_Interpreter_decl.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: wiesner
 */

#ifndef MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP
#define MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Operator_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_MLParameterListInterpreter_fwd.hpp"

#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"

#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_PgPFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_GenericRFactory_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_IfpackSmoother_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_HierarchyHelpers_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_UCAggregationFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

namespace MueLu {

  /*!
    @class MLParameterListInterpreter class.
    @brief Class that accepts ML-style parameters and builds a MueLu preconditioner.
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class MLParameterListInterpreter : public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> { 
#undef MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    MLParameterListInterpreter() : nullspace_(NULL), blksize_(1) { }

    //! Constructor.
    //! @param paramList: parameter list with ML parameters
    //! @param factoryList: vector with RCP of FactoryBase objects
    //!
    //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
    //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
    //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    MLParameterListInterpreter(Teuchos::ParameterList & paramList,std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));

    //! Constructor.
    //! @param xmlFileName: file name for XML file with ML parameters
    //! @param factoryList: vector with RCP of FactoryBase objects
    //!
    //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
    //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
    //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    MLParameterListInterpreter(const std::string & xmlFileName,std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));
    
    //! Destructor.
    virtual ~MLParameterListInterpreter() { }

    //@}

    //@{

    void SetParameterList(const Teuchos::ParameterList & paramList);

    //@}

    //@{

    //! Setup Hierarchy object
    virtual void SetupHierarchy(Hierarchy & H) const;

    //@}

    //@{

    //! @name static helper functions translating parameter list to factories
    //! @brief static helper functions that also can be used from outside for translating ML parameters into MueLu objects
    //@{

    //! Read coarse solver options and build the corresponding smoother factory
    static RCP<SmootherFactory> GetCoarsestSolverFactory(const Teuchos::ParameterList & params, const RCP<FactoryBase> & AFact = Teuchos::null);

    //! Read smoother options and build the corresponding smoother factory
    static RCP<SmootherFactory> GetSmootherFactory(const Teuchos::ParameterList & params, int level, const RCP<FactoryBase> & AFact = Teuchos::null);

    //@}


    //! @name Handling of additional user-specific transfer factories
    //@{
    /*! @brief Add transfer factory in the end of list of transfer factories for RAPFactory.

    This allows the user to add user-specific factories to the MueLu Hierarchy. The idea is to be able
    to add some factories that write out some debug information etc. which are not handled by the ML
    Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    */
    void AddTransferFactory(const RCP<FactoryBase>& factory);

    //! Returns number of transfer factories.
    size_t NumTransferFactories() const;
    //@}

  private:

    //! nullspace can be embedded in the ML parameter list
    int     nullspaceDim_;
    double* nullspace_; //TODO: replace by Teuchos::ArrayRCP<>

    //! export aggregates
    bool    bExportAggregates_; //!< if set to true an AggregationExportFactory is used to export aggregation information (default = false)

    //! list of user-defined transfer Factories
    //! We use this vector to add some special user-given factories to the Hierarchy (RAPFactory)
    //! This way the user can extend the standard functionality of the MLParameterListInterpreter beyond the
    //! capabibilities of ML.
    std::vector<RCP<FactoryBase> > TransferFacts_;

    //@{ Operator configuration

    //! Setup Operator object
    virtual void SetupOperator(Operator & Op) const;

    //! Operator configuration storage
    int blksize_;

    //@}

  }; // class MLParameterListInterpreter

} // namespace MueLu

#define MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP */
