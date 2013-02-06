/*
 * MueLu_AdaptiveSaMLParamterListInterpreter_decl.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: tobias
 */

#ifndef MUELU_ADAPTIVESAMLPARAMTERLISTINTERPRETER_DECL_HPP_
#define MUELU_ADAPTIVESAMLPARAMTERLISTINTERPRETER_DECL_HPP_

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_AdaptiveSaMLParameterListInterpreter_fwd.hpp"

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
#include "MueLu_CoupledAggregationFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_MLParameterListInterpreter_fwd.hpp"

namespace MueLu {

  /*
    Utility that from an existing Teuchos::ParameterList creates a new list, in
    which level-specific parameters are replaced with sublists.

    Currently, level-specific parameters that begin with "smoother:"
    or "aggregation:" are placed in sublists. Coarse options are also placed
    in a coarse list.

    Example:
    Input:
    smoother: type (level 0) = symmetric Gauss-Seidel
    smoother: sweeps (level 0) = 1
    Output:
    smoother: list (level 0) ->
    smoother: type = symmetric Gauss-Seidel
    smoother: sweeps = 1
  */
  // This function is a copy of ML_CreateSublists to avoid dependency on ML
  // Throw exception on error instead of exit()
  //void CreateSublists(const ParameterList &List, ParameterList &newList);


  /*!
    @class AdaptiveSAMLParameterListInterpreter class.
    @brief Class that accepts ML-style parameters and builds a MueLu preconditioner.
    This interpreter uses the same default values as ML. This allows to compare ML/MueLu results
  */

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AdaptiveSaMLParameterListInterpreter : public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    AdaptiveSaMLParameterListInterpreter() : nullspace_(NULL), blksize_(1) { }

    //! Constructor.
    //! @param paramList: parameter list with ML parameters
    //! @param nspVector: MultiVector with fine-level nullspace approximation
    //! @param factoryList: vector with RCP of FactoryBase objects
    //!
    //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
    //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
    //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    AdaptiveSaMLParameterListInterpreter(Teuchos::ParameterList & paramList,Teuchos::RCP<MultiVector> nspVector, std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));

    //! Constructor.
    //! @param xmlFileName: file name for XML file with ML parameters
    //! @param factoryList: vector with RCP of FactoryBase objects
    //!
    //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
    //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
    //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    AdaptiveSaMLParameterListInterpreter(const std::string & xmlFileName,std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));

    //! Destructor.
    virtual ~AdaptiveSaMLParameterListInterpreter() { }

    //@}

    //@{

    void SetParameterList(const Teuchos::ParameterList & paramList);

    //@}

    //@{

    //! Setup Hierarchy object
    virtual void SetupHierarchy(Hierarchy & H) const;

    //@}

    //@{

    //! @name Handling of additional user-specific transfer factories
    //@{
    /*! @brief Add transfer factory in the end of list of transfer factories for RAPFactory.

    This allows the user to add user-specific factories to the MueLu Hierarchy. The idea is to be able
    to add some factories that write out some debug information etc. which are not handled by the ML
    Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
    */
    void AddTransferFactory(const RCP<FactoryBase> & factory);

    //! Returns number of transfer factories.
    size_t NumTransferFactories() const;
    //@}

  private:

    //! build multigrid hierarchy for improving nullspace
    //! use ML settings that are also used for the final full multigrid
    //! hierarchy. In contrary to the final multigrid hierarchy use
    //! only nonsmoothed transfer operators (safe time of prolongator smoothing)
    //! and cheap level smoothers (no direct solver on coarsest level).
    void SetupInitHierarchy(Hierarchy & H) const;

    //! internal routine to add a new factory manager used for the initialization phase
    void AddInitFactoryManager(int startLevel, int numDesiredLevel, RCP<FactoryManagerBase> manager) {
      const int lastLevel = startLevel + numDesiredLevel - 1;
      if (init_levelManagers_.size() < lastLevel + 1) init_levelManagers_.resize(lastLevel + 1);

      for(int iLevel = startLevel; iLevel <= lastLevel; iLevel++) {
        init_levelManagers_[iLevel] = manager;
      }
    }

    //! Used in SetupInitHierarchy() to access levelManagers_
    //! Inputs i=-1 and i=size() are allowed to simplify calls to hierarchy->Setup()
    Teuchos::Ptr<FactoryManagerBase> InitLvlMngr(int levelID, int lastLevelID) const {

      // Please not that the order of the 'if' statements is important.

      if (levelID == -1)                    return Teuchos::null; // when this routine is called with levelID == '-1', it means that we are processing the finest Level (there is no finer level)
      if (levelID == lastLevelID+1)         return Teuchos::null; // when this routine is called with levelID == 'lastLevelID+1', it means that we are processing the last level (ie: there is no nextLevel...)

      if (0       == init_levelManagers_.size()) {                     // default factory manager.
        // the default manager is shared across levels, initialized only if needed and deleted with the HierarchyManager.
        static RCP<FactoryManagerBase> defaultMngr = rcp(new FactoryManager());
        return defaultMngr();
      }
      if (levelID >= init_levelManagers_.size()) return init_levelManagers_[init_levelManagers_.size()-1](); // last levelManager is used for all the remaining levels.

      return init_levelManagers_[levelID](); // throw exception if out of bound.
    }

    //! nullspace can be embedded in the ML parameter list
    int     nullspaceDim_;
    double* nullspace_;

    //! export aggregates
    bool    bExportAggregates_; //!< if set to true an AggregationExportFactory is used to export aggregation information (default = false)

    //! list of user-defined transfer Factories
    //! We use this vector to add some special user-given factories to the Hierarchy (RAPFactory)
    //! This way the user can extend the standard functionality of the MLParameterListInterpreter beyond the
    //! capabibilities of ML.
    std::vector<RCP<FactoryBase> > TransferFacts_;

    //! list of levelManagers for adaptive smoothed aggregation
    //! initialization phase
    Array<RCP<FactoryManagerBase> > init_levelManagers_;

    //@{ Matrix configuration

    //! Setup Matrix object
    //! overloaded from HierarchyManager to set nDofsPerNode
    virtual void SetupMatrix(Matrix & Op) const;

    //! Matrix configuration storage
    int blksize_;
    //@}

  }; // class AdaptiveSaMLParameterListInterpreter

} // namespace MueLu

#define MUELU_ADAPTIVESAMLPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_ADAPTIVESAMLPARAMTERLISTINTERPRETER_DECL_HPP_ */
