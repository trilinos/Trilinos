// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP
#define MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP

#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_Operator_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyManager.hpp"
#include "MueLu_MLParameterListInterpreter_fwd.hpp"

#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"

#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_PgPFactory_fwd.hpp"
#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_GenericRFactory_fwd.hpp"
#include "MueLu_SmootherPrototype_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_IfpackSmoother_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

#if defined(HAVE_MUELU_ISORROPIA) && defined(HAVE_MPI)
#include "MueLu_RepartitionHeuristicFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_RebalanceTransferFactory_fwd.hpp"
#include "MueLu_IsorropiaInterface_fwd.hpp"
#include "MueLu_RebalanceAcFactory_fwd.hpp"
#include "MueLu_RebalanceMapFactory_fwd.hpp"
#endif

#ifdef HAVE_MUELU_DEPRECATED_CODE
#ifdef MueLu_SHOW_DEPRECATED_WARNINGS
#warning "The header file MueLu_MLParameterListInterpreter.hpp is deprecated"
#endif
#else
#error "The header file MueLu_MLParameterListInterpreter.hpp is deprecated"
#endif

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
void CreateSublists(const ParameterList& List, ParameterList& newList);

/*!
  @class MLParameterListInterpreter class.
  @brief Class that accepts ML-style parameters and builds a MueLu preconditioner.
  This interpreter uses the same default values as ML. This allows to compare ML/MueLu results

  The parameter list is validated only if the package ML is available and parameter "ML validate parameter list" is true.
  TODO: A warning is issued if ML is not available
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class MLParameterListInterpreter : public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
#undef MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  MLParameterListInterpreter()
    : nullspace_(NULL)
    , blksize_(1) {}

  //! Constructor.
  //! @param paramList: parameter list with ML parameters
  //! @param[in] comm  (RCP<Teuchos::Comm<int> >): Optional RCP of a Teuchos communicator  (default: Teuchos::null)
  //! @param factoryList: vector with RCP of FactoryBase objects
  //!
  //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
  //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
  //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
  MLParameterListInterpreter(Teuchos::ParameterList& paramList, Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::null, std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));

  //! Constructor.
  //! @param xmlFileName: file name for XML file with ML parameters
  //! @param factoryList: vector with RCP of FactoryBase objects
  //!
  //! The factories in factoryList allow the user to add user-specific factories to the MueLu Hierarchy.
  //! The idea is to be able to add some factories that write out some debug information etc. which are not handled by the ML
  //! Parameter List itself. See information about the RAPFactory::AddTransferFactory method, too!
  MLParameterListInterpreter(const std::string& xmlFileName, std::vector<RCP<FactoryBase> > factoryList = std::vector<RCP<FactoryBase> >(0));

  //! Destructor.
  virtual ~MLParameterListInterpreter() = default;

  //@}

  //@{

  void SetParameterList(const Teuchos::ParameterList& paramList);

  //@}

  //@{

  //! Setup Hierarchy object
  virtual void SetupHierarchy(Hierarchy& H) const;

  //@}

  //@{

  //! @name static helper functions translating parameter list to factories
  //! @brief static helper functions that also can be used from outside for translating ML parameters into MueLu objects
  //@{

  //! Read smoother options and build the corresponding smoother factory
  // @param AFact: Factory used by smoother to find 'A'
  static RCP<SmootherFactory> GetSmootherFactory(const Teuchos::ParameterList& paramList, const RCP<FactoryBase>& AFact = Teuchos::null);

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
  int nullspaceDim_;
  double* nullspace_;  // TODO: replace by Teuchos::ArrayRCP<>

  //! coordinates can be embedded in the ML parameter list
  double* xcoord_;
  double* ycoord_;
  double* zcoord_;

  //! list of user-defined transfer Factories
  //! We use this vector to add some special user-given factories to the Hierarchy (RAPFactory)
  //! This way the user can extend the standard functionality of the MLParameterListInterpreter beyond the
  //! capabibilities of ML.
  std::vector<RCP<FactoryBase> > TransferFacts_;

  //@{ Matrix configuration

  //! Setup Operator object
  virtual void SetupOperator(Operator& Op) const;

  //! Matrix configuration storage
  int blksize_;

  //@}

};  // class MLParameterListInterpreter

}  // namespace MueLu

#define MUELU_MLPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_MLPARAMETERLISTINTERPRETER_DECL_HPP */
