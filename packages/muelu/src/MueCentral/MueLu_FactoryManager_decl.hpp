// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FACTORYMANAGER_DECL_HPP
#define MUELU_FACTORYMANAGER_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryManagerBase.hpp"

#include "MueLu_AmalgamationFactory_fwd.hpp"
#include "MueLu_AggregateQualityEstimateFactory_fwd.hpp"
#include "MueLu_CoalesceDropFactory_fwd.hpp"
#include "MueLu_CoarseMapFactory_fwd.hpp"
#include "MueLu_ConstraintFactory_fwd.hpp"
#include "MueLu_DirectSolver_fwd.hpp"
#include "MueLu_InitialBlockNumberFactory_fwd.hpp"
#include "MueLu_LineDetectionFactory_fwd.hpp"
#include "MueLu_NullspaceFactory_fwd.hpp"
#include "MueLu_PatternFactory_fwd.hpp"
#include "MueLu_RAPFactory_fwd.hpp"
#include "MueLu_RepartitionHeuristicFactory_fwd.hpp"
#include "MueLu_RepartitionFactory_fwd.hpp"
#include "MueLu_SaPFactory_fwd.hpp"
#include "MueLu_ScaledNullspaceFactory_fwd.hpp"
#include "MueLu_SmootherFactory_fwd.hpp"
#include "MueLu_StructuredAggregationFactory_fwd.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_TransPFactory_fwd.hpp"
#include "MueLu_TrilinosSmoother_fwd.hpp"
#include "MueLu_UncoupledAggregationFactory_fwd.hpp"
#include "MueLu_ZoltanInterface_fwd.hpp"
#include "MueLu_InterfaceMappingTransferFactory_fwd.hpp"
#include "MueLu_InterfaceAggregationFactory_fwd.hpp"

#include "MueLu_CoalesceDropFactory_kokkos_fwd.hpp"
#include "MueLu_TentativePFactory_kokkos_fwd.hpp"

namespace MueLu {

/*!
  @class FactoryManager class.
  @brief This class specifies the default factory that should generate some data on a Level if the data does not exist and
  the generating factory has not been specified.

  Consider the following example.

  @code
    RCP<SingleLevelFactory> Afact;
    Level currentLevel;
    RCP<Matrix> thisLevelA;
    thisLevelA = currentLevel.Get<Matrix>("A", Afact.get());
  @endcode

  @todo If Afact is null (actually, Teuchos::null), then the FactoryManager associated with currentLevel will determine whether a default factory has
  been specified for creating A.  If "yes", then that factory will be called, A will be stored in currentLevel, and an RCP will be returned by
  the Get call. If "no", then the FactoryManager will <b>throw an exception indicating that it does not know how to generate A</b>.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FactoryManager : public FactoryManagerBase {
#undef MUELU_FACTORYMANAGER_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructor/Destructors
  //@{

  //! @brief Constructor.
  FactoryManager();

  //! Constructor used by HierarchyFactory (temporary, will be removed)
  FactoryManager(const std::map<std::string, RCP<const FactoryBase> >& factoryTable);

  //! Destructor.
  virtual ~FactoryManager();

  //@}

  //! @name Get/Set functions.
  //@{

  /*! @brief Set Factory

      Register the factory that should generate data if said factory is not specified in the request.

      @param[in] name of variable
      @param[in] factory that generates the data
  */
  void SetFactory(const std::string& varName, const RCP<const FactoryBase>& factory);

  /*! @brief Get factory associated with a particular data name.

     @param[in] varName name of variable.

  */
  const RCP<const FactoryBase> GetFactory(const std::string& varName) const;

  /*! @brief Get factory associated with a particular data name (NONCONST version)

     @param[in] varName name of variable.

  */
  const RCP<FactoryBase> GetFactoryNonConst(const std::string& varName);

  //! Check
  // Return true if Factory associated with varName is registered
  bool hasFactory(const std::string& varName) const;

  //!
  const RCP<const FactoryBase> GetDefaultFactory(const std::string& varName) const;

  void SetKokkosRefactor(const bool useKokkos) {
    useKokkos_ = useKokkos;
  }

  bool GetKokkosRefactor() const { return useKokkos_; }

  //@}

  void Clean() const { defaultFactoryTable_.clear(); }

#ifdef HAVE_MUELU_DEBUG
  void ResetDebugData() const;
#endif

  void Print() const;

 private:
  //! @name Helper functions
  //@{

  /*! Add a factory to the default factory list and return it. This helper function is used by GetDefaultFactory()

   @todo TODO factory->setObjectLabel("Default " + varName + "Factory");
  */

  const RCP<const FactoryBase> SetAndReturnDefaultFactory(const std::string& varName, const RCP<const FactoryBase>& factory) const;
  //@}

  /*! @brief User-defined factories.
   *
   * User may overwrite default behaviour. The user provided factories are stored in a separate table. When we try to determine
   * which factory generates the data, this table is searched first.

    Note: we distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
  */
  std::map<std::string, RCP<const FactoryBase> > factoryTable_;

  /*! @brief Table that holds default factories.

    -# We distinguish 'user defined factory' and 'default factory' to allow the deallocation of default factories separately.
    -# <tt>defaultFactoryTable_</tt> is mutable because default factories are only added to the list when they are requested
    to avoid allocation of unused factories.
  */
  mutable std::map<std::string, RCP<const FactoryBase> > defaultFactoryTable_;

  //! Whether or not to use kokkos factories.
  bool useKokkos_;

};  // class

}  // namespace MueLu

#define MUELU_FACTORYMANAGER_SHORT
#endif  // MUELU_FACTORYMANAGER_DECL_HPP
