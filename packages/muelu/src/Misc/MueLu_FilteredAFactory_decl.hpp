// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_FILTEREDAFACTORY_DECL_HPP
#define MUELU_FILTEREDAFACTORY_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FilteredAFactory_fwd.hpp"

#include "MueLu_LWGraph.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
namespace MueLu {

/*!
  @class FilteredAFactory class.
  @brief Factory for building filtered matrices using filtered graphs.
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class FilteredAFactory : public SingleLevelFactoryBase {
#undef MUELU_FILTEREDAFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  FilteredAFactory() {}

  //! Destructor.
  virtual ~FilteredAFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds filtered matrix and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}
 private:
  void BuildReuse(const Matrix& A, const LWGraph& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const;
  void BuildNew(const Matrix& A, const LWGraph& G, const bool lumping, double dirichletThresh, Matrix& filteredA) const;
  void BuildNewUsingRootStencil(const Matrix& A, const LWGraph& G, double dirichletThresh, Level& currentLevel, Matrix& filteredA, bool use_spread_lumping, double DdomAllowGrowthRate, double DdomCap) const;
  void ExperimentalLumping(const Matrix& A, Matrix& filteredA, double rho, double rho2) const;

};  // class FilteredAFactory

}  // namespace MueLu

#define MUELU_FILTEREDAFACTORY_SHORT
#endif  // MUELU_FILTEREDAFACTORY_DECL_HPP
