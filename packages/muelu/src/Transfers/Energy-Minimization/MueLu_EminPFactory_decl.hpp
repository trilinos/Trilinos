// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_EMINPFACTORY_DECL_HPP
#define MUELU_EMINPFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_StridedMapFactory_fwd.hpp>

#include "MueLu_EminPFactory_fwd.hpp"

#include "MueLu_CGSolver_fwd.hpp"
#include "MueLu_Constraint_fwd.hpp"
#include "MueLu_GMRESSolver_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_SolverBase_fwd.hpp"
#include "MueLu_SteepestDescentSolver_fwd.hpp"

namespace MueLu {

/*!
  @class EminPFactory class.
  @brief Factory for building Energy Minimization prolongators.
  @ingroup MueLuTransferClasses
  */

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class EminPFactory : public PFactory {
#undef MUELU_EMINPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! @brief Constructor.
  EminPFactory() {}

  //! Destructor.
  virtual ~EminPFactory() {}

  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! @name Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds energy minimization prolongator and returns it in <tt>coarseLevel</tt>.
    */
  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

};  // class EminPFactory

}  // namespace MueLu

#define MUELU_EMINPFACTORY_SHORT
#endif  // MUELU_EMINPFACTORY_DECL_HPP
