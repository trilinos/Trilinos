// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_RAPSHIFTFACTORY_DECL_HPP
#define MUELU_RAPSHIFTFACTORY_DECL_HPP

#include <string>

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory2_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_RAPShiftFactory_fwd.hpp"

#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {
/*!
  @class RAPShiftFactory
  @brief Factory for building coarse grid matrices, when the matrix
         is of the form K+a*M. Useful when you want to change the shift
         variable ("a") at every level. Each level must store the stiffness
         matrix K and mass matrix M separately.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RAPShiftFactory : public TwoLevelFactoryBase {
#undef MUELU_RAPSHIFTFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  RAPShiftFactory();

  virtual ~RAPShiftFactory() {}

  //@}

  //! @name Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &fineLevel, Level &coarseLevel) const;
  //@}

  //! @name Handling of user-defined transfer factories
  //@{

  //! Indicate that the restriction operator action should be implicitly defined by the transpose of the prolongator.
  void SetImplicitTranspose(bool const &implicit) {
    implicitTranspose_ = implicit;
  }

  void SetShifts(std::vector<Scalar> &shifts) {
    shifts_.clear();
    shifts_ = shifts;
  }

  //@}

  //@{
  /*! @brief Add transfer factory in the end of list of transfer factories in RepartitionAcFactory.

  Transfer factories are derived from TwoLevelFactoryBase and project some data from the fine level to
  the next coarser level.
  */
  void AddTransferFactory(const RCP<const FactoryBase> &factory);

  // TODO add a function to remove a specific transfer factory?

  //! Returns number of transfer factories.
  size_t NumTransferFactories() const { return transferFacts_.size(); }

  //@}

 private:
  //! If true, the action of the restriction operator action is implicitly defined by the transpose of the prolongator.
  bool implicitTranspose_;

  //! list of user-defined transfer Factories
  std::vector<RCP<const FactoryBase> > transferFacts_;

  // vector of shifting terms
  std::vector<Scalar> shifts_;

};  // class RAPShiftFactory

}  // namespace MueLu

#define MUELU_RAPSHIFTFACTORY_SHORT
#endif  // MUELU_RAPSHIFTFACTORY_DECL_HPP
