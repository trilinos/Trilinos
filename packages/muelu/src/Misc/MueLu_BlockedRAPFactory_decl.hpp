// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BLOCKEDRAPFACTORY_DECL_HPP
#define MUELU_BLOCKEDRAPFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_RAPFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"
#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

namespace MueLu {
/*!
  @class BlockedRAPFactory
  @brief Factory for building coarse matrices.
*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class BlockedRAPFactory : public TwoLevelFactoryBase {
#undef MUELU_BLOCKEDRAPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  BlockedRAPFactory();

  virtual ~BlockedRAPFactory() = default;
  //@}

  //! @name Input
  //@{

  RCP<const ParameterList> GetValidParameterList() const override;

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const override;

  //@}

  //! @name Build methods.
  //@{
  void Build(Level &fineLevel, Level &coarseLevel) const override;
  //@}

  //! @name Handling of user-defined transfer factories
  //@{

  //! Indicate that zero entries on the diagonal of Ac shall be repaired (i.e. if A(i,i) == 0.0 set A(i,i) = 1.0)
  void SetRepairZeroDiagonal(bool const &repair) {
    repairZeroDiagonals_ = repair;
    if (repair) checkAc_ = true;  // make sure that plausibility check is performed. Otherwise SetRepairZeroDiagonal(true) has no effect.
  }

  //! Indicate that a simple plausibility check shall be done for Ac after building RAP
  void SetPlausibilityCheck(bool const &check) {
    checkAc_ = check;
  }

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
  //! @name internal plausibility check methods
  //! checks main diagonal entries of (0,0) block. Does not affect entries in (1,1) block!
  static void CheckMainDiagonal(RCP<BlockedCrsMatrix> &bAc, bool repairZeroDiagonals = false);

  //! If true, perform a basic plausibility check on Ac (default = false)
  //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
  bool checkAc_;

  //! If true, the CheckMainDiagonal routine automatically repairs zero entries on main diagonal (default = false)
  //! i.e. if A(i,i) == 0.0 set A(i,i) = 1.0
  //! note, that the repairZeroDiagonals_ flag only is valid for checkAc_ == true
  bool repairZeroDiagonals_;

  //@{

  //! list of user-defined transfer Factories
  std::vector<RCP<const FactoryBase> > transferFacts_;

  //@}

};  // class BlockedRAPFactory

}  // namespace MueLu

#define MUELU_BLOCKEDRAPFACTORY_SHORT
#endif  // MUELU_BLOCKEDRAPFACTORY_DECL_HPP
