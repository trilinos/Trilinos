// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REBALANCETRANSFERFACTORY_DECL_HPP
#define MUELU_REBALANCETRANSFERFACTORY_DECL_HPP

#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include "Xpetra_Vector_fwd.hpp"
#include "Xpetra_VectorFactory_fwd.hpp"
#include "Xpetra_MultiVector_fwd.hpp"
#include "Xpetra_MultiVectorFactory_fwd.hpp"
#include "Xpetra_Import_fwd.hpp"
#include "Xpetra_ImportFactory_fwd.hpp"

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_RebalanceTransferFactory_fwd.hpp"

#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Types.hpp"

namespace MueLu {

/*!
  @class RebalanceTransferFactory class.
  @brief Applies permutation to grid transfer operators.
  @ingroup MueLuTransferClasses
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RebalanceTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_REBALANCETRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  RebalanceTransferFactory();

  //! Destructor.
  virtual ~RebalanceTransferFactory();

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Input
  //@{

  /*! @brief Specifies the data that this class needs, and the factories that generate that data.

      If the Build method of this class requires some data, but the generating factory is not specified in DeclareInput, then this class
      will fall back to the settings in FactoryManager.
  */
  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  //! Build an object with this factory.
  void Build(Level& fineLevel, Level& coarseLevel) const;

  //@}

};  // class RebalanceTransferFactory

}  // namespace MueLu

#define MUELU_REBALANCETRANSFERFACTORY_SHORT
#endif  // MUELU_REBALANCETRANSFERFACTORY_DECL_HPP
