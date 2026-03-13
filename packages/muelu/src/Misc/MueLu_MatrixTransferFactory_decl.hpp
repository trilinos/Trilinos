// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_MATRIXTRANSFERFACTORY_DECL_HPP
#define MUELU_MATRIXTRANSFERFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "Xpetra_Matrix_fwd.hpp"
#include "MueLu_MatrixTransferFactory_fwd.hpp"

namespace MueLu {

/*!
  @class MatrixTransferFactory class.
  @brief Class for restricting a Matrix from a finer to a coarser level.

  This is to be used in conjunction with Muelu::RAPFactory::AddTransferFactory().
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class MatrixTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_MATRIXTRANSFERFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.

     @param vectorName The name of the quantity to be restricted.
     @param restrictionName The name of the restriction Matrix.

     The operator associated with <tt>projectionName</tt> will be applied to the MultiVector associated with
     <tt>vectorName</tt>.
  */
  MatrixTransferFactory() = default;

  //! Destructor.
  virtual ~MatrixTransferFactory() = default;

  RCP<const ParameterList> GetValidParameterList() const;

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
  void Build(Level &fineLevel, Level &coarseLevel) const;

  //@}

};  // class MatrixTransferFactory

}  // namespace MueLu

#define MUELU_MATRIXTRANSFERFACTORY_SHORT
#endif  // MUELU_MATRIXTRANSFER_FACTORY_DECL_HPP
