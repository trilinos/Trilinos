// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DECL_HPP
#define MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"

#include "MueLu_ToggleCoordinatesTransferFactory_fwd.hpp"

namespace MueLu {

/*!
  @class ToggleCoordinatesTransferFactory class.
  @brief Class for transferring coordinates from a finer level to a coarser one

*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class ToggleCoordinatesTransferFactory : public TwoLevelFactoryBase {
#undef MUELU_TOGGLECOORDINATESTRANSFERFACTORY_SHORT
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
  ToggleCoordinatesTransferFactory()
    : hasDeclaredInput_(false) {}

  //! Destructor.
  virtual ~ToggleCoordinatesTransferFactory() {}

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

  //! @name Get/Set functions
  //@{

  /*! @brief Add a coordinate transfer factory in the end of list of coordinate transfer factories */
  void AddCoordTransferFactory(const RCP<const FactoryBase> &factory);

  //! Returns number of coordinate transfer factories.
  size_t NumCoordTransferFactories() const { return coordFacts_.size(); }
  //@}

 private:
  //! list of user-defined transfer coordinate factories which provide coordinates on the coarse level!
  mutable std::vector<RCP<const FactoryBase> > coordFacts_;

  mutable bool hasDeclaredInput_;
};  // class ToggleCoordinatesTransferFactory

}  // namespace MueLu

#define MUELU_TOGGLECOORDINATESTRANSFERFACTORY_SHORT
#endif  // MUELU_TOGGLECOORDINATESTRANSFER_FACTORY_DECL_HPP
