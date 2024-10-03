// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP
#define MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Types.hpp"

#include <Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_Core.hpp>

#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_RegionRFactory_kokkos_fwd.hpp"

namespace MueLu {

/*!
  @class RegionRFactory_kokkos class
  @brief Factory that builds a restriction operator for region multigrid
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class RegionRFactory_kokkos : public TwoLevelFactoryBase {
#undef MUELU_REGIONRFACTORY_KOKKOS_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  using real_type                  = typename Teuchos::ScalarTraits<SC>::coordinateType;
  using realvaluedmultivector_type = typename Xpetra::MultiVector<real_type, LO, GO, NO>;
  using execution_space            = typename Node::execution_space;
  using memory_space               = typename Node::memory_space;
  using device_type                = Kokkos::Device<execution_space, memory_space>;
  using intTupleView               = typename Kokkos::View<int[3], device_type>;
  using LOTupleView                = typename Kokkos::View<LO[3], device_type>;

  //! @name Constructors/Destructors.
  //@{

  //! Default Constructor
  RegionRFactory_kokkos() = default;

  //! Destructor
  virtual ~RegionRFactory_kokkos() = default;
  //@}

  //! Input
  //@{
  RCP<const ParameterList> GetValidParameterList() const;

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;

  void Build3D(const int numDimensions,
               Array<LO>& lFineNodesPerDim,
               const RCP<Matrix>& A,
               const RCP<realvaluedmultivector_type>& fineCoordinates,
               RCP<Matrix>& R,
               RCP<realvaluedmultivector_type>& coarseCoordinates,
               Array<LO>& lCoarseNodesPerDim) const;

  //@}

};  // class RegionRFactory_kokkos

}  // namespace MueLu

#define MUELU_REGIONRFACTORY_KOKKOS_SHORT
#endif  // MUELU_REGIONRFACTORY_KOKKOS_DECL_HPP
