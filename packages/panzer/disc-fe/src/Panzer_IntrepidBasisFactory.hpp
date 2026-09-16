// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_INTREPID_BASIS_FACTORY_H
#define PANZER_INTREPID_BASIS_FACTORY_H

#include <sstream>
#include <string>

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Teuchos_RCP.hpp"
#include "Intrepid2_Basis.hpp"

#include "Shards_CellTopology.hpp"


namespace panzer {


  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName() 
                               after trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename DeviceType, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<DeviceType,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                       const shards::CellTopology & cell_topology);

  /** \brief Creates an Intrepid2::Basis object given the basis, order and cell topology.

      \param[in] basis_type The name of the basis.
      \param[in] basis_order The order of the polynomial used to construct the basis.
      // TODO BWR Is the cell_topology documentation below correct?
      \param[in] cell_topology Cell topology for the basis.  Taken from shards::CellTopology::getName() after
                               trimming the extended basis suffix.

      To be backwards compatible, this method takes deprecated
      descriptions and transform it into a valid type and order.  For
      example "Q1" is transformed to basis_type="HGrad",basis_order=1.

      \returns A newly allocated panzer::Basis object.
  */
  template <typename DeviceType, typename OutputValueType, typename PointValueType>
  Teuchos::RCP<Intrepid2::Basis<DeviceType,OutputValueType,PointValueType> >
  createIntrepid2Basis(const std::string basis_type, int basis_order,
                       const Teuchos::RCP<const shards::CellTopology> & cell_topology);
}

#endif
