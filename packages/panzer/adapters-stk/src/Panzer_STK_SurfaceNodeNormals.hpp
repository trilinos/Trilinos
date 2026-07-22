// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_SurfaceNormals_hpp__
#define __Panzer_STK_SurfaceNormals_hpp__

#include "PanzerAdaptersSTK_config.hpp"
#include "Teuchos_RCP.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include <unordered_map>
#include <string>
#include <vector>

namespace panzer_stk {

  class STK_Interface;
  

  /** \brief Computes the normals for all nodes associated with a sideset surface

      Computes the node normals for a given side set.  This
      computation will use ghosted element contributions for all
      surfaces associated with the node, but will only use elements
      associated with the sideset.  So if two sidesets adjoin at the
      seam, the nodes on the seam will NOT have the same node normal.
      A simple future addition would be to allow for a list of
      sidesets to be provided to allow for the union computation.

      \param[out] normals Map of the node normals (including ghosted nodes).  Key is the global id from the stk mesh node entity id and value is vector of the size of the parent element dimension containing the normal vector components.  Vector components will be resized on calling this method.
      \param[in] mesh (Required) Panzer stk mesh 
      \param[in] sidesetName (Required) Name of the sideset that the normals will be computed on
      \param[in] elementBlockName (Required) Name of the element block that the outward facing normals will be computed on
      \param[in] out (Optional) The ostream used for serial debug output on print process only.  If non-null this will print debug info.
      \param[in] pout (Optional) The ostream used for parallel debug output by all processes.  If non-null this will print debug info.
  */
  void computeSidesetNodeNormals(std::unordered_map<unsigned,std::vector<double> >& normals,
				 const Teuchos::RCP<const panzer_stk::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* out = NULL,
				 std::ostream* pout = NULL);

  /** \brief Computes the normals for all nodes associated with a sideset surface

      Computes the node normals for a given side set.  This
      computation will use ghosted element contributions for all
      surfaces associated with the node, but will only use elements
      associated with the sideset.  So if two sidesets adjoin at the
      seam, the nodes on the seam will NOT have the same node normal.
      A simple future addition would be to allow for a list of
      sidesets to be provided to allow for the union computation.

      \param[out] normals Map of the node normals (including ghosted nodes).  Key is the panzer_stk::STK_Interface local element id and value is a multidimensional array of the size of the number of nodes times the parent element dimension containing the normal vector components.  Vector components will be resized on calling this method.
      \param[in] mesh (Required) Panzer stk mesh 
      \param[in] sidesetName (Required) Name of the sideset that the normals will be computed on
      \param[in] elementBlockName (Required) Name of the element block that the outward facing normals will be computed on
      \param[in] out (Optional) The ostream used for serial debug output on print process only.  If non-null this will print debug info.
      \param[in] pout (Optional) The ostream used for parallel debug output by all processes.  If non-null this will print debug info.
  */
  void computeSidesetNodeNormals(std::unordered_map<std::size_t,Kokkos::DynRankView<double,PHX::Device> >& elementToNormalMap,
				 const Teuchos::RCP<const panzer_stk::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* out = NULL,
				 std::ostream* pout = NULL);


}

#endif
