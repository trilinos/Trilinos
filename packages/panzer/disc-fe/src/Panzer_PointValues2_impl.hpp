// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_POINT_VALUES2_IMPL_HPP
#define PANZER_POINT_VALUES2_IMPL_HPP

#include "Intrepid2_CellTools.hpp"

#include "Panzer_CommonArrayFactories.hpp"

// ***********************************************************
// * Evaluation and SetupArrays are NOT specialized
// ***********************************************************

namespace panzer {

  template <typename Scalar>
  void PointValues2<Scalar>::
  setupArrays(const Teuchos::RCP<const PointRule> & pr)
  {
    MDFieldArrayFactory af(prefix_, ddims_, alloc_arrays_);

    point_rule = pr;
    
    int num_nodes = point_rule->topology->getNodeCount();
    int num_cells = point_rule->workset_size;
    int num_space_dim = point_rule->spatial_dimension;

    if (point_rule->isSide()) {
       TEUCHOS_ASSERT(false); // not implemented!!!!
    }

    int num_points = point_rule->num_points;

    coords_ref = af.template buildStaticArray<Scalar,IP,Dim>("coords_ref",num_points, num_space_dim);

    node_coordinates = af.template buildStaticArray<Scalar,Cell,NODE,Dim>("node_coordinates",num_cells, num_nodes, num_space_dim);
    
    jac = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac",num_cells, num_points, num_space_dim,num_space_dim);
    jac_inv = af.template buildStaticArray<Scalar,Cell,IP,Dim,Dim>("jac_inv",num_cells, num_points, num_space_dim,num_space_dim);
    jac_det = af.template buildStaticArray<Scalar,Cell,IP>("jac_det",num_cells, num_points);
    
    point_coords = af.template buildStaticArray<Scalar,Cell,IP,Dim>("point_coords",num_cells, num_points, num_space_dim);
  }

  template <typename Scalar>
  void PointValues2<Scalar>::
  evaluateValues()
  {
    if (point_rule->isSide()) {
       TEUCHOS_ASSERT(false); // not implemented!!!!
    }
    
    Intrepid2::CellTools<PHX::exec_space> cell_tools;
    
    cell_tools.setJacobian(jac.get_view(), coords_ref.get_view(), node_coordinates.get_view(), *(point_rule->topology));
    cell_tools.setJacobianInv(jac_inv.get_view(), jac.get_view());
    cell_tools.setJacobianDet(jac_det.get_view(), jac.get_view());
    
    // IP coordinates
    cell_tools.mapToPhysicalFrame(point_coords.get_view(), coords_ref.get_view(), node_coordinates.get_view(), *(point_rule->topology));
  }

  template <typename Scalar>
  template <typename CoordinateArray>
  void PointValues2<Scalar>::
  copyNodeCoords(const CoordinateArray& in_node_coords)
  {
    // copy cell node coordinates
    {
      size_type num_cells = in_node_coords.extent(0);
      size_type num_nodes = in_node_coords.extent(1);
      size_type num_dims = in_node_coords.extent(2);
     
      for (size_type cell = 0; cell < num_cells;  ++cell)
	for (size_type node = 0; node < num_nodes; ++node)
	  for (size_type dim = 0; dim < num_dims; ++dim)
	    node_coordinates(cell,node,dim) = in_node_coords(cell,node,dim);
    }
  }

  template <typename Scalar>
  template <typename CoordinateArray>
  void PointValues2<Scalar>::
  copyPointCoords(const CoordinateArray& in_point_coords)
  {
    // copy reference point values
    {
      size_type num_points = in_point_coords.extent(0);
      size_type num_dims = in_point_coords.extent(1);
     
      for (size_type point = 0; point < num_points; ++point)
        for (size_type dim = 0; dim < num_dims; ++dim)
          coords_ref(point,dim) = in_point_coords(point,dim);
    }
  }

}

#endif
