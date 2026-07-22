// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_POINT_VALUES2_IMPL_HPP
#define PANZER_POINT_VALUES2_IMPL_HPP

#include "Intrepid2_CellTools.hpp"

#include "Panzer_CommonArrayFactories.hpp"

// ***********************************************************
// * Evaluation and SetupArrays are NOT specialized
// ***********************************************************

namespace panzer {

  

  template<typename Scalar,typename CoordinateArray>
  struct CopyNodeCoords {
    using size_type = typename panzer::PointValues2<Scalar>::size_type;
    const CoordinateArray source_;
    PHX::MDField<Scalar,Cell,NODE,Dim> target_;
    CopyNodeCoords(const CoordinateArray in_source,
                   PHX::MDField<Scalar,Cell,NODE,Dim> in_target)
      : source_(in_source),target_(in_target) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type cell,const size_type node,const size_type dim) const {
      target_(cell,node,dim) = source_(cell,node,dim);
    }
  };

  template<typename Scalar,typename CoordinateArray>
  struct CopyPointCoords {
    using size_type = typename panzer::PointValues2<Scalar>::size_type;
    const CoordinateArray source_;
    PHX::MDField<Scalar,IP,Dim> target_;
    CopyPointCoords(const CoordinateArray in_source,
                    PHX::MDField<Scalar,IP,Dim> in_target
                    )
      :source_(in_source),target_(in_target) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type pt,const size_type dim) const {
      target_(pt,dim) = source_(pt,dim);
    }
  };

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
  evaluateValues(const int in_num_cells)
  {
    if (point_rule->isSide()) {
       TEUCHOS_ASSERT(false); // not implemented!!!!
    }
    
    const int num_cells = in_num_cells < 0 ? (int) jac.extent(0) : in_num_cells;
    const auto cell_range = std::pair<int,int>(0,num_cells);
    auto s_jac = Kokkos::subview(jac.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    auto s_jac_inv = Kokkos::subview(jac_inv.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_node_coordinates = Kokkos::subview(node_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_point_coords = Kokkos::subview(point_coords.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    Intrepid2::CellTools<PHX::exec_space> cell_tools;

    cell_tools.setJacobian(s_jac, coords_ref.get_view(), s_node_coordinates, *(point_rule->topology));
    cell_tools.setJacobianInv(s_jac_inv, s_jac);
    cell_tools.setJacobianDet(s_jac_det, s_jac);
    
    // IP coordinates
    cell_tools.mapToPhysicalFrame(s_point_coords, coords_ref.get_view(), s_node_coordinates, *(point_rule->topology));
  }

  template <typename Scalar>
  template <typename CoordinateArray>
  void PointValues2<Scalar>::
  copyNodeCoords(const CoordinateArray& in_node_coords)
  {
    // copy cell node coordinates
    {
      const size_type num_cells = in_node_coords.extent(0);
      const size_type num_nodes = in_node_coords.extent(1);
      const size_type num_dims = in_node_coords.extent(2);

      Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<3>> policy({0,0,0},{num_cells,num_nodes,num_dims});
      Kokkos::parallel_for("PointValues2::copyNodeCoords",policy,panzer::CopyNodeCoords<Scalar,CoordinateArray>(in_node_coords,node_coordinates));
      PHX::Device::execution_space().fence();
    }
  }

  template <typename Scalar>
  template <typename CoordinateArray>
  void PointValues2<Scalar>::
  copyPointCoords(const CoordinateArray& in_point_coords)
  {
    // copy reference point values
    {
      const size_type num_points = in_point_coords.extent(0);
      const size_type num_dims = in_point_coords.extent(1);

      Kokkos::MDRangePolicy<PHX::Device::execution_space,Kokkos::Rank<2>> policy({0,0},{num_points,num_dims});
      Kokkos::parallel_for("PointValues2::copyPointCoords",policy,panzer::CopyPointCoords<Scalar,CoordinateArray>(in_point_coords,coords_ref));
      PHX::Device::execution_space().fence();
    }
  }

}

#endif
