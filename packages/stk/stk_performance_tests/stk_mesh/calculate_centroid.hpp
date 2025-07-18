// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef calculate_centroid_hpp
#define calculate_centroid_hpp

#include <vector>

#include "stk_mesh/base/NgpMesh.hpp"
#include "stk_mesh/base/NgpField.hpp"
#include "stk_mesh/base/GetNgpField.hpp"
#include "stk_mesh/base/GetNgpMesh.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/ForEachEntity.hpp"

namespace stk {
namespace performance_tests {

/**
   num_nodes: number of nodes connected to the element
   elem_node_coords: array of length num_nodes*3, containing the
   coordinates for an element's nodes
   elem_centroid: array of length 3
*/
template<class T>
inline
void
calculate_centroid_3d(
  size_t                num_nodes,
  const T *        elem_node_coords,
  T *              elem_centroid)
{
  //compute the element-centroid:
  for(size_t n = 0; n < num_nodes; ++n) {
    elem_centroid[0] += elem_node_coords[n*3 + 0];
    elem_centroid[1] += elem_node_coords[n*3 + 1];
    elem_centroid[2] += elem_node_coords[n*3 + 2];
  }
  elem_centroid[0] /= num_nodes;
  elem_centroid[1] /= num_nodes;
  elem_centroid[2] /= num_nodes;
}

/**
   num_elements: number of element
   num_nodes: number of nodes connected to the element
   elem_node_coords: array of length num_nodes*3, containing the
   coordinates for an element's nodes
   elem_centroid: array of length 3
*/
template<class T>
inline
void
calculate_centroid_3d(
  size_t                num_elements,
  size_t                num_nodes,
  const T *        elem_node_coords,
  T *              elem_centroid)
{
  //compute the element-centroid:
  for (size_t i = 0; i < num_elements; ++i) {
    for(size_t n = 0; n < num_nodes; ++n) {
      elem_centroid[i*3 + 0] += elem_node_coords[i*num_nodes*3 + n*3 + 0];
      elem_centroid[i*3 + 1] += elem_node_coords[i*num_nodes*3 + n*3 + 1];
      elem_centroid[i*3 + 2] += elem_node_coords[i*num_nodes*3 + n*3 + 2];
    }
    elem_centroid[i*3 + 0] /= num_nodes;
    elem_centroid[i*3 + 1] /= num_nodes;
    elem_centroid[i*3 + 2] /= num_nodes;
  }
}

template <typename CoordFieldType>
void calculate_centroid(const stk::mesh::NgpMesh &ngpMesh, const CoordFieldType &ngpCoords, const stk::mesh::Selector &sel, stk::mesh::NgpField<double> &ngpCentroid)
{
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, sel,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                 {
                                   stk::mesh::NgpMesh::ConnectedNodes nodes = ngpMesh.get_nodes(stk::topology::ELEM_RANK, elem);
                                   const unsigned numComponents = ngpCentroid.get_num_components_per_entity(elem);
                                   if (numComponents == 0) {
                                     return;
                                   }

                                   stk::mesh::EntityFieldData<double> elemCentroid = ngpCentroid(elem);
                                   elemCentroid[0] = 0.0;
                                   elemCentroid[1] = 0.0;
                                   elemCentroid[2] = 0.0;

                                   for (size_t i = 0; i < nodes.size(); i++) {
                                     stk::mesh::FastMeshIndex nodeIndex = ngpMesh.fast_mesh_index(nodes[i]);
                                     const stk::mesh::EntityFieldData<double> nodeCoords = ngpCoords(nodeIndex);
 
                                     elemCentroid[0] += nodeCoords[0];
                                     elemCentroid[1] += nodeCoords[1];
                                     elemCentroid[2] += nodeCoords[2];
                                   }

                                   elemCentroid[0] /= nodes.size();
                                   elemCentroid[1] /= nodes.size();
                                   elemCentroid[2] /= nodes.size();
                                 });
  ngpCentroid.modify_on_device();
}

template <typename CoordFieldType>
void calculate_centroid_using_coord_field(const stk::mesh::BulkData &bulk, const stk::mesh::Selector& selector, stk::mesh::FieldBase &centroid)
{
  const stk::mesh::FieldBase& coords = *bulk.mesh_meta_data().coordinate_field();
  stk::mesh::NgpField<double>& ngpCentroid = stk::mesh::get_updated_ngp_field<double>(centroid);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  stk::mesh::NgpField<double>& ngpCoords = stk::mesh::get_updated_ngp_field<double>(coords);
  calculate_centroid(ngpMesh, ngpCoords, selector, ngpCentroid);

  ngpCentroid.sync_to_host();
}

template <typename CoordFieldType>
void calculate_centroid_using_coord_field(const stk::mesh::BulkData &bulk, stk::mesh::FieldBase &centroid)
{
  calculate_centroid_using_coord_field<CoordFieldType>(bulk, bulk.mesh_meta_data().locally_owned_part(), centroid);
}

inline
void calc_centroid(const double** coordData, int numNodes, double* centroid)
{
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;

  for(int i=0; i<numNodes; ++i) {
    centroid[0] += coordData[i][0];
    centroid[1] += coordData[i][1];
    centroid[2] += coordData[i][2];
  }

  centroid[0] /= numNodes;
  centroid[1] /= numNodes;
  centroid[2] /= numNodes;
}

inline void calculate_centroid_using_host_coord_fields(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase& centroid)
{
  stk::mesh::Selector selector(centroid);
  selector &= bulk.mesh_meta_data().locally_owned_part();
  const stk::mesh::FieldBase& coords = *bulk.mesh_meta_data().coordinate_field();

  constexpr size_t maxNumNodes = 8;
  const double*elemNodeCoords[maxNumNodes];

  auto centroidCalculator = [&](stk::mesh::Entity elem, const stk::mesh::Entity* nodes, size_t numNodes)
  {
      STK_ThrowAssertMsg(numNodes <= maxNumNodes, "numNodes("<<numNodes<<") must be <= maxNumNodes("<<maxNumNodes<<")");

      for(size_t i = 0; i < numNodes; i++) {
        elemNodeCoords[i] = reinterpret_cast<double*>(stk::mesh::field_data(coords, nodes[i]));
      }

      double* centroidData = reinterpret_cast<double*>(stk::mesh::field_data(centroid, elem));
      calc_centroid(elemNodeCoords, numNodes, centroidData);
  };

  stk::mesh::for_each_entity_run_with_nodes(bulk, stk::topology::ELEM_RANK, selector, centroidCalculator);
}

template <typename CentroidFieldType>
inline
std::vector<double> get_centroid_average_from_host(stk::mesh::BulkData &bulk, CentroidFieldType& centroid,
                                                   const stk::mesh::Selector& selector)
{
  std::vector<double> average = {0, 0, 0};
  size_t numElems = 0;
  stk::mesh::Selector fieldSelector(centroid);
  fieldSelector &= selector;
  auto centroidData = centroid.template data<stk::mesh::ReadOnly>();

  for (const stk::mesh::Bucket *bucket : bulk.get_buckets(stk::topology::ELEM_RANK, fieldSelector)) {
    for (stk::mesh::Entity elem : *bucket) {
      auto centroidValues = centroidData.entity_values(elem);
      for (stk::mesh::ComponentIdx component : centroidValues.components()) {
        average[component] += centroidValues(component);
      }
      numElems++;
    }
  }

  if (numElems > 0) {
    for(size_t dim = 0; dim < 3; dim++) {
      average[dim] /= numElems;
    }
  }

  return average;
}

inline
std::vector<double> get_centroid_average_from_device(stk::mesh::BulkData &bulk, stk::mesh::Field<double> &centroid,
                                                     const stk::mesh::Selector& selector)
{
  stk::mesh::NgpField<double>& ngpField = stk::mesh::get_updated_ngp_field<double>(centroid);
  stk::mesh::NgpMesh& ngpMesh = stk::mesh::get_updated_ngp_mesh(bulk);

  typedef Kokkos::View<double*, stk::ngp::MemSpace> DeviceAverageView;
  typedef typename DeviceAverageView::HostMirror HostAverageView;

  DeviceAverageView deviceAverageView("averageVew", 4);
  HostAverageView hostAverageView = Kokkos::create_mirror_view(deviceAverageView);

  Kokkos::deep_copy(deviceAverageView, 0.0);

  stk::mesh::Selector fieldSelector(centroid);
  fieldSelector &= selector;
  stk::mesh::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, fieldSelector,
                                 KOKKOS_LAMBDA(const stk::mesh::FastMeshIndex& elem)
                                 {
                                   for(size_t dim = 0; dim < 3; dim++) {
                                     Kokkos::atomic_add(&deviceAverageView(dim), ngpField(elem, dim));
                                   }
                                   Kokkos::atomic_inc(&deviceAverageView(3));
                                 });

  Kokkos::deep_copy(hostAverageView, deviceAverageView);

  std::vector<double> average = {0, 0, 0};
  unsigned numElems = hostAverageView(3);

  for(size_t dim = 0; dim < 3; dim++) {
    average[dim] = numElems > 0 ? hostAverageView(dim) / numElems : 0.0;
  }

  return average;
}

} // namespace performance_tests
} // namespace stk

#endif
