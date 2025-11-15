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
#include "stk_util/parallel/ParallelReduce.hpp"

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

template <typename CentroidFieldType>
inline
std::vector<double> get_centroid_average_from_host(stk::mesh::BulkData &bulk, CentroidFieldType& centroid,
                                                   const stk::mesh::Selector& selector)
{
  std::vector<double> localAverage(4);
  stk::mesh::Selector fieldSelector(centroid);
  fieldSelector &= selector;
  auto centroidData = centroid.template data<>();

  for (const stk::mesh::Bucket *bucket : bulk.get_buckets(stk::topology::ELEM_RANK, fieldSelector)) {
    for (stk::mesh::Entity elem : *bucket) {
      auto centroidValues = centroidData.entity_values(elem);
      for (stk::mesh::ComponentIdx component : centroidValues.components()) {
        localAverage[component] += centroidValues(component);
      }
      localAverage[3] += 1;
    }
  }

  std::vector<double> globalAverage(4);
  stk::all_reduce_sum(bulk.parallel(), localAverage.data(), globalAverage.data(), 4);

  const unsigned numElems = globalAverage[3];

  std::vector<double> average(3);
  if (numElems > 0) {
    for (size_t dim = 0; dim < 3; ++dim) {
      average[dim] = globalAverage[dim] / numElems;
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
  typedef typename DeviceAverageView::host_mirror_type HostAverageView;

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

  std::vector<double> globalAverage(4);
  stk::all_reduce_sum(bulk.parallel(), hostAverageView.data(), globalAverage.data(), 4);

  unsigned numElems = globalAverage[3];

  std::vector<double> average(3);
  for(size_t dim = 0; dim < 3; ++dim) {
    average[dim] = numElems > 0 ? globalAverage[dim] / numElems : 0.0;
  }

  return average;
}

} // namespace performance_tests
} // namespace stk

#endif
