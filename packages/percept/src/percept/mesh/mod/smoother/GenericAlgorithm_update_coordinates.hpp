// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GenericAlgorithm_update_coordinates_hpp
#define GenericAlgorithm_update_coordinates_hpp

#include <percept/MeshType.hpp>

namespace percept {

template<typename MeshType>
class ReferenceMeshSmootherConjugateGradientImpl;

  template<typename MeshType>
  struct GenericAlgorithm_update_coordinates
  {
    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<MeshType>;
    using This = GenericAlgorithm_update_coordinates<MeshType>;

    RefMeshSmoother *m_rms;
    PerceptMesh *m_eMesh;

    typename MeshType::MTSelector on_locally_owned_part;
    typename MeshType::MTSelector on_globally_shared_part;
    int spatialDim;
    typename MeshType::MTField *coord_field;
    typename MeshType::MTField *coord_field_current;
    typename MeshType::MTField *coord_field_lagged;
    typename MeshType::MTField *cg_s_field;
    typename MeshType::MTField *cg_edge_length_field;


    std::vector<typename MeshType::MTNode> nodes;

    Double dmax;

    Double dmax_relative;

    Double /*double*/ alpha;
    bool m_calc_deltas;

    GenericAlgorithm_update_coordinates(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in);

    KOKKOS_INLINE_FUNCTION
    void operator()(const int64_t& index) const;

    void run(bool calc_deltas=false);

  };

  KOKKOS_INLINE_FUNCTION
  Double
  device_safe_max_H(Double a, Double b)
  {
      return (a>b ? a : b);
  }

} // namespace percept

#endif

