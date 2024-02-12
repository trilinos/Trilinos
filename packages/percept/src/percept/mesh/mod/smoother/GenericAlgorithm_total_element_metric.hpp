// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef GenericAlgorithm_total_element_metric_hpp
#define GenericAlgorithm_total_element_metric_hpp

#include <percept/MeshType.hpp>

#include <percept/PerceptMesh.hpp>

namespace percept {


template<typename MeshType>
class SmootherMetricImpl;
class PerceptMesh;

///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////

    template<typename MeshType>
    struct GenericAlgorithm_total_element_metric
    {
      using This = GenericAlgorithm_total_element_metric<MeshType>;

      SmootherMetricImpl<MeshType> *m_metric;

      PerceptMesh *m_eMesh;

      typename MeshType::MTSelector on_locally_owned_part;
      typename MeshType::MTSelector on_globally_shared_part;
      int spatialDim;
      typename MeshType::MTField *coord_field;
      typename MeshType::MTField *coord_field_current;
      typename MeshType::MTField *coord_field_lagged;

      typename MeshType::MTField *cg_s_field;

      std::vector<typename MeshType::MTElement> elements;
      std::vector<int> element_invalid_flags;
      std::vector<const typename MeshType::MTCellTopology *> topos;
      bool valid;
      size_t * num_invalid;
      Double mtot;
      size_t n_invalid;

      GenericAlgorithm_total_element_metric(SmootherMetricImpl<MeshType> *metric, PerceptMesh *eMesh, bool valid_in, size_t *num_invalid_in, Double mtot_in, size_t n_invalid_in);

      void run(unsigned iBlock);

      void updateState(SmootherMetricImpl<MeshType> *metric, PerceptMesh *eMesh, bool valid_in, size_t *num_invalid_in, Double mtot_in, size_t n_invalid_in)
      {m_metric = metric; m_eMesh = eMesh; valid = valid_in; num_invalid = num_invalid_in; mtot = mtot_in; n_invalid = n_invalid_in;}

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, Double& mtot_loc) const
      {
        const_cast<This *>(this)->operator()(index, mtot_loc);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, Double& mtot_loc);

    };
} // namespace percept

#endif

