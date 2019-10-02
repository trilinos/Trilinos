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
//#include <percept/structured/StructuredGridRefinerDef.hpp>
#include <percept/structured/StructuredBlock.hpp>
#include <Kokkos_Vector.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/structured/BlockStructuredGrid.hpp>


namespace percept {


template<typename MeshType>
class SmootherMetricImpl;
class PerceptMesh;

///////////////////////////////////////////////////////////////////////

    template<typename MetricType>
    struct SGridGenericAlgorithm_total_element_metric
    {//computes total metric of a structured mesh on a per block basis
      mutable MetricType m_metric;

      std::vector< Kokkos::View<int*, DataLayout , MemSpace > > element_invalid_flags_per_block;
      Kokkos::View<int*, DataLayout , MemSpace > flags_iterator;

      StructuredGrid::MTField *m_coord_field_current;
      StructuredGrid::MTField *m_coord_field_original;

      StructuredGrid::MTField::Array4D m_coords_current_iterator;
      StructuredGrid::MTField::Array4D m_coords_original_iterator;


      Double mtot;
      size_t n_invalid;
      bool valid;
      unsigned m_iBlock;

      std::vector<SGridSizes> block_sizes;
      std::vector<std::array<unsigned int,3>> loop_orderings;
      SGridSizes block_sizes_iterator;
      Kokkos::Array<unsigned int,3> loop_orderings_iterator;

      SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, Double mtot_in, size_t n_invalid_in);

      void run(unsigned iBlock);

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned& index, Double& mtot_loc) const;
    };

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
      Kokkos::View<int*, SecondaryDataLayout , SecondaryMemSpace > element_invalid_flags;
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

	struct element_invalid_flags_reducer
	{
    	Kokkos::View<int*, DataLayout , MemSpace > flags;
		unsigned n_invalid;

		void reduce()
		{
			if(!flags.size())
			{
				printf("flag view is empty, aborting reduction\n");
				return;
			}

			unsigned local = n_invalid;
			Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0,flags.size()),*this,local);
			n_invalid = local;
		}


	    KOKKOS_INLINE_FUNCTION
	    void operator()(const unsigned& index, unsigned& local_n_invalid) const
	    {
	      const_cast<element_invalid_flags_reducer *>(this)->operator()(index, local_n_invalid);
	    }

	    KOKKOS_INLINE_FUNCTION
	    void operator()(const unsigned& index, unsigned& local_n_invalid)
	    {
	    	local_n_invalid += flags(index);
	    }
	};

} // namespace percept

#endif

