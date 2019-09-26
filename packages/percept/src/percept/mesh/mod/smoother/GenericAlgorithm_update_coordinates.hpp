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

using Array4D = viewType;

template<typename MeshType>
struct A4DMD{
	unsigned blkIndex;
	unsigned numNodes;
	Kokkos::View</*typename MeshType::MTNode**/unsigned**, DataLayout, MemSpace> blockNodes;

	A4DMD(){
		blkIndex =0;
		numNodes=0;
	}

	~A4DMD(){}

};

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

    std::vector< Kokkos::View<Double*,DataLayout,MemSpace> > dmax_candidates; //only safely usable on on bsg's
    Kokkos::View<Double*,DataLayout,MemSpace> dmax_cands_this_block;
    Double dmax;

    std::vector< Kokkos::View<Double*,DataLayout,MemSpace> > dmax_relative_candidates; //only safely usable on on bsg's
    Kokkos::View<Double*,DataLayout,MemSpace> dmax_rel_cands_this_block;
    Double dmax_relative;

    Double /*double*/ alpha;
    bool m_calc_deltas;

    GenericAlgorithm_update_coordinates(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in);

    std::list<A4DMD<MeshType>> blockMetaDatas;
	Kokkos::View</*typename MeshType::MTNode**/unsigned**, DataLayout, MemSpace> nodesThisBlock; //will be used to loop over a particular block's nodes.
	Array4D cfc;
	Array4D cfl;
	Array4D cgsf;
	Array4D cgelf;

    KOKKOS_INLINE_FUNCTION
    void operator()(const int64_t& index) const;

    void run(bool calc_deltas=false);

  };

  struct simplified_gatm_1_BSG{ //uses A4D's directly

    using RefMeshSmoother =     ReferenceMeshSmootherConjugateGradientImpl<StructuredGrid>;
	  using This = simplified_gatm_1_BSG;
	  using Array4D = viewType;

	  RefMeshSmoother * m_rms; //adding smoother pointer a member didn't effect performance
	  PerceptMesh * m_eMesh;

    StructuredGrid::MTSelector on_locally_owned_part;
    StructuredGrid::MTSelector on_globally_shared_part;
    int spatialDim;
    MTSGridField * coord_field;
	  MTSGridField * coord_field_current;
	  MTSGridField * coord_field_lagged;
	  MTSGridField * cg_s_field;


	  std::vector<StructuredGrid::MTNode> nodes;

	  Double alpha;

	  std::list<A4DMD<StructuredGrid>> blockMetaDatas;
	  Kokkos::View<StructuredGrid::MTNode*, DataLayout, MemSpace> nodesThisBlock;

	  Array4D cfc;
	  Array4D cfl;
	  Array4D cgsf;



	  Array4D cfcB;
	  Array4D cflB;
	  Array4D cgsfB;
	  Kokkos::View< StructuredCellIndex*, DataLayout , MemSpace > nodesV;
	  unsigned numNodes;
	  unsigned numBlocks;

	  simplified_gatm_1_BSG(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in);

  	KOKKOS_INLINE_FUNCTION
  	void
  	operator()(int64_t index) const
  	{
  		int i  = nodesV(index)[0];
  		int j  = nodesV(index)[1];
  		int k  = nodesV(index)[2];

  		for (int iDim=0; iDim < spatialDim; iDim++)
  	        {
  			  Double dt = alpha*cgsfB(i,j,k,iDim);
  	  	  	  cfcB(i,j,k,iDim) += dt;
  	        }
  	}

	void run();
  };

  template<typename T>
  struct max_scanner
  {
      Kokkos::View<T*,DataLayout,MemSpace> candidates;

      KOKKOS_INLINE_FUNCTION
      void
      init(T&interimMax) const
      {//init function still not working correctly. Need to find out how to use it properly
          interimMax = candidates(0);
      }

      KOKKOS_INLINE_FUNCTION void
      join (volatile T& dst,
      const volatile T& src) const
      {
          if (dst < src) {
              dst = src;
          }
      }

      KOKKOS_INLINE_FUNCTION
      void
      operator()(const int64_t index, T&interimMax) const
      {
        if(interimMax<candidates(index))
            interimMax=candidates(index);
      }


      max_scanner(Kokkos::View<T*,DataLayout,MemSpace> candidates_in)
      {
          candidates=candidates_in;
      }

      T find_max()
      {
          T max;
          Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace>(0,candidates.size()),*this,max);
          return max;
      }

  };

  KOKKOS_INLINE_FUNCTION
  Double
  device_safe_max_H(Double a, Double b)
  {
      return (a>b ? a : b);
  }

} // namespace percept

#endif

