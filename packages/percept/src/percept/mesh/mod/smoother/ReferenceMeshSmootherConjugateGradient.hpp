// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef ReferenceMeshSmootherConjugateGradient_hpp
#define ReferenceMeshSmootherConjugateGradient_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/gradient_functors.hpp>
#include <percept/mesh/mod/smoother/ReferenceMeshSmootherBase.hpp>
#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>
#include <percept/PerceptUtils.hpp>
#include <percept/mesh/mod/smoother/get_alpha_0_refmesh.hpp>
#include <percept/mesh/mod/smoother/get_edge_len_avg.hpp>

#define DEBUG_RMSCG_PRINT 0
#define RMSCG_PRINT(a) do { if (DEBUG_RMSCG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define RMSCG_PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define RMSCG_PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)

struct GenericAlgorithm_update_coordinates;

namespace percept {

struct SGridBoundarySelector : public StructuredGrid::MTSelector {
	PerceptMesh *m_eMesh;
	std::shared_ptr<BlockStructuredGrid> m_bsg;

	SGridBoundarySelector(PerceptMesh *eMesh) : m_eMesh(eMesh)
	{
	    m_bsg = m_eMesh->get_block_structured_grid();
	}

	bool operator()(StructuredCellIndex& index)
	{   //generalized boundary selector, should work with an arbitrary block structured mesh
		unsigned iblock = index[3];

		std::shared_ptr<StructuredBlock> sblock = m_bsg->m_sblocks[iblock];
		bool isBC = false;
		for(unsigned iBC = 0;iBC < sblock->m_boundaryConditions.size() ;iBC++)
		{
		    Ioss::BoundaryCondition bc = sblock->m_boundaryConditions[iBC];
		    std::array<unsigned, 3> beg;
		    std::array<unsigned, 3> end;

		    for(unsigned ijk=0;ijk<3;ijk++)
		    {
		        beg[ijk] = (unsigned)std::min(bc.m_rangeBeg[ijk],bc.m_rangeEnd[ijk]);
		        end[ijk] = (unsigned)std::max(bc.m_rangeBeg[ijk],bc.m_rangeEnd[ijk]);

		        beg[ijk]--;
		        end[ijk]--;
		    }
		    isBC = ( beg[0]<=index[0] && index[0]<=end[0] ) &&
                    ( beg[1]<=index[1] && index[1]<=end[1] ) &&
                    ( beg[2]<=index[2] && index[2]<=end[2] );
		    if(isBC)
		        break;
		}

		return isBC;
	}
};

    /// A Jacobian based optimization smoother, e.g. 1/A - 1/W, W/A - I, etc. (A = local current Jacobian, W is for original mesh)
    /// Conjugate-gradient version, element-based metrics

    template<typename MeshType>
    class ReferenceMeshSmootherConjugateGradientImpl : public ReferenceMeshSmootherBaseImpl<MeshType> {

    public:

      using Base = ReferenceMeshSmootherBaseImpl<MeshType>;
      using Base::m_eMesh;
      using Base::gradNorm;
      using Base::m_current_position;
      using Base::m_delta;
      using Base::m_weight;
      using Base::m_nweight;

      using Base::m_scale;

      using Base::m_dmax;
      using Base::m_dmax_relative;
      using Base::m_dnew;
      using Base::m_dold;
      using Base::m_d0;
      using Base::m_dmid;
      using Base::m_dd;
      using Base::m_alpha;
      using Base::m_alpha_0;
      using Base::m_grad_norm;
      using Base::m_grad_norm_scaled;
      using Base::m_total_metric;
      using Base::m_stage;
      using Base::m_iter;

      using Base::m_num_invalid;
      using Base::m_global_metric;
      using Base::m_untangled;
      using Base::m_num_nodes;

      using Base::m_coord_field_original;
      using Base::m_coord_field_projected;
      using Base::m_coord_field_current;
      using Base::m_coord_field_lagged;

      using Base::m_metric;

      using Base::m_use_ref_mesh;
      using Base::m_do_animation;
      using Base::m_meshGeometry;

      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length
      ReferenceMeshSmootherConjugateGradientImpl(PerceptMesh *eMesh,
//                                             typename MeshType::MTSelector *boundary_selector=0,
                                             STKMesh::MTSelector *stk_select=0,
                                             StructuredGrid::MTSelector *sgrid_select=0,
                                             typename MeshType::MTMeshGeometry *meshGeometry=0,
                                             int inner_iterations = 100,
                                             double grad_norm =1.e-8,
                                             int parallel_iterations = 20);


      virtual ~ReferenceMeshSmootherConjugateGradientImpl() {}

    public:
      double m_max_edge_length_factor;

      GenericAlgorithm_update_coordinates  <MeshType> m_coord_updater;
      GenericAlgorithm_total_element_metric<MeshType> m_metric_computinator;
      SGridGenericAlgorithm_total_element_metric<HexMeshSmootherMetric> sgrid_metric;
      sGrid_GenericAlgorithm_get_gradient_1 m_sgrid_gradient;
      SGrid_GenericAlgorithm_get_alpha_0 m_sgrid_get_alpha_0;
      sGrid_GenericAlgorithm_get_edge_lengths sgrid_gels;

    public:

      virtual void get_gradient();

      void debug_print(double alpha);

      virtual double run_one_iteration();

      virtual Double total_metric( Double alpha, double multiplicative_edge_scaling, bool& valid, size_t *num_invalid=0);
      virtual void update_node_positions( Double alpha);
      virtual bool check_convergence();

      Double line_search(bool& restarted, double mfac_mult=1.0);

      Double get_alpha_0();
      void get_edge_lengths(PerceptMesh * eMesh);


      void get_surface_normals(PerceptMesh * eMesh);
      void snap_nodes();

      double nodal_metric(typename MeshType::MTNode node, double alpha, double *coord_current, double *cg_d,  bool& valid );
      void nodal_gradient(typename MeshType::MTNode node, double alpha, double *coord_current, double *cg_d,  bool& valid, double *ng);

    public:
      // helpers
      virtual Double metric(typename MeshType::MTElement entity, bool& valid)
      {
        return Base::m_metric->metric(entity,valid);
      }
      Double nodal_edge_length_ave(typename MeshType::MTNode node);
    };

    //using ReferenceMeshSmootherConjugateGradient =  ReferenceMeshSmootherConjugateGradientImpl<STKMesh>;

  }//percept
//#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradientDef.hpp>
//#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradientSpec.hpp>

#endif
#endif
