// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/geometry/kernel/GeometryKernel.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherBase.hpp>
#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <stdio.h>
#include <ctime>

#include "mpi.h"

  namespace percept {

  template <typename MeshType>
  ReferenceMeshSmootherBaseImpl<MeshType>::
  ReferenceMeshSmootherBaseImpl(PerceptMesh *eMesh,
                            STKMesh::MTSelector *stk_select,
                            typename MeshType::MTMeshGeometry *meshGeometry,
                            int inner_iterations,
                            double grad_norm,
                            int parallel_iterations)
    : MeshSmootherImpl<MeshType>(eMesh, stk_select,meshGeometry, inner_iterations, grad_norm, parallel_iterations),
      m_scale(0),
      m_dmax(0),
      m_dmax_relative(0),
      m_dnew(0), m_dold(0), m_d0(1), m_dmid(0), m_dd(0), m_alpha(0), m_alpha_0(0), m_grad_norm(0), m_grad_norm_scaled(0),
      m_total_metric(0),
      m_stage(0),
      m_iter(0),
      m_num_invalid(0), m_global_metric(std::numeric_limits<double>::max()), m_untangled(false), m_num_nodes(0),
      m_use_ref_mesh(true), m_do_animation(0)
  {
    if (eMesh->getProperty("ReferenceMeshSmootherBase_do_anim") != "")
      {
        m_do_animation = toInt(eMesh->getProperty("ReferenceMeshSmootherBase_do_anim"));
        std::cout << "m_do_animation= " << m_do_animation << std::endl;
      }


    std::string fileName = "smootherlog" + eMesh->getProperty("in_filename") + ".txt";

    if(Base::m_eMesh->get_rank() == 0) {

      myFile.open(fileName.c_str(), std::ios::out);

      myFile.precision(6);
      myFile.setf(std::ios_base::left,std::ios_base::adjustfield);
      myFile.setf(std::ios_base::scientific);

      unsigned width = 14;
      myFile << std::setw(8) << "%iter";
      myFile << std::setw(width) << "global_metric";
      myFile << std::setw(width) << "invalid_elems";
      myFile << std::setw(width) << "dmax";
      myFile << std::setw(width) << "dmax_rel";
      myFile << std::setw(width) << "grad/g0";
      myFile << std::setw(width) << "gradScaled";
      myFile << std::setw(width) << "alpha";
      myFile << std::setw(width) << "alpha_0";
      myFile << std::setw(6) << "stage";
      myFile << std::setw(width) << "untangled_elems";
      myFile << std::setw(width) << "noNodes";
      myFile << std::endl;
    }
  }

  template <typename MeshType>
    ReferenceMeshSmootherBaseImpl<MeshType>::
	~ReferenceMeshSmootherBaseImpl()
  {}

  template<>
    void ReferenceMeshSmootherBaseImpl<STKMesh>::sync_fields(int /*iter*/)
    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(m_eMesh->get_coordinates_field());
      fields.push_back(m_coord_field_projected);
      fields.push_back(m_coord_field_lagged);
      fields.push_back(m_coord_field_current);
      stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    }

    template<>
    void ReferenceMeshSmootherBaseImpl<STKMesh>::check_project_all_delta_to_tangent_plane(stk::mesh::FieldBase *field)
    {
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );

      int nc = m_eMesh->get_spatial_dim();

      // node loop
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_nodes_in_bucket = bucket.size();

              for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                {
                  stk::mesh::Entity node = bucket[i_node];
                  std::pair<bool,int> fixed = this->get_fixed_flag(node);
                  if (fixed.first)
                    {
                      continue;
                    }

                  double dproj[3];
                  double *cg_g = m_eMesh->field_data(field, node);
                  for (int ii=0; ii < nc; ++ii)
                    dproj[ii] = cg_g[ii];

                  if (fixed.second == MS_SURFACE)
                    {
                      double normal[3];
                      project_delta_to_tangent_plane(node, dproj, normal);
                      double dot=0.0, sc=0.0;
                      for (int ii=0; ii < nc; ++ii)
                        {
                          dot += dproj[ii]*normal[ii];
                          sc += dproj[ii]*dproj[ii];
                        }
                      sc = std::sqrt(sc);
                      if (sc > 1.e-10)
                        VERIFY_OP_ON(std::fabs(dot), <=, 1.e-6*sc, "bad norm");
                    }
                }
            }
        }
    }

    template<typename MeshType>
    bool ReferenceMeshSmootherBaseImpl<MeshType>::check_convergence()
    {
      throw std::runtime_error("not implemented");
      return false;
    }

    template<typename MeshType>
    double ReferenceMeshSmootherBaseImpl<MeshType>::run_one_iteration()
    {
      throw std::runtime_error("not implemented");
      return 0.0;
    }

    template<typename MeshType>
    void anim(PerceptMesh *eMesh, int& anim_step)
    {
    }

    template<>
    void anim<STKMesh>(PerceptMesh *eMesh, int& anim_step)
    {
      std::ostringstream fileid_ss;
      fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

      std::string oname = "anim_all.e";
      if (anim_step > 0) oname += "-s" + fileid_ss.str();
      eMesh->save_as(oname);
      ++anim_step;
    }

    template<typename MeshType>
    void set_geometry_info(typename MeshType::MTMeshGeometry *meshGeom)
    {
    }

    template<>
    void set_geometry_info<STKMesh>(typename STKMesh::MTMeshGeometry *meshGeometry)
    {
#if defined(STK_PERCEPT_HAS_GEOMETRY)
      if (meshGeometry && meshGeometry->geomKernel)
        {
          meshGeometry->geomKernel->m_useFoundFaceMap = true;
        }
#endif
    }

    template<typename MeshType>
    void ReferenceMeshSmootherBaseImpl<MeshType>::print_iteration_status(const int num_invalid_0)
	{
    	if (0 != Base::m_eMesh->get_rank()) return;

    	// m_iter=0 write header

    	// write columns
        unsigned width = 14;
        myFile << std::setw(8) << m_iter;
        myFile << std::setw(width) << m_global_metric;
        myFile << std::setw(width) << num_invalid_0;
        myFile << std::setw(width) << m_dmax;
        myFile << std::setw(width) << m_dmax_relative;
        myFile << std::setw(width) << std::sqrt(m_dnew/m_d0);
        myFile << std::setw(width) << m_grad_norm_scaled;
        myFile << std::setw(width) << m_alpha;
        myFile << std::setw(width) << m_alpha_0;
        myFile << std::setw(6) << m_stage;
        myFile << std::setw(width) << m_untangled;
        myFile << std::setw(width) << m_num_nodes;
        myFile << std::endl;
	}

    template<typename MeshType>
    void ReferenceMeshSmootherBaseImpl<MeshType>::run_algorithm()
    {
      PerceptMesh *eMesh = Base::m_eMesh;
      set_local_field_ptrs(eMesh);

      m_num_nodes = eMesh->get_number_nodes();

      eMesh->copy_field("coordinates_lagged", "coordinates_NM1");

      m_anim_step = 0;

      // untangle
      SmootherMetricUntangleImpl<MeshType> untangle_metric(eMesh);

      // shape
      SmootherMetricShapeB1Impl<MeshType> shape_b1_metric(eMesh, m_use_ref_mesh);

        {
          eMesh->nodal_field_axpbypgz(1.0, "coordinates_N", 0.0, "coordinates_NM1", 0.0, "coordinates");
          size_t num_invalid = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);
          if (!eMesh->get_rank())
            std::cout << "MeshSmoother: number of invalid elements before untangling = " << num_invalid << std::endl;

          m_num_invalid = num_invalid;
          m_untangled = (m_num_invalid == 0);

          int iter_all=0;

          int do_anim = m_do_animation; // = frequency of anim writes
          if (do_anim)
            {
              anim<MeshType>(eMesh, m_anim_step);
            }

          int nstage = get_num_stages();
          for (int stage = 0; stage < nstage; stage++)
            {
              set_geometry_info<MeshType>(Base::m_meshGeometry);

              m_stage = stage;
              if (stage==0)
                {
                  m_metric = &untangle_metric;
                }
              else
                {
                  int num_invalid_1 = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);
                  VERIFY_OP_ON(num_invalid_1, ==, 0, "Invalid elements exist for start of stage 2, quiting");

                  if (!eMesh->get_rank())
                    std::cout << "MeshSmoother: number of invalid elements after  untangling = " << num_invalid_1 << std::endl;

                  m_metric = &shape_b1_metric;
                }

              for (int iter = 0; iter < Base::innerIter; ++iter, ++iter_all)
                {
                  m_iter = iter;
                  m_num_invalid = MeshSmootherImpl<MeshType>::parallel_count_invalid_elements(eMesh);
                  int num_invalid_0 = m_num_invalid;
                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }

                  m_global_metric = run_one_iteration();

                  sync_fields(iter);
                  bool conv = check_convergence();

                  print_iteration_status(num_invalid_0);

                  if (do_anim)
                    {
                      if (iter_all % do_anim == 0)
                        {
                          anim<MeshType>(eMesh, m_anim_step);
                        }
                    }

                  if (!m_untangled && m_num_invalid == 0)
                    {
                      m_untangled = true;
                    }
                  if (conv && m_untangled)
                    {
                      break;
                    }
                }
            }

          eMesh->copy_field("coordinates_lagged", "coordinates");
        }
    }

    template<>
	void ReferenceMeshSmootherBaseImpl<STKMesh>::set_local_field_ptrs(PerceptMesh *eMesh)
	{
    	m_coord_field_original  = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_NM1");
    	m_coord_field_projected = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_N");
    	m_coord_field_lagged    = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");
    	m_coord_field_current   = eMesh->get_coordinates_field();
	}

    template class ReferenceMeshSmootherBaseImpl<STKMesh>;
  }

#endif
