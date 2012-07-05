#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;

    void PMMParallelReferenceMeshSmoother::run_one_iteration( Mesh* mesh, MeshDomain *domain,
                                                              MsqError& err )
    {
      //std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk PMMParallelReferenceMeshSmoother::run_one_iteration start..." << std::endl;

      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field->field_state(stk::mesh::StateNP1);
      stk::mesh::FieldBase *coord_field_projected = coord_field->field_state(stk::mesh::StateN);
      stk::mesh::FieldBase *coord_field_original  = coord_field->field_state(stk::mesh::StateNM1);

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      int spatialDim = eMesh->get_spatial_dim();

      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))  // this is where we do part selection
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    m_current_position[&node] = Vector(coord_current, coord_current+spatialDim);
                    m_delta[&node] = Vector(spatialDim,0.0);
                  }
              }
          }
      }
      
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        Vector centroid_current(spatialDim, 0.0);
        Vector centroid_projected(spatialDim, 0.0);
        Vector centroid_original(spatialDim, 0.0);

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))  // this is where we do part selection
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();
                    centroid_current.assign(spatialDim, 0.0);
                    centroid_projected.assign(spatialDim, 0.0);
                    centroid_original.assign(spatialDim, 0.0);
                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = * elem_nodes[ inode ].entity();
                        double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                        double *coord_projected = PerceptMesh::field_data(coord_field_projected, node);
                        double *coord_original = PerceptMesh::field_data(coord_field_original, node);

                        for (int i=0; i < spatialDim; i++)
                          {
                            centroid_current[i] += coord_current[i]/((double)num_node);
                            centroid_projected[i] += coord_projected[i]/((double)num_node);
                            centroid_original[i] += coord_original[i]/((double)num_node);
                          }
                      }
                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = * elem_nodes[ inode ].entity();
                        double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                        double *coord_projected = PerceptMesh::field_data(coord_field_projected, node);
                        double *coord_original = PerceptMesh::field_data(coord_field_original, node);
                        Vector& delta = m_delta[&node];
                        for (int i=0; i < spatialDim; i++)
                          {
                            double alpha_prev = 0.0;
                            //double alpha_prev = m_alpha_prev;
                            double coord_base = coord_original[i]*(1.0-alpha_prev) + alpha_prev*coord_projected[i];
                            double centroid_base = centroid_original[i]*(1.0-alpha_prev) + alpha_prev*centroid_projected[i];
                            double new_pos = coord_base + (centroid_current[i] - centroid_base);
                            delta[i] += (new_pos - coord_current[i])/((double)num_node);
                          }
                      }
                  }
              }
          }
      }

      m_dmax = 0.0;
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))  // this is where we do part selection
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    if (fixed)
                      {
                        //std::cout << "tmp srk found fixed= " << node << std::endl;
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    //m_current_position[&node] = Vector(coord_current, coord_current+spatialDim);
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        m_dmax = std::max(std::abs(delta[i]), m_dmax);
                        coord_current[i] += delta[i];  // if not fixed
                      }
                  }
              }
          }
      }

      //if (!get_parallel_rank()) 
      //std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelReferenceMeshSmoother: running shape improver... done \n" << std::endl;

      MSQ_ERRRTN(err);
    }

    void PMMParallelReferenceMeshSmoother::run_wrapper( Mesh* mesh,
                                                        ParallelMesh* pmesh,
                                                        MeshDomain* domain,
                                                        Settings* settings,
                                                        QualityAssessor* qa,
                                                        MsqError& err )
    {
      std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk PMMParallelReferenceMeshSmoother innerIter= " << innerIter << " parallelIterations= " << parallelIterations << std::endl;

      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelReferenceMeshSmoother: running shape improver... \n" << std::endl;

      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field->field_state(stk::mesh::StateNP1);
      stk::mesh::FieldBase *coord_field_projected = coord_field->field_state(stk::mesh::StateN);
      stk::mesh::FieldBase *coord_field_original  = coord_field->field_state(stk::mesh::StateNM1);

      double alphas[] = {0.0,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0};
      int nalpha = sizeof(alphas)/sizeof(alphas[0]);
      
      for (int outer = 0; outer < nalpha; outer++)
        {
          double alpha = alphas[outer];
          m_alpha = alpha;
          m_alpha_prev = alpha;
          if (outer > 0) m_alpha_prev = alphas[outer-1];

          // set current state and evaluate mesh validity (current = alpha*project + (1-alpha)*original)
          eMesh->nodal_field_axpbypgz(alpha, coord_field_projected, (1.0-alpha), coord_field_original, 0.0, coord_field_current);

          int num_invalid = PMMParallelShapeImprover::count_invalid_elements(*mesh, domain);
          
          if (!get_parallel_rank()) 
            std::cout << "\ntmp srk PMMParallelReferenceMeshSmoother num_invalid current= " << num_invalid << " for outer_iter= " << outer
                      << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : " OK")
                      << std::endl;

          for (int iter = 0; iter < innerIter; iter++)
            {
              //
              int num_invalid = PMMParallelShapeImprover::count_invalid_elements(*mesh, domain);
//               if (!get_parallel_rank() && num_invalid) 
//                 std::cout << "\ntmp srk PMMParallelReferenceMeshSmoother num_invalid current= " << num_invalid 
//                           << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : "OK")
//                           << std::endl;
              run_one_iteration(mesh, domain, err);
              std::cout << "tmp srk iter= " << iter << " dmax= " << m_dmax << " num_invalid= " << num_invalid << std::endl;
              //sync_fields();
              //check_convergence();
              if (m_dmax < gradNorm)
                break;
            }
        }

      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelReferenceMeshSmoother: running shape improver... done \n" << std::endl;

      MSQ_ERRRTN(err);
    }


  }
}


#endif
#endif
