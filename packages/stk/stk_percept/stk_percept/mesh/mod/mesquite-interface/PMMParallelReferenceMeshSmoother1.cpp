#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/JacobianUtil.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;

    double PMMParallelReferenceMeshSmoother1::metric(stk::mesh::Entity& element)
    {
      JacobianUtil jacA, jacW;
      double A_ = 0.0, W_ = 0.0; // current and reference detJ
      jacA(A_, *m_eMesh, element, m_coord_field_current);
      jacW(W_, *m_eMesh, element, m_coord_field_original);
      double val=0.0;
      for (int i=0; i < jacA.m_num_nodes; i++)
        {
          double Ai = jacA.mMetrics[i];
          double Wi = jacW.mMetrics[i];
          Ai = (Ai < 0 ? -1 : 1)*(std::abs(Ai) + 1.e-10);
          Wi = (Wi < 0 ? -1 : 1)*(std::abs(Wi) + 1.e-10);

          double v = (1.0/Ai - 1.0/Wi);
          val += v*v;
        }
      val = std::sqrt(val);
      return val;
    }

    void PMMParallelReferenceMeshSmoother1::run_one_iteration( Mesh* mesh, MeshDomain *domain,
                                                              MsqError& err )
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();

      // node loop, initialize delta to 0, etc, for all nodes including non-locally owned
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (on_locally_owned_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    m_current_position[&node] = Vector(coord_current, coord_current+spatialDim);
                    m_delta[&node] = Vector(spatialDim, 0.0);
                    m_weight[&node] = Vector(spatialDim, 0.0);
                    m_nweight[&node] = Vector(spatialDim, 0.0);
                  }
              }
          }
      }
      
      // element loop: compute deltas
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (on_locally_owned_part(**k))  
            // loop over all elements
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    double edge_length_ave = m_eMesh->edge_length_ave(element);

                    double metric_0 = metric(element);

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = * elem_nodes[ inode ].entity();

                        bool fixed = pmm->get_fixed_flag(&node);
                        if (fixed)
                          continue;

                        double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                        
                        Vector& delta = m_delta[&node];
                        Vector& weight = m_weight[&node];
                        Vector& nweight = m_nweight[&node];
                        
                        double eps = 1.e-5;
                        double eps1 = eps*edge_length_ave;

                        for (int i=0; i < spatialDim; i++)
                          {
                            coord_current[i] += eps1;
                            double mp = metric(element);
                            coord_current[i] -= 2.0*eps1;
                            double mm = metric(element);
                            coord_current[i] += eps1;
                            double dd = (mp - mm)/(2*eps1);
                            delta[i] += dd;
                            weight[i] += edge_length_ave;
                            nweight[i] += 1.0;
                            
                            {
                              coord_current[i] -= dd*eps;
                              double m1=metric(element);
                              coord_current[i] += dd*eps;
                              if (metric_0 > 1.e-6)
                                {
                                  if (!(m1 < metric_0))
                                    {
                                      std::cout << "bad grad" << " m1-metric_0 = " << (m1-metric_0) << std::endl;
                                    }
                                  VERIFY_OP_ON(m1, <, metric_0, "bad gradient");
                                }
                            }
                          }

                      }
                  }
              }
          }
      }

      // node loop: get scaling factor
      //double scale = std::numeric_limits<double>::max();
      double scale = 1.e-10; // give it a lower limit
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    Vector& delta = m_delta[&node];
                    Vector& weight = m_weight[&node];
                    Vector& nweight = m_nweight[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = std::abs(delta[i] / (weight[i]/nweight[i]));
                        scale = std::max(scale, dt);
                      }
                  }
              }
          }
        std::cout << "tmp srk scale= " << scale;
        scale = m_max_edge_length_factor/scale;
        std::cout << " scale2= " << scale << std::endl;

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        // negate since we want to reduce the metric
                        delta[i] = -delta[i] * scale;
                      }
                  }
              }
          }
      }

      /// line search
      double alpha = 1.0;
      {
        double metric_0 = total_metric(mesh, 0.0, 1.0);
        double metric=0.0;
        double sigma=0.95;
        while (((metric = total_metric(mesh, alpha, 1.0)) > sigma*metric_0) && (alpha > 1.e-16))
          {
            std::cout << "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << std::endl;
            alpha = alpha / 2.0;
          }
        if (metric > sigma*metric_0)
          {
            throw std::runtime_error("can't reduce metric");
          }
      }

      m_dmax = 0.0;
      // node loop: update node positions
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            // update local and globally shared 
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    if (fixed)
                      {
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    //m_current_position[&node] = Vector(coord_current, coord_current+spatialDim);
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha*delta[i];
                        m_dmax = std::max(std::abs(dt), m_dmax);
                        coord_current[i] += dt;  
                      }
                  }
              }
          }
      }

      MSQ_ERRRTN(err);
    }

    double PMMParallelReferenceMeshSmoother1::total_metric(Mesh *mesh, double alpha, double multiplicative_edge_scaling)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      stk::mesh::FieldBase *coord_field = m_eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;
      stk::mesh::FieldBase *coord_field_lagged  = m_eMesh->get_field("coordinates_lagged");
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = m_eMesh->get_spatial_dim();

      double mtot = 0.0;

      // cache coordinates
      m_eMesh->copy_field(coord_field_lagged, coord_field_current);

      // node loop
      {
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            // update local and globally shared 
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    if (fixed)
                      {
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha * multiplicative_edge_scaling * delta[i];
                        coord_current[i] += dt;
                      }
                  }
              }
          }
      }

      // element loop
      {
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];
                    mtot += metric(element);
                  }
              }
          }
      }

      // reset coordinates
      m_eMesh->copy_field(coord_field_current, coord_field_lagged);

      stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceSum<1>( & mtot ) );

      return mtot;
      
    }
  }
}


#endif
#endif
