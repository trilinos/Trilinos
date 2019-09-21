// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT) 


#include <percept/mesh/mod/smoother/ReferenceMeshSmoother3.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/Math.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>
#include <limits>

#include "mpi.h"
#include <cstdio>

#define DEBUG_PRINT 0
#define PRINT(a) do { if (DEBUG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << " "; } while(0)

  namespace percept {

    //static double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());

#if 0
    bool ReferenceMeshSmoother3::check_convergence()
    {
      if (m_stage == 0 && (m_dnew == 0.0 || m_total_metric == 0.0))
        {
          return true; // for untangle
        }
      if (m_stage == 0 && m_num_invalid == 0 && (m_dmax < gradNorm || m_grad_norm_scaled < gradNorm))
        {
          return true;
        }
      if (m_num_invalid == 0 && (m_grad_norm_scaled < gradNorm || (m_iter > 0 && m_dmax < gradNorm && m_dnew < gradNorm*gradNorm*m_d0)))
        {
          return true;
        }      
      return false;
    }

    /// gets a global scale factor so that local gradient*scale is approximately the size of the local mesh edges
    /// also uses the reference mesh to compute a local scaled gradient norm for convergence checks
    void ReferenceMeshSmoother3::get_scale()
    {
      PerceptMesh *eMesh = m_eMesh;
      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();

      m_scale = 1.e-10;

      // element loop
      {
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmoother::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];

                    const stk::mesh::PairIterRelation elem_nodes = element.relations( stk::topology::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    double edge_length_ave = m_eMesh->edge_length_ave(element, m_coord_field_original);

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        stk::mesh::Entity node = * elem_nodes[ inode ].entity();

                        bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                        VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
                        bool fixed = this->get_fixed_flag(&node);
                        if (fixed || isGhostNode)
                          continue;

                        double *cg_g = stk::mesh::field_data(cg_g_field, node);
                        
                        for (int idim=0; idim < spatialDim; idim++)
                          {
                            m_scale = std::max(m_scale, std::abs(cg_g[idim])/edge_length_ave);
                          }
                      }
                  }
              }
          }
      }

      {
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & m_scale ) );
        m_scale = (m_scale < 1.0) ? 1.0 : 1.0/m_scale;
        PRINT("tmp srk m_scale= " << m_scale);
      }

      // node loop
      m_grad_norm_scaled = 0.0;
      double gn=0.0;
      double el=1.e+10, es1=0, es2=0;
      //double pw=std::max(m_metric->length_scaling_power() - 1.0,0.0);

      {
        const stk::mesh::BucketVector & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity node = bucket[i_node];
                    bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                    VERIFY_OP_ON(isGhostNode, ==, false, "hmmmm");
                    bool fixed = this->get_fixed_flag(&node);
                    if (fixed || isGhostNode)
                      continue;

                    double edge_length_ave = nodal_edge_length_ave(node);
                    double *cg_g = stk::mesh::field_data(cg_g_field, node);
                        
                    double sum=0.0;
                    for (int idim=0; idim < spatialDim; idim++)
                      {
                        sum += cg_g[idim]*cg_g[idim];
                      }
                    sum = std::sqrt(sum);
                    gn = std::max(gn, sum);
                    double s1 = sum;
                    //sum = std::pow(sum, m_metric->length_scaling_power());
                    double s2 = sum;
                    sum /= edge_length_ave;
                    //sum /= (pw != 0 ? std::pow(edge_length_ave, pw) : edge_length_ave);
                    //m_grad_norm_scaled = std::max(m_grad_norm_scaled, sum);
                    if (sum > m_grad_norm_scaled)
                      {
                        m_grad_norm_scaled = sum;
                        el = edge_length_ave;
                        es1 = s1;
                        es2 = s2;
                      }
                  }
              }
          }
      }

      {
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & m_grad_norm_scaled ) );
        PRINT("tmp srk m_grad_norm_scaled= " << m_grad_norm_scaled << " gn= " << gn << " el= " << el << " es1= " << es1 << " es2=  " << es2);
      }

    }
#endif

    double ReferenceMeshSmoother3::run_one_iteration()
    {
      PerceptMesh *eMesh = m_eMesh;

      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_g");
      stk::mesh::FieldBase *cg_d_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_d");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      bool total_valid=true;

      if (1) //m_iter == 0)
        {
          m_dmax = 0.0;

          //PRINT_1("tmp srk get_gradient at m_iter=0");
          get_gradient(mesh, domain);

          // d = -g
          eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_d_field);

          m_dnew = eMesh->nodal_field_dot(cg_d_field, cg_d_field);
          PRINT("tmp srk m_dnew = " << m_dnew);

          // FIXME
          if (0)
            {
              double coord_mag = eMesh->nodal_field_dot(m_coord_field_current, m_coord_field_current);
              if (!m_eMesh->get_rank()) printf("tmp srk m_dnew[%d] = %30.10g m_dmid= %30.10g coord_mag= %30.10g\n",  m_iter, m_dnew, m_dmid, coord_mag);
            }

          /// d0 = dnew
          if (m_iter==0) m_d0 = m_dnew;
          m_grad_norm = std::sqrt(m_dnew);
        }

      double metric_check = total_metric(mesh, 0.0, 1.0, total_valid);
      m_total_metric = metric_check;
      if (check_convergence() || metric_check == 0.0)
        {
          PRINT_1( "tmp srk already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " m_d0= " << m_d0 << " m_grad_norm_scaled= " << m_grad_norm_scaled);
          //update_node_positions
          return total_metric(mesh,0.0,1.0, total_valid);
        }

      /// line search
      double cfl_factor = 0.2;
      double alpha = m_scale*cfl_factor;
      double metric_0 = total_metric(mesh, 0.0, 1.0, total_valid);
      double metric=0.0;
      metric = total_metric(mesh, alpha, 1.0, total_valid);
      //converged = (metric < metric_0 + mfac);
      //      if (m_untangled) converged = converged && total_valid;
      PRINT_1(  "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 )
                    << " m_untangled = " << m_untangled
                    << " total_valid= " << total_valid );

      /// x = x + alpha*d
      m_alpha = alpha;
      update_node_positions(mesh, alpha);
      //PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha= " << alpha);

      if (DEBUG_PRINT)
        {
          bool total_valid_0=true;
          total_metric(mesh, 0.0, 1.0, total_valid_0);
          if (m_stage != 0) VERIFY_OP_ON(total_valid_0, ==, true, "bad mesh after update_node_positions...");
        }

      double tm = total_metric(mesh,0.0,1.0, total_valid);

      return tm;
    }

  }

#endif
