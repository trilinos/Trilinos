#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)


#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother2.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/JacobianUtil.hpp>
#include <stk_percept/math/Math.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>

#include "mpi.h"
#include <cstdio>

#define DEBUG_PRINT 0
#define PRINT(a) do { if (DEBUG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << " "; } while(0)

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;


    double PMMParallelReferenceMeshSmoother2::run_one_iteration( Mesh* mesh, MeshDomain *domain,
                                                              MsqError& err )
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();

      int spatialDim = m_eMesh->get_spatial_dim();

      stk::mesh::FieldBase *cg_g_field    = eMesh->get_field("cg_g");
      stk::mesh::FieldBase *cg_r_field    = eMesh->get_field("cg_r");
      stk::mesh::FieldBase *cg_d_field    = eMesh->get_field("cg_d");
      stk::mesh::FieldBase *cg_s_field    = eMesh->get_field("cg_s");

      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      bool total_valid=true;

      bool reduced_metric=false;

      m_dmax = 0.0;
      get_gradient(mesh, domain);

      {
        /// r = -g
        eMesh->nodal_field_axpby(-1.0, cg_g_field, 0.0, cg_r_field);
        /// s = r  (allows for preconditioning later s = M^-1 r)
        eMesh->copy_field(cg_s_field, cg_r_field);
        /// d = s
        eMesh->copy_field(cg_d_field, cg_s_field);
        /// dnew = r.d
        m_dnew = eMesh->nodal_field_dot(cg_r_field, cg_d_field);
      }

      double metric_orig = total_metric(mesh, 0.0, 1.0, total_valid);
      m_total_metric = metric_orig;
      if (check_convergence() || metric_orig == 0.0)
        {
          PRINT_1( "tmp srk already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " m_d0= " << m_d0 );
          //update_node_positions
          return total_metric(mesh,0.0,1.0, total_valid);
        }

      // node loop: local line search
      double alpha_min = 1.e+30;
      double alpha_max = 0.0;
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            // update local and globally shared 
            //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            if (on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    bool fixed = pmm->get_fixed_flag(&node);
                    bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                    if (fixed || isGhostNode)
                      {
                        continue;
                      }

                    //double edge_length_ave = m_eMesh->edge_length_ave(element);
                    double edge_length_ave = nodal_edge_length_ave(node);

                    double *coord_current = PerceptMesh::field_data(m_coord_field_current, node);
                    double *cg_d = PerceptMesh::field_data(cg_d_field, node);
                    double *cg_g = PerceptMesh::field_data(cg_g_field, node);

                    double local_scale = 0.0;
                    for (int i=0; i < spatialDim; i++)
                      {
                        local_scale = std::max(local_scale, std::abs(cg_g[i])/edge_length_ave);
                      }

                    local_scale = (local_scale < 1.0) ? 1.0 : 1.0/local_scale;
                    //PRINT("tmp srk node= " << node.identifier() << " iter= " << m_iter << " local_scale= " << local_scale);

                    /// line search

                    //get_nodal_gradient(node, cg_g_local);
                    double norm_gradient2 = 0.0;
                    for (int i=0; i < spatialDim; i++)
                      {
                        norm_gradient2 += cg_g[i]*cg_g[i];
                      }

                    double alpha = local_scale;  // m_scale
                    if (std::sqrt(norm_gradient2) > edge_length_ave*1.e-8)
                    {
                      bool local_valid=true, local_valid_0=true, local_valid_1=true;
                      double metric_0 = nodal_metric(node, 0.0, coord_current, cg_d, local_valid_0);
                      double metric=0.0;
                      //double sigma=0.95;
                      double tau = 0.5;
                      double c0 = 1.e-4;

                      double armijo_offset_factor = c0*norm_gradient2;
                      bool converged = false;
                      while (!converged)
                        {
                          metric = nodal_metric(node, alpha, coord_current, cg_d, local_valid);

                          //converged = (metric > sigma*metric_0) && (alpha > 1.e-16);
                          double mfac = alpha*armijo_offset_factor;
                          converged = metric == 0.0 || (metric < metric_0 + mfac);
                          if (m_untangled) converged = converged && local_valid;
                          PRINT(  "tmp srk node= " << node.identifier() << " iter= " << m_iter
                                  << " alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac) 
                                  << " m_untangled = " << m_untangled << " norm_gradient2= " << norm_gradient2
                                  << " local_valid= " << local_valid );
                          if (!converged)
                            alpha *= tau;
                          if (alpha < std::max(1.e-6*m_scale, 1.e-16))
                            {
                              PRINT_1(  "tmp srk not conv node= " << node.identifier() << " iter= " << m_iter
                                      << " alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac) 
                                      << " m_untangled = " << m_untangled << " norm_gradient2= " << norm_gradient2
                                      << " local_valid= " << local_valid    );

                              break;
                            }
                        }
                      //if (metric > sigma*metric_0)
                      if (!converged)
                        {
                          double metric_0 = nodal_metric(node, 0.0, coord_current, cg_d, local_valid_0);
                          double metric_1 = nodal_metric(node, 1.e-6, coord_current, cg_d, local_valid_1);
                          PRINT_1(  "tmp srk can't reduce metric, node= " << node.identifier() << " iter= " << m_iter
                                  << " alpha= " << alpha << " metric_0= " << metric_0 << " metric[1.e-6]= " << metric_1 << " diff[want neg]= " << metric_1 - metric_0
                                  << "\n m_untangled = " << m_untangled << " norm_gradient2= " << norm_gradient2
                                    << " local_valid= " << local_valid 
                                    << " local_valid_0= " << local_valid_0 
                                    << " local_valid_1= " << local_valid_1 
                                    << " cg_d= " << cg_d[0] << " " << cg_d[1] << " edge_length_ave= " << edge_length_ave);
                          alpha = 0.0;
                        }
                      else
                        {
                          reduced_metric = true;
                          double a1 = alpha/2.;
                          double a2 = alpha;
                          double f0 = metric_0, 
                            f1 = nodal_metric(node, a1, coord_current, cg_d, total_valid), 
                            f2 = nodal_metric(node, a2, coord_current, cg_d, total_valid);
                          double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
                          double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
                          if (std::fabs(den) > 1.e-10)
                            {
                              double alpha_quadratic = num/den;
                              if (alpha_quadratic < 2.0*alpha)
                                {
                                  double fm = nodal_metric(node, alpha_quadratic, coord_current, cg_d, total_valid);
                                  //if (fm < f2 && (!m_untangled || total_valid))
                                  if (fm < f2)
                                    {
                                      alpha = alpha_quadratic;
                                      PRINT( "tmp srk alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 );
                                    }
                                } 
                            }
                        }
                    }
                    alpha_min = std::min(alpha_min, alpha);
                    alpha_max = std::max(alpha_max, alpha);
                    // scale cg_d by local alpha
                    for (int i=0; i < spatialDim; i++)
                      {
                        cg_d[i] *= alpha;
                      }
                  } // node in bucket loop
              }
          } // bucket loop
      }

      /// check for valid mesh
      double alpha_global = 1.0;
      if (m_stage == 1)
      {
        bool converged = false;
        while (!converged)
          {
            int num_invalid=0;
            double metric_new = total_metric(mesh, alpha_global, 1.0, total_valid, &num_invalid);
            PRINT_1(  "tmp srk global alpha = " << alpha_global << " metric_new= " << metric_new << " total_valid= " << total_valid << " num_invalid= " << num_invalid);
            if (total_valid)
              {
                converged = true;
                break;
              }
            else
              {
                alpha_global *= 0.5;
              }
            if (alpha_global < 1.e-10)
              {
                break;
              }
          }
        //if (metric > sigma*metric_0)
        if (!converged)
          {
            bool local_valid_0=true, local_valid_1=true;
            double metric_0 = total_metric(mesh,0.0, 1.0, local_valid_0);
            double metric_1 = total_metric(mesh, 1.e-6, 1.0, local_valid_1);
            PRINT_1(  "tmp srk can't get valid mesh... " 
                      " metric_0= " << metric_0 << " metric[1.e-6]= " << metric_1 << " diff[want neg]= " << metric_1 - metric_0
                      << "\n m_untangled = " << m_untangled 
                      << " local_valid_0= " << local_valid_0 
                      << " local_valid_1= " << local_valid_1 );
          }
      }

      /// x = x + alpha*d
      m_dmax=0.0;
      update_node_positions(mesh, alpha_global);
      double metric_new = total_metric(mesh,0.0, 1.0, total_valid);
      PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha_min= " << alpha_min << " alpha_max= " << alpha_max << " total_valid= " << total_valid);
      if (m_stage == 1 && !total_valid)
        {
          throw std::runtime_error("created bad mesh");
        }

      //if (!reduced_metric || metric_new >= metric_orig) 
      if (!reduced_metric )
        {
          PRINT_1( "can't reduce metric, metric_new= " << metric_new << " metric_0 = " << metric_orig << " reduced_metric= " << reduced_metric);
          throw std::runtime_error("can't reduce metric");
        }

      return metric_new;
    }


  }
}


#endif

