#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/JacobianUtil.hpp>
#include <stk_percept/math/Math.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;
    const bool do_tot_test = false;
    bool do_print_elem_val = false;

    double PMMParallelReferenceMeshSmoother1::metric(stk::mesh::Entity& element)
    {
      JacobianUtil jacA, jacW;
      //jacA.m_scale_to_unit = true;

      double A_ = 0.0, W_ = 0.0; // current and reference detJ
      jacA(A_, *m_eMesh, element, m_coord_field_current);
      jacW(W_, *m_eMesh, element, m_coord_field_original);
      double val=0.0, val_shape=0.0, val_untangle=0.0;
      double A_tot=0, W_tot=0;
      MsqMatrix<3,3> ident; 
      ident.identity();

      for (int i=0; i < jacA.m_num_nodes; i++)
        {
          double Ai = jacA.mMetrics[i];
          double Wi = jacW.mMetrics[i];
          if (Ai < 0)
            {
              //std::cout << "Ai= " << Ai << std::endl;
            }
          //VERIFY_OP_ON(Ai, >, 0, "Ai < 0");
          //VERIFY_OP_ON(Wi, >, 0, "Wi < 0");
          //Ai = (Ai < 0 ? -1 : 1)*(std::fabs(Ai) + 1.e-10);
          //Wi = (Wi < 0 ? -1 : 1)*(std::fabs(Wi) + 1.e-10);

          A_tot += Ai;
          W_tot += Wi;
          //double v = (1.0/Ai - 1.0/Wi);
          //double v = (Ai - Wi);
          //val += v*v;
          //val += v;
          double shape_metric = 0.0;
          if (std::fabs(Ai) > 1.e-10)
            {
              shape_metric = sqr_Frobenius(jacW.mJ[i]*inverse(jacA.mJ[i]) - ident);
              //shape_metric = sqr_Frobenius(jacW.mJ[i]-jacA.mJ[i]);
            }
          //shape_metric = sqr_Frobenius(jacW.mJ[i]-jacA.mJ[i]);
          //shape_metric = (1.0/Ai - 1.0/Wi)*(1.0/Ai - 1.0/Wi);
          //shape_metric = (Ai - Wi);
          //shape_metric = shape_metric*shape_metric;
          double metric_switch = Math::heavy_smooth(-Ai,0.0,1.e-6);
          double untangle_metric = -Ai*metric_switch;
          //untangle_metric = (Ai < 0 ? -100.0*Ai : Ai*Ai*Ai/100.0);
          //untangle_metric = (Ai < 0 ? -Ai : Wi/Ai);
          //untangle_metric = Ai*Ai/(Wi*Wi);
          double beta = 0.05*Wi;
          double temp_var = Ai;
          temp_var -= beta;
          double fval=0.0;
          if(temp_var<0.0){
            fval=std::fabs(temp_var)-temp_var;
          }
          untangle_metric = fval*fval;

          //untangle_metric = (Ai-1)*(Ai-1);
          //untangle_metric = (Ai < 0 ? std::fabs(Ai*Ai*Ai)/(Wi*Wi*Wi) : 0);
          //untangle_metric = (Ai < 0 ? (Ai*Ai)/(Wi*Wi) : 0);
          //untangle_metric = Ai*Ai/(Wi*Wi);
          //untangle_metric = shape_metric;
          //untangle_metric = -Ai;
          //val += untangle_metric + shape_metric * (1.0 - metric_switch);
          //val_shape += shape_metric*(1.0 - metric_switch);
          val_shape += shape_metric;
          val_untangle += untangle_metric;
        }
      //val = std::sqrt(val);
      //val = (A_tot - W_tot);
      //val = std::fabs(A_tot - W_tot);
      //val = (1.0/A_tot - 1.0/W_tot);
      //val = val*val;
      //val = Math::my_max_hi(val_untangle, val_shape, 1.e-6);
      //val = m_num_invalid ? val_untangle : val_shape;
      val = !m_untangled ? val_untangle : val_shape;
      //val = val_shape;
      //val = val_untangle;
      return val;
    }

    double PMMParallelReferenceMeshSmoother1::run_one_iteration( Mesh* mesh, MeshDomain *domain,
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

        int num_fixed=0;
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (on_locally_owned_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();
                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];

                    bool fixed = pmm->get_fixed_flag(&node);
                    if (fixed) ++num_fixed;

                    double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                    m_current_position[&node] = Vector(coord_current, coord_current+spatialDim);
                    m_delta[&node] = Vector(spatialDim, 0.0);
                    m_weight[&node] = Vector(spatialDim, 0.0);
                    m_nweight[&node] = Vector(spatialDim, 0.0);
                  }
              }
          }
        std::cout << "num_fixed= " << num_fixed << std::endl;
      }

      // element loop: compute deltas; all elements (even ghosts) contribute to the non-ghosted nodes
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (on_locally_owned_part(**k))  
            // loop over all elements
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                std::cout << "tmp srk num_elements_in_bucket= " << num_elements_in_bucket << std::endl;

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

                        bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                        bool fixed = pmm->get_fixed_flag(&node);
                        if (fixed || isGhostNode)
                          continue;

                        double *coord_current = PerceptMesh::field_data(coord_field_current, node);
                        
                        Vector& delta = m_delta[&node];
                        Vector& weight = m_weight[&node];
                        Vector& nweight = m_nweight[&node];
                        
                        double eps = 1.e-6;
                        double eps1 = eps*edge_length_ave;

                        for (int idim=0; idim < spatialDim; idim++)
                          {
                            coord_current[idim] += eps1;
                            double mp = metric(element);
                            coord_current[idim] -= 2.0*eps1;
                            double mm = metric(element);
                            coord_current[idim] += eps1;
                            double dd = (mp - mm)/(2*eps1);
                            delta[idim] += dd;
                            weight[idim] += edge_length_ave;
                            nweight[idim] += 1.0;
                            
                            // FIXME 
                            if (0)
                            {
                              coord_current[idim] -= dd*eps;
                              double m1=metric(element);
                              coord_current[idim] += dd*eps;
                              if (metric_0 > 1.e-6)
                                {
                                  if (!(m1 < metric_0*(1.0+eps)))
                                    {
                                      std::cout << "bad grad" << " m1-metric_0 = " << (m1-metric_0) << std::endl;
                                    }
                                  VERIFY_OP_ON(m1, <, metric_0*(1.0+eps), "bad gradient");
                                }
                            }
                          }

                      }
                  }
              }
          }
      }

      // node loop: get scaling factor; only non-ghost nodes
      if (m_iter == 0)
        m_scale = 1.e-10; // give it a lower limit

      double norm_gradient=0.0;
      {
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->node_rank() );

        if (m_iter == 0)
          {
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
                            double dt = std::fabs(delta[i] / (weight[i]/nweight[i]));
                            m_scale = std::max(m_scale, dt);
                          }
                      }
                  }
              }
            std::cout << "tmp srk scale= " << m_scale;
            if (m_scale > m_max_edge_length_factor)
              m_scale = m_max_edge_length_factor/m_scale;
            else 
              m_scale = 1.0;
          }

        std::cout << "iter= " << m_iter << " scale2= " << m_scale << std::endl;

        m_dmax = 0.0;
        int num_n=0;
        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    ++num_n;
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        // negate since we want to reduce the metric
                        norm_gradient += delta[i]*delta[i];
                        delta[i] = -delta[i] * m_scale;
                        m_dmax = std::max(std::fabs(delta[i]), m_dmax);
                      }
                  }
              }
          }
      }
      norm_gradient = std::sqrt(norm_gradient);
      if (m_dmax < gradNorm && m_untangled) 
        {
          std::cout << "tmp srk already converged m_dmax= " << m_dmax << " m_scale= " << m_scale << " gradNorm= " << gradNorm << " norm_gradient= " << norm_gradient << std::endl;
          //update_node_positions
          return total_metric(mesh,0.0,1.0);
        }

      /// line search
      double alpha = 1.0;
      {
        double metric_1 = total_metric(mesh, 1.e-6, 1.0);
        double metric_0 = total_metric(mesh, 0.0, 1.0);
        std::cout << "tmp srk dmax before line search= " << m_dmax << " norm_gradient= " << norm_gradient << std::endl;
        std::cout << "tmp srk " << " metric_0= " << metric_0 << " metric(1.e-6) = " << metric_1 << " diff= " << metric_1-metric_0 << std::endl;
        metric_1 = total_metric(mesh, -1.e-6, 1.0);
        std::cout << "tmp srk " << " metric_0= " << metric_0 << " metric(-1.e-6)= " << metric_1 << " diff= " << metric_1-metric_0 << std::endl;
        double metric=0.0;
        //double sigma=0.95;
        double tau = 0.5;
        double beta = 1.e-4;
        double armijo_offset_factor = beta*norm_gradient*norm_gradient*m_scale;
        bool converged = false;
        while (!converged)
          {
            metric = total_metric(mesh, alpha, 1.0);
            //converged = (metric > sigma*metric_0) && (alpha > 1.e-16);
            double mfac = alpha*armijo_offset_factor;
            converged = (metric < metric_0 + mfac);
            std::cout << "tmp srk alpha= " << alpha << " metric_0= " << metric_0 << " metric= " << metric << " diff= " << metric - (metric_0 + mfac) << std::endl;
            if (!converged)
              alpha = alpha * tau;
            if (alpha < 1.e-6)
              break;
          }
        //if (metric > sigma*metric_0)
        if (!converged)
          {
            std::cout << "can't reduce metric= " << metric << " metric_0 + armijo_offset " << metric_0+alpha*armijo_offset_factor <<  std::endl;
            do_print_elem_val = true;
            metric_1 = total_metric(mesh, 1.e-6, 1.0);
            metric_0 = total_metric(mesh, 0.0, 1.0);
            do_print_elem_val = false;
            std::cout << "tmp srk " << " metric_0= " << metric_0 << " metric(1.e-6) = " << metric_1 << " diff= " << metric_1-metric_0 << std::endl;
            metric_1 = total_metric(mesh, -1.e-6, 1.0);
            std::cout << "tmp srk " << " metric_0= " << metric_0 << " metric(-1.e-6)= " << metric_1 << " diff= " << metric_1-metric_0 << std::endl;

            throw std::runtime_error("can't reduce metric");
          }
        else
          {
            double a1 = alpha/2.;
            double a2 = alpha;
            double f0 = metric_0, f1 = total_metric(mesh, a1, 1.0), f2 = total_metric(mesh, a2, 1.0);
            double den = 2.*(a2*(-f0 + f1) + a1*(f0 - f2));
            double num = a2*a2*(f1-f0)+a1*a1*(f0-f2);
            if (std::fabs(den) > 1.e-10)
              {
                double alpha_quadratic = num/den;
                if (alpha_quadratic < 2*alpha)
                  {
                    double fm=total_metric(mesh, alpha_quadratic, 1.0);
                    if (fm < f2)
                      {
                        alpha = alpha_quadratic;
                        std::cout << "tmp srk alpha_quadratic= " << alpha_quadratic << " alpha= " << a2 << std::endl;
                      }
                  } 
              }
          }
      }

      update_node_positions(mesh, alpha);

      //MSQ_ERRRTN(err);
      return total_metric(mesh,0.0,1.0);
      
    }
    
    void PMMParallelReferenceMeshSmoother1::update_node_positions(Mesh* mesh, double alpha)
    {
      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
      int spatialDim = eMesh->get_spatial_dim();

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
                    bool isGhostNode = !(on_locally_owned_part(node) || on_globally_shared_part(node));
                    if (fixed || isGhostNode)
                      {
                        continue;
                      }

                    double *coord_current = PerceptMesh::field_data(m_coord_field_current, node);
                    Vector& delta = m_delta[&node];
                    for (int i=0; i < spatialDim; i++)
                      {
                        double dt = alpha*delta[i];
                        m_dmax = std::max(std::fabs(dt), m_dmax);
                        coord_current[i] += dt;  
                      }
                  }
              }
          }
      }
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
                        //double dt = alpha * multiplicative_edge_scaling * delta[i];
                        double dt = alpha * delta[i];
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
                    double mm = metric(element);
                    if (do_print_elem_val) std::cout << "element= " << element.identifier() << " metric= " << mm << std::endl;
                    mtot += mm;
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
