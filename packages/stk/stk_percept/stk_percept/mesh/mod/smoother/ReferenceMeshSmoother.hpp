/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef ReferenceMeshSmoother_hpp
#define ReferenceMeshSmoother_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__)

#include <stk_percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <stk_percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <boost/unordered_map.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace stk {
  namespace percept {

    /// A smoother based on a reference mesh - tries to make the new mesh the same local size as original

    class ReferenceMeshSmoother : public MeshSmoother {

    public:

      typedef std::vector<double> Vector;
      typedef boost::unordered_map<stk::mesh::Entity , Vector  > NodeMap;


      ReferenceMeshSmoother(PerceptMesh *eMesh,
                            stk::mesh::Selector *boundary_selector=0,
                            MeshGeometry *meshGeometry=0,
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                            int parallel_iterations = 20)
        : MeshSmoother(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations),
          m_scale(0),
          m_dmax(0),
          m_dnew(0), m_dold(0), m_d0(0), m_dmid(0), m_dd(0), m_alpha(0), m_grad_norm(0), m_grad_norm_scaled(0),
          m_total_metric(0),
          m_stage(0),
          m_omega(0),
          m_omega_prev(0),
          m_iter(0),
          m_num_invalid(0), m_global_metric(std::numeric_limits<double>::max()), m_untangled(false), m_num_nodes(0),
          m_use_ref_mesh(true)
      {}


    protected:

      virtual void run_algorithm();
      virtual double run_one_iteration();

      void sync_fields(int iter=0);
      virtual bool check_convergence();

      template<typename T>
      void check_equal(T& val)
      {
        T global_min = val, global_max=val;
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & global_max ) );
        stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , ReduceMax<1>( & global_max ) );
        VERIFY_OP_ON( global_max, ==, val , "bad parallel val");
        VERIFY_OP_ON( global_min, ==, val , "bad parallel val");
        VERIFY_OP_ON( global_max, ==, global_min , "bad parallel val");
      }

    protected:
      NodeMap m_current_position;
      NodeMap m_delta;
      NodeMap m_weight;
      NodeMap m_nweight;

      double m_scale;
      double m_dmax;
      double m_dnew, m_dold, m_d0, m_dmid, m_dd, m_alpha, m_grad_norm, m_grad_norm_scaled;
      double m_total_metric;
      int m_stage;
      double m_omega;
      double m_omega_prev;
      int m_iter;

      int m_num_invalid;
      double m_global_metric;
      bool m_untangled;
      int m_num_nodes;

      stk::mesh::FieldBase *m_coord_field_original;
      stk::mesh::FieldBase *m_coord_field_projected;
      stk::mesh::FieldBase *m_coord_field_current;
      stk::mesh::FieldBase *m_coord_field_lagged;

      SmootherMetric *m_metric;
    public:
      bool m_use_ref_mesh;
    };


  }
}

#endif
#endif
