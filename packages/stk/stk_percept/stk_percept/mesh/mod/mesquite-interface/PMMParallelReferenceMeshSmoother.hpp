/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother_hpp
#define PMMParallelReferenceMeshSmoother_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMSmootherMetric.hpp>
#include <boost/unordered_map.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    /// A weighted Laplace smother - tries to make the new mesh the same local size as original
    class PMMParallelReferenceMeshSmoother : public PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper {
     
    public:  

      typedef std::vector<double> Vector;
      typedef boost::unordered_map<stk::mesh::Entity *, Vector  > NodeMap;
      
        
      PMMParallelReferenceMeshSmoother(int inner_iterations = 100,
                                       double cpu_time = 0.0, 
                                       double grad_norm =1.e-8,
                                       int parallel_iterations = 20)
        : PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper(inner_iterations, cpu_time, grad_norm, parallel_iterations),
          m_num_invalid(0), m_global_metric(std::numeric_limits<double>::max()), m_untangled(false)
      {}


    protected:

      void run_wrapper( Mesh* mesh,
                        ParallelMesh* pmesh,
                        MeshDomain* domain,
                        Settings* settings,
                        QualityAssessor* qa,
                        MsqError& err );

      virtual double run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                                      MsqError& err );

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

      int count_invalid_elements(PerceptMesh *eMesh);

    protected:
      NodeMap m_current_position;
      NodeMap m_delta;
      NodeMap m_weight;
      NodeMap m_nweight;
      double m_dmax;
      double m_dnew, m_dold, m_d0, m_dmid, m_dd, m_alpha, m_grad_norm, m_scaled_grad_norm;
      double m_total_metric;
      int m_stage;
      double m_omega;
      double m_omega_prev;
      int m_iter;
      int m_num_invalid;
      double m_global_metric;
      bool m_untangled;

      PerceptMesquiteMesh *m_pmm;
      PerceptMesh *m_eMesh;

      stk::mesh::FieldBase *m_coord_field_original;
      stk::mesh::FieldBase *m_coord_field_projected;
      stk::mesh::FieldBase *m_coord_field_current;
      stk::mesh::FieldBase *m_coord_field_lagged;

      PMMSmootherMetric *m_metric;
    };


  }
}

#endif
#endif
