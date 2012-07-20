/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother_hpp
#define PMMParallelReferenceMeshSmoother_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <boost/unordered_map.hpp>

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
      bool check_convergence();
      
    protected:
      NodeMap m_current_position;
      NodeMap m_delta;
      NodeMap m_weight;
      NodeMap m_nweight;
      double m_dmax;
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

    };


  }
}

#endif
#endif
#endif
