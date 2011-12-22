#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>


namespace stk {
  namespace percept {
    using namespace Mesquite;

    // this is a modified copy of Mesquite::LaplaceWrapper's run_wrapper method

    void PMMLaplaceSmoother1::run_wrapper( Mesh* mesh,
                                           ParallelMesh* pmesh,
                                           MeshDomain* geom,
                                           Settings* settings,
                                           QualityAssessor* qa,
                                           MsqError& err )
    {
      //std::cout << "tmp srk start PMMLaplaceSmoother1::run_wrapper..." << std::endl;
      if (get_cpu_time_limit() <= 0.0 && get_vertex_movement_limit_factor() <= 0.0 && get_iteration_limit() <= 0) {
        MSQ_SETERR(err)("No termination criterion set.  "
                        "PMMLaplaceSmoother1 will run forever.", 
                        MsqError::INVALID_STATE);
        return;
      }
  
      IdealWeightInverseMeanRatio qa_metric;
      qa->add_quality_assessment( &qa_metric );
  
      LaplacianSmoother& smoother = m_smoother;
      TerminationCriterion outer, inner;
      if (get_cpu_time_limit() > 0.0)
        outer.add_cpu_time( get_cpu_time_limit() );
      if (get_iteration_limit() > 0)
        outer.add_iteration_limit( get_iteration_limit() );
      if (is_culling_enabled() && get_vertex_movement_limit_factor() > 0.0) {
        inner.cull_on_absolute_vertex_movement_edge_length( get_vertex_movement_limit_factor() );
        smoother.set_inner_termination_criterion( &inner );
      }
      else if (get_vertex_movement_limit_factor() > 0.0) {
        outer.add_absolute_vertex_movement_edge_length( get_vertex_movement_limit_factor() );
      }
      smoother.set_outer_termination_criterion( &outer );
  
      InstructionQueue q;
      q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
      q.set_master_quality_improver( &smoother, err ); MSQ_ERRRTN(err);
      q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
      q.run_common( mesh, pmesh, geom, settings, err ); MSQ_ERRRTN(err);
      //std::cout << "tmp srk start PMMLaplaceSmoother1::run_wrapper...done" << std::endl;
    }
  }
}


#endif
#endif
