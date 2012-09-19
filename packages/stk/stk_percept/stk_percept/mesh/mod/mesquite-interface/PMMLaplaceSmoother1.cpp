#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)


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
      std::cout << "tmp srk start PMMLaplaceSmoother1::run_wrapper... get_iteration_limit() = " << get_iteration_limit() << std::endl;
      if (get_cpu_time_limit() <= 0.0 && get_vertex_movement_limit_factor() <= 0.0 && get_iteration_limit() <= 0) {
        MSQ_SETERR(err)("No termination criterion set.  "
                        "PMMLaplaceSmoother1 will run forever.", 
                        MsqError::INVALID_STATE);
        return;
      }
  
      IdealWeightInverseMeanRatio qa_metric;
      qa->add_quality_assessment( &qa_metric );
  
      LaplacianSmoother& smoother = m_smoother;
      TerminationCriterion outer("<type:laplace_outer>"), inner("<type:laplace_inner>");
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
      std::cout << "tmp srk start PMMLaplaceSmoother1::run_wrapper...done" << std::endl;
    }

    void PMMLaplaceSmoother1::run(Mesquite::Mesh &mesh, Mesquite::MeshDomain &domain, bool always_smooth, int debug)
    {
      Mesquite::ParallelMesh *pmesh = dynamic_cast<Mesquite::ParallelMesh *>(&mesh);

      if (debug)
        {
          Mesquite::MsqDebug::enable(1);
          if (debug > 1) Mesquite::MsqDebug::enable(2);
          if (debug > 2) Mesquite::MsqDebug::enable(3);
        }

      Mesquite::MsqError mErr;

      int num_invalid = 0;
      bool check_quality=true;
      if (check_quality)
        {
          num_invalid = PMMShapeImprover::count_invalid_elements(mesh, pmesh, domain);
          std::cout << "tmp srk PMMLaplaceSmoother1 num_invalid before= " << num_invalid 
                    << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : " ")
                    << std::endl;
        }

      if (num_invalid || always_smooth)
        {
          std::cout << "tmp srk PMMLaplaceSmoother1 running laplace smoother..." << std::endl;
          if (pmesh)
            this->run_instructions(pmesh, &domain, mErr);
          else
            this->run_instructions(&mesh, &domain, mErr);
          std::cout << "tmp srk PMMLaplaceSmoother1 running laplace smoother...done" << std::endl;
          if (check_quality)
            {
              num_invalid = PMMShapeImprover::count_invalid_elements(mesh, pmesh, domain);
              std::cout << "tmp srk PMMLaplaceSmoother1 num_invalid after= " << num_invalid << " " 
                        << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
                            " SUCCESS: smoothed and removed invalid elements ")
                        << std::endl;
            }
        }
    }


  }
}


#endif
