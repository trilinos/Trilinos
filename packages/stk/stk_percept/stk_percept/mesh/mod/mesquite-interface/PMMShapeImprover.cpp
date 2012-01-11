#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>

namespace stk {
  namespace percept {
    using namespace Mesquite;

    void PMMShapeImprover::PMMShapeImprovementWrapper::run_wrapper( Mesh* mesh,
                                                                    ParallelMesh* pmesh,
                                                                    MeshDomain* domain,
                                                                    Settings* settings,
                                                                    QualityAssessor* qa,
                                                                    MsqError& err )
    {
      // Define an untangler
      UntangleBetaQualityMetric untangle_metric( untBeta );
      LPtoPTemplate untangle_func( 2, &untangle_metric );
      ConjugateGradient untangle_solver( &untangle_func );
      //SteepestDescent untangle_solver( &untangle_func );
      //TerminationCriterion untangle_inner("<type:inner>"), untangle_outer("<type:outer>");
      TerminationCriterion untangle_inner, untangle_outer;
      untangle_solver.use_global_patch();
      untangle_inner.add_absolute_quality_improvement( 0.0 );
      //untangle_inner.add_absolute_gradient_L2_norm( gradNorm );
      //untangle_inner.add_absolute_successive_improvement( successiveEps );
      //untangle_inner.add_relative_successive_improvement( 1.e-6 );
      untangle_inner.add_iteration_limit( 20 );
      untangle_inner.write_iterations("untangle.gpt", err);
      untangle_inner.add_untangled_mesh();

      untangle_outer.add_absolute_quality_improvement( 0.0 );
      untangle_outer.add_iteration_limit( pmesh ? parallelIterations : 1 );

      std::cout << "tmp srk pmesh= " << pmesh << std::endl;
      untangle_solver.set_inner_termination_criterion( &untangle_inner );
      untangle_solver.set_outer_termination_criterion( &untangle_outer );
      //exit(123);

      // define shape improver
      IdealWeightInverseMeanRatio inverse_mean_ratio;
      inverse_mean_ratio.set_averaging_method( QualityMetric::LINEAR );
      LPtoPTemplate obj_func( 2, &inverse_mean_ratio );
      //FeasibleNewton shape_solver( &obj_func );
      ConjugateGradient shape_solver( &obj_func );
      //TerminationCriterion term_inner("<type:inner>"), term_outer("<type:outer>");
      TerminationCriterion term_inner, term_outer;
      shape_solver.use_global_patch();
      qa->add_quality_assessment( &inverse_mean_ratio );
      term_inner.add_absolute_gradient_L2_norm( gradNorm );
      //!term_inner.add_relative_successive_improvement( successiveEps );
      term_inner.add_iteration_limit( 50 );
      term_inner.write_iterations("shape.gpt", err);

      term_outer.add_iteration_limit( pmesh ? parallelIterations : 1 );
      //term_outer.add_absolute_quality_improvement( 1.e-6 );
      term_outer.add_absolute_gradient_L2_norm( gradNorm );
      //!term_outer.add_relative_successive_improvement( successiveEps );

      shape_solver.set_inner_termination_criterion( &term_inner );
      shape_solver.set_outer_termination_criterion( &term_outer );

      // Apply CPU time limit to untangler
      if (maxTime > 0.0)
        untangle_inner.add_cpu_time( maxTime );
  
      Timer totalTimer;

      // Run untangler
      std::cout << "\ntmp srk PMMShapeImprovementWrapper: running untangler...\n " << std::endl;
      bool use_untangle_wrapper = false;
      if (use_untangle_wrapper)
        {
          UntangleWrapper uw;
          //uw.set_untangle_metric(UntangleWrapper::BETA);
          uw.run_instructions(mesh, domain, err);
        }
      else
        {
          InstructionQueue q1;
          q1.set_master_quality_improver( &untangle_solver, err ); MSQ_ERRRTN(err);
          q1.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
          q1.run_common( mesh, pmesh, domain, settings, err ); 
        }
      std::cout << "\ntmp srk PMMShapeImprovementWrapper: running untangler... done\n " << std::endl;
      std::cout << "\ntmp srk PMMShapeImprovementWrapper: MsqError after untangler: " << err << std::endl;

      bool check_quality_after_untangler = true;
      if (check_quality_after_untangler)
        {
          int num_invalid = count_invalid_elements(*mesh, *domain);
          std::cout << "\ntmp srk PMMShapeImprover num_invalid after untangler= " << num_invalid << " " 
                    << (num_invalid ? " ERROR still have invalid elements after Mesquite untangle" : 
                        " SUCCESS: untangled invalid elements ")
                    << std::endl;
          if (num_invalid) return;
        }
      if (m_do_untangle_only) return;
      MSQ_ERRRTN(err);
  

      // If limited by CPU time, limit next step to remaning time
      if (maxTime > 0.0) {
        double remaining = maxTime - totalTimer.since_birth();
        if (remaining <= 0.0 ){
          MSQ_DBGOUT(2) << "Optimization is terminating without perfoming shape improvement." << std::endl;
          remaining = 0.0;
        }
        term_inner.add_cpu_time( remaining );
      }
  
      // Run shape improver
      InstructionQueue q2;
      std::cout << "\ntmp srk PMMShapeImprovementWrapper: running shape improver... \n" << std::endl;
      q2.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
      q2.set_master_quality_improver( &shape_solver, err ); MSQ_ERRRTN(err);
      q2.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
      q2.run_common( mesh, pmesh, domain, settings, err ); 
      std::cout << "\ntmp srk PMMShapeImprovementWrapper: running shape improver... done \n" << std::endl;
      MSQ_ERRRTN(err);
    }


    int PMMShapeImprover::count_invalid_elements(Mesh &mesh, MeshDomain &domain)
    {
      MsqError err;
      InstructionQueue q;
  
#if 1      
      IdealWeightInverseMeanRatio metric;
      metric.set_averaging_method( QualityMetric::LINEAR );
#else
      // Set up barrier metric to see if mesh contains inverted elements
      TShapeB1 mu_b;
      IdealShapeTarget w_ideal;
      TQualityMetric metric( &w_ideal, &mu_b );
#endif  

      // Check for inverted elements in the mesh
      QualityAssessor inv_check( &metric );
      //inv_check.disable_printing_results();
      q.add_quality_assessor( &inv_check, err );  MSQ_ERRZERO(err);
      Settings settings;
      q.run_common( &mesh, 0, &domain, &settings, err ); MSQ_ERRZERO(err);
      //q.remove_quality_assessor( 0, err ); MSQ_ERRZERO(err);
      const QualityAssessor::Assessor* inv_b = inv_check.get_results( &metric );
      int num_invalid = inv_b->get_invalid_element_count();
      return num_invalid;
    }

    void PMMShapeImprover::run(Mesquite::Mesh &mesh, Mesquite::MeshDomain &domain, bool always_smooth, int debug)
    {
#ifdef USE_CALLGRIND
      CALLGRIND_START_INSTRUMENTATION
        CALLGRIND_TOGGLE_COLLECT
#endif
        if (debug)
          {
            Mesquite::MsqDebug::enable(1);
            if (debug > 1) Mesquite::MsqDebug::enable(2);
            if (debug > 2) Mesquite::MsqDebug::enable(3);
          }

      Mesquite::ParallelMesh *pmesh = dynamic_cast<Mesquite::ParallelMesh *>(&mesh);
      std::cout << "tmp srk PMMShapeImprover::run: pmesh= " << pmesh << std::endl;

      Mesquite::MsqError mErr;
      int num_invalid = 0;
      bool check_quality=true;
      if (check_quality)
        {
          num_invalid = count_invalid_elements(mesh, domain);
          std::cout << "\ntmp srk PMMShapeImprover num_invalid before= " << num_invalid 
                    << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : 
                        (!always_smooth ? "WARNING: no smoothing requested since always_smooth=false" : " "))
                    << std::endl;
        }

      if (num_invalid || always_smooth)
        {
          bool use_canned_wrapper = false;
          if (use_canned_wrapper)
            {
              Mesquite::ShapeImprovementWrapper siw(mErr);
              if (pmesh)
                siw.run_instructions(pmesh, &domain, mErr);
              else
                siw.run_instructions(&mesh, &domain, mErr);
            }
          else
            {
              int  msq_debug             = 2; // 1,2,3 for more debug info
              bool always_smooth_local   = false;
              bool do_laplace            = false;
              bool do_jacobi             = true;

              // Define a Laplace smoother
              if (do_laplace)
                {
                  int num_laplace_iter = 1;
                  PMMLaplaceSmoother1 ls(num_laplace_iter);
                  if (do_jacobi) ls.get_smoother().do_jacobi_optimization();
                  ls.run(mesh, domain, always_smooth_local, msq_debug);
                }

              bool do_untangle_only = false;
              PMMShapeImprovementWrapper siw(mErr);
              siw.m_do_untangle_only = do_untangle_only;
              if (pmesh)
                siw.run_instructions(pmesh, &domain, mErr);
              else
                siw.run_instructions(&mesh, &domain, mErr);
            }

          std::cout << "\ntmp srk PMMShapeImprover: MsqError after ShapeImprovementWrapper: " << mErr << std::endl;

          if (check_quality)
            {
              num_invalid = count_invalid_elements(mesh, domain);
              std::cout << "\ntmp srk PMMShapeImprover num_invalid after= " << num_invalid << " " 
                        << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
                            " SUCCESS: smoothed and removed invalid elements ")
                        << std::endl;
            }

          MSQ_ERRRTN(mErr);

        }
#ifdef USE_CALLGRIND
      CALLGRIND_TOGGLE_COLLECT
        CALLGRIND_STOP_INSTRUMENTATION
#endif
        }

  }
}


#endif
#endif
