#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother2.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother3.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMSmootherMetric.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMLaplaceSmoother1.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>

#include "PMeanPTemplate.hpp"
#include "TQualityMetric.hpp"
#include "AddQualityMetric.hpp"

#include "TShapeB1.hpp"
#include "TShapeNB1.hpp"

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {
    using namespace Mesquite;

    // this is a sample (unfinished) implementation of how we might use Mesquite (needs to have global
    // updates of ConjugateGradient global quantities, etc, to make it produce parallel/serial consistency)
    double PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper::run_one_iteration( Mesh* mesh, MeshDomain *domain,
                                                                                          MsqError& err )
    {
      std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper::run_one_iteration start..." << std::endl;

      // define shape improver
      IdealWeightInverseMeanRatio inverse_mean_ratio;
      inverse_mean_ratio.set_averaging_method( QualityMetric::LINEAR );
      LPtoPTemplate obj_func( 2, &inverse_mean_ratio );

      ConjugateGradient shape_solver( &obj_func );
      TerminationCriterion term_inner("<type:shape_inner>"), term_outer("<type:shape_outer>");
      term_inner.write_iterations("shape.gpt", err);

      shape_solver.use_global_patch();
      //!!! qa->add_quality_assessment( &inverse_mean_ratio );

      //!term_inner.add_relative_successive_improvement( successiveEps );

      //term_inner.add_absolute_gradient_L2_norm( gradNorm );
      //term_inner.add_absolute_vertex_movement(0.0);
      term_inner.add_iteration_limit( 1 );

      //term_outer.add_absolute_gradient_L2_norm( gradNorm );
      //term_outer.add_absolute_vertex_movement(0.0);
      term_outer.add_iteration_limit( 1 );

      //term_outer.add_absolute_quality_improvement( 1.e-6 );
      //!term_outer.add_relative_successive_improvement( successiveEps );

      shape_solver.set_inner_termination_criterion( &term_inner );
      shape_solver.set_outer_termination_criterion( &term_outer );

      Timer totalTimer;

      // Run shape improver
      InstructionQueue q2;
      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper: running shape improver... \n" << std::endl;

      QualityAssessor qa_check( &inverse_mean_ratio );

      Settings settings;

      q2.add_quality_assessor( &qa_check, err ); 
      q2.set_master_quality_improver( &shape_solver, err );
      q2.run_common( mesh, 0, domain, &settings, err ); 

      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper: running shape improver... done \n" << std::endl;

      return 0;
    }

    void PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper::run_wrapper( Mesh* mesh,
                                                                    ParallelMesh* pmesh,
                                                                    MeshDomain* domain,
                                                                    Settings* settings,
                                                                    QualityAssessor* qa,
                                                                    MsqError& err )
    {
#if 0
      std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper innerIter= " << innerIter << " parallelIterations= " << parallelIterations << std::endl;

      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper: running shape improver... \n" << std::endl;

      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();
      stk::mesh::FieldBase *coord_field = eMesh->get_coordinates_field();
      stk::mesh::FieldBase *coord_field_current   = coord_field;
      stk::mesh::FieldBase *coord_field_projected = eMesh->get_field("coordinates_N"); 
      stk::mesh::FieldBase *coord_field_original  = eMesh->get_field("coordinates_NM1");

      //double alphas[] = {0.0,0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0};
      //double alphas[] = {0.001,0.01,0.1,0.2,0.4,0.6,0.8,1.0};
      double alphas[] = {0.001};
      int nalpha = sizeof(alphas)/sizeof(alphas[0]);
      
      for (int outer = 0; outer < nalpha; outer++)
        {
          double alpha = alphas[outer];

          // set current state and evaluate mesh validity
          eMesh->nodal_field_axpbypgz(alpha, coord_field_projected, (1.0-alpha), coord_field_original, 0.0, coord_field_current);

          int num_invalid = parallel_count_invalid_elements(eMesh);
          if (!get_parallel_rank()) 
            std::cout << "\ntmp srk PMMParallelShapeImprover num_invalid current= " << num_invalid 
                      << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : "OK")
                      << std::endl;

#if 0
          for (int iter = 0; iter < innerIter; iter++)
            {
              //
              int num_invalid = parallel_count_invalid_elements(eMesh);
              if (!get_parallel_rank()) 
                std::cout << "\ntmp srk PMMParallelShapeImprover num_invalid current= " << num_invalid 
                          << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : "OK")
                          << std::endl;
              run_one_iteration(mesh, err);
              sync_fields();
              check_convergence();
            }
#endif
        }

      //if (!get_parallel_rank()) 
      std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelShapeImprovementWrapper: running shape improver... done \n" << std::endl;

      MSQ_ERRRTN(err);
#endif
    }


    /// preferred for parallel
    int PMMParallelShapeImprover::parallel_count_invalid_elements(PerceptMesh *eMesh)
    {
      PMMSmootherMetricUntangle utm(eMesh);

      int num_invalid=0;
      // element loop
      {
        stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
        const std::vector<stk::mesh::Bucket*> & buckets = eMesh->get_bulk_data()->buckets( eMesh->element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (PerceptMesquiteMesh::select_bucket(**k, eMesh) && on_locally_owned_part(**k))  
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];
                    bool valid=true;
                    utm.metric(element, valid);
                    if (!valid)
                      ++num_invalid;
                  }
              }
          }
      }
      stk::all_reduce( MPI_COMM_WORLD, stk::ReduceSum<1>( &num_invalid ) );
      return num_invalid;
    }


    int PMMParallelShapeImprover::count_invalid_elements(Mesh &mesh, MeshDomain *domain)
    {
      MsqError err;
      InstructionQueue q;
      int num_invalid = 0;
      VERIFY_OP_ON(get_parallel_size(), ==, 1, "not ready for parallel; use PerceptMesh form of count_invalid_elements");

      if (1)
        {
          IdealWeightInverseMeanRatio metric;
          //metric.set_averaging_method( QualityMetric::LINEAR );

          // Check for inverted elements in the mesh
          QualityAssessor inv_check( &metric );
          inv_check.disable_printing_results();
          q.add_quality_assessor( &inv_check, err );  MSQ_ERRZERO(err);
          Settings settings;
          q.run_common( &mesh, 0, domain, &settings, err ); MSQ_ERRZERO(err);
          const QualityAssessor::Assessor* inv_b = inv_check.get_results( &metric );
          num_invalid = inv_b->get_invalid_element_count();
        }
      else
        {
          // Set up barrier metric to see if mesh contains inverted elements
          TShapeB1 mu_b;
          IdealShapeTarget w_ideal;
          TQualityMetric barrier( &w_ideal, &mu_b );
  
          // Check for inverted elements in the mesh
          QualityAssessor inv_check( &barrier );
          inv_check.disable_printing_results();
          q.add_quality_assessor( &inv_check, err ); MSQ_ERRZERO(err);
          Settings settings;
          q.run_common( &mesh, 0, domain, &settings, err ); MSQ_ERRZERO(err);
          const QualityAssessor::Assessor* inv_b = inv_check.get_results( &barrier );
          num_invalid = inv_b->get_invalid_element_count();
        }

      stk::all_reduce( MPI_COMM_WORLD, stk::ReduceSum<1>( &num_invalid ) );
      
      return num_invalid;
    }

    void PMMParallelShapeImprover::run(Mesquite::Mesh &mesh, Mesquite::MeshDomain *domain, bool always_smooth, int debug)
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
      if (!get_parallel_rank()) std::cout << "tmp srk PMMParallelShapeImprover::run: pmesh= " << pmesh << std::endl;

      PerceptMesquiteMesh *pmm = dynamic_cast<PerceptMesquiteMesh *>(&mesh);
      PerceptMesh *eMesh = pmm->getPerceptMesh();

      Mesquite::MsqError mErr;
      int num_invalid = parallel_count_invalid_elements(eMesh);
      if (!get_parallel_rank()) 
        std::cout << "\ntmp srk PMMParallelShapeImprover num_invalid before= " << num_invalid 
                      << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : 
                          (!always_smooth ? "WARNING: no smoothing requested since always_smooth=false" : " "))
                      << std::endl;
      //if (num_invalid) throw std::runtime_error("PMMParallelShapeImprover can't start from invalid mesh...");

      if (always_smooth)
        {
          //int  msq_debug             = debug; // 1,2,3 for more debug info

          bool do_untangle_only = false;
          std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk innerIter= " << innerIter << " parallelIterations= " << parallelIterations << std::endl;
          //PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper siw(innerIter, 0.0, gradNorm, parallelIterations);
          //PMMParallelReferenceMeshSmoother siw(innerIter, 0.0, gradNorm, parallelIterations);
          //PMMParallelReferenceMeshSmoother1 siw(0.05, innerIter, 0.0, gradNorm, parallelIterations);
          PMMParallelReferenceMeshSmoother1 siw(0.05, innerIter, 0.0, gradNorm, parallelIterations);
          siw.m_do_untangle_only = do_untangle_only;
          siw.run_instructions(&mesh, domain, mErr);

          //if (!get_parallel_rank()) 
          std::cout << "\nP[" << get_parallel_rank() << "] tmp srk PMMParallelShapeImprover: MsqError after ShapeImprovementWrapper: " << mErr << std::endl;

          num_invalid = parallel_count_invalid_elements(eMesh);
          //if (!get_parallel_rank()) 
          std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk PMMParallelShapeImprover num_invalid after= " << num_invalid << " " 
                    << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
                        " SUCCESS: smoothed and removed invalid elements ")
                    << std::endl;
          MPI_Barrier( MPI_COMM_WORLD );
          std::cout << "\nP[" << Mesquite::get_parallel_rank() << "] tmp srk after barrier" << std::endl;
        }

          MSQ_ERRRTN(mErr);

#ifdef USE_CALLGRIND
      CALLGRIND_TOGGLE_COLLECT
        CALLGRIND_STOP_INSTRUMENTATION
#endif
    }


  }
}


#endif
