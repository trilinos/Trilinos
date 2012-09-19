#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeSizeOrientImprover.hpp>

#include <Mesquite.hpp>
#include <ViscousCFDTetShapeWrapper.hpp>

#include <QualityAssessor.hpp>
#include <InstructionQueue.hpp>
#include <TagVertexMesh.hpp>
#include <TrustRegion.hpp>
#include <FeasibleNewton.hpp>
#include <TerminationCriterion.hpp>

#include <PMeanPTemplate.hpp>
#include <TQualityMetric.hpp>
#include <AddQualityMetric.hpp>

#include <TShapeB1.hpp>
#include <TShapeNB1.hpp>
#include <TShapeSizeOrientB1.hpp>
#include <TShapeSizeOrientNB1.hpp>

#include <IdealShapeTarget.hpp>
#include <RefMeshTargetCalculator.hpp>
#include <ReferenceMesh.hpp>
#include <TetDihedralWeight.hpp>
#include <RemainingWeight.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    void PMMShapeSizeOrientImprover::run_wrapper( Mesh* mesh,
                                                  ParallelMesh* pmesh,
                                                  MeshDomain* domain,
                                                  Settings* settings,
                                                  QualityAssessor* qa,
                                                  MsqError& err )
    {
      InstructionQueue q;
  
      // Set up barrier metric to see if mesh contains inverted elements
      TShapeB1 mu_b;
      IdealShapeTarget w_ideal;
      TQualityMetric barrier( &w_ideal, &mu_b );
  
      // Check for inverted elements in the mesh
      QualityAssessor inv_check( &barrier );
      //inv_check.disable_printing_results();
      q.add_quality_assessor( &inv_check, err );  MSQ_ERRRTN(err);
      q.run_common( mesh, pmesh, domain, settings, err ); MSQ_ERRRTN(err);
      q.remove_quality_assessor( 0, err ); MSQ_ERRRTN(err);
      const QualityAssessor::Assessor* inv_b = inv_check.get_results( &barrier );
      //const bool use_barrier = (0 == inv_b->get_invalid_element_count());
      std::cout << "tmp srk PMMShapeSizeOrientImprover::run_wrapper get_invalid_element_count= " 
                << inv_b->get_invalid_element_count() << std::endl;
  
      // Create remaining metric instances
      TShapeNB1 mu;
      TShapeSizeOrientNB1 mu_o;
      TShapeSizeOrientB1 mu_ob;
  
      // Select which target metrics to use
      //TMetric *mu_p, *mu_op;
      //if (use_barrier) {
      //  mu_p = &mu_b;
      //  mu_op = &mu_ob;
      //}
      //else {
      //  mu_p = &mu;
      //  mu_op = &mu_o;
      //}
  
      // Set up target and weight calculators
      std::string altCoordName = "msq_jacobi_temp_coords"; // FIXME
      bool should_clean_up_tag_data = true; // default
      TagVertexMesh init_mesh( err, pmesh ? (Mesh*)pmesh : mesh, should_clean_up_tag_data, altCoordName );  MSQ_ERRRTN(err);
      ReferenceMesh ref_mesh( &init_mesh );
      RefMeshTargetCalculator w_init( &ref_mesh );

      //TetDihedralWeight c_dihedral( &ref_mesh, dCutoff, aVal );
      //RemainingWeight c_remaining( &c_dihedral );
  
      // Create objective function
      //       TQualityMetric metric1( &w_ideal, &c_dihedral,  mu_p  );
      //       TQualityMetric metric2( &w_init,  &c_remaining, mu_op );
      //       AddQualityMetric of_metric( &metric1, &metric2, err );  MSQ_ERRRTN(err);
      TQualityMetric of_metric( &w_init, &mu_o );
      PMeanPTemplate obj_func( 1.0, &of_metric );
  
      // Create optimizer
      //TrustRegion solver( &obj_func );
      FeasibleNewton solver( &obj_func );
      TerminationCriterion term, ptc;
      term.add_iteration_limit( iterationLimit );
      term.add_absolute_vertex_movement( maxVtxMovement );
      ptc.add_iteration_limit( pmesh ? parallelIterations : 1 );
      solver.set_inner_termination_criterion( &term );
      solver.set_outer_termination_criterion( &ptc );
  
      // Create instruction queue
      //qa->add_quality_assessment( &metric1 );
      //qa->add_quality_assessment( &metric2 );
      qa->add_quality_assessment( &of_metric );
      q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);
      q.set_master_quality_improver( &solver, err ); MSQ_ERRRTN(err);
      q.add_quality_assessor( qa, err ); MSQ_ERRRTN(err);

      // Optimize mesh
      q.run_common( mesh, pmesh, domain, settings, err ); MSQ_CHKERR(err);  
    }

  }
}
#endif
