/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMShapeSizeOrientImprover_hpp
#define PMMShapeSizeOrientImprover_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)


#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <InstructionQueue.hpp>
#include <Wrapper.hpp>
#include <ShapeImprovementWrapper.hpp>

#include <IdealWeightInverseMeanRatio.hpp> 
#include <QualityAssessor.hpp>
#include <InstructionQueue.hpp>
#include <TerminationCriterion.hpp>
#include <MsqError.hpp>

#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>


namespace stk {
  namespace percept {

    using namespace Mesquite;

    class PMMShapeSizeOrientImprover : public Wrapper
    {
    private:
      double dCutoff, aVal;
      int iterationLimit;
      int parallelIterations;
      double maxVtxMovement;

      void run_wrapper( Mesh* mesh,
                        ParallelMesh* pmesh,
                        MeshDomain* geom,
                        Settings* settings,
                        QualityAssessor* qa,
                        MsqError& err );

    public:
  
      /**
       *\param max_vertex_movement  Termination optimization if no vertex is moved
       *                            by more than this distance in the previous solver
       *                            step.
       *\param a                    Coefficient for target metric weight
       *\param d_prime              Dihedral handle cut-off for target metric weight
       *\param max_iterations       Termination optimizaiton after this many solver 
       *                            steps.
       */
      PMMShapeSizeOrientImprover( double max_vertex_movement,
                                  double a = 0.4395, 
                                  double d_prime = 135,
                                  int max_iterations = 50,
                                  int parallel_iterations = 10 )
        : dCutoff(d_prime), 
          aVal(a), 
          iterationLimit( max_iterations ),
          parallelIterations( parallel_iterations ),
          maxVtxMovement( max_vertex_movement )
      {}

      void run(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain, int debug=0)
      {
        if (debug)
          {
            Mesquite::MsqDebug::enable(1);
            if (debug > 1) Mesquite::MsqDebug::enable(2);
            if (debug > 2) Mesquite::MsqDebug::enable(3);
          }
        Mesquite::MsqError mErr;

        //this->set_iteration_limit(1);
        Mesquite::ParallelMesh *pmesh = dynamic_cast<Mesquite::ParallelMesh *>(&mesh);
        if (pmesh)
          this->run_instructions(pmesh, &domain, mErr);
        else
          this->run_instructions(&mesh, &domain, mErr);

        //siw.run_instructions(&mesh, &domain, mErr);
      }
    };

  }
}

#endif
#endif
