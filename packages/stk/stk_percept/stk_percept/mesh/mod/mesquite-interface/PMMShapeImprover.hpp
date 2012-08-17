/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMShapeImprover_hpp
#define PMMShapeImprover_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <InstructionQueue.hpp>
#include <Settings.hpp>
#include <ShapeImprovementWrapper.hpp>
#include <UntangleWrapper.hpp>

#include <IdealWeightInverseMeanRatio.hpp> 
#include <QualityAssessor.hpp>
#include <TerminationCriterion.hpp>

#include <TQualityMetric.hpp>
#include <TShapeB1.hpp>
#include <IdealShapeTarget.hpp>

#include <MsqTimer.hpp>
#include <UntangleBetaQualityMetric.hpp>
#include <LPtoPTemplate.hpp>
#include <ConjugateGradient.hpp>
#include <SteepestDescent.hpp>
#include <FeasibleNewton.hpp>

#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>

#undef USE_CALLGRIND
//#define USE_CALLGRIND
#ifdef USE_CALLGRIND
#include "/usr/netpub/valgrind-3.6.0/include/valgrind/callgrind.h"
#endif

namespace stk {
  namespace percept {

    using namespace Mesquite;

    const double DEF_UNT_BETA = 1e-8;
    const double DEF_SUC_EPS = 1e-4;

    class PMMShapeImprover
    {
    public:
      class PMMShapeImprovementWrapper : public Wrapper {
     
      public:  
        
        //Constructor sets the instructions in the queue.
        PMMShapeImprovementWrapper(int inner_iterations = 100,
                                   double cpu_time = 0.0, 
                                   double grad_norm =1.e-8,
                                   int parallel_iterations = 20)
          : innerIter(inner_iterations),
            maxTime(cpu_time), 
            gradNorm(grad_norm),
            untBeta(DEF_UNT_BETA),
            successiveEps(DEF_SUC_EPS),
            parallelIterations(parallel_iterations),
            m_do_untangle_only(false)
        {}


      protected:

        void run_wrapper( Mesh* mesh,
                          ParallelMesh* pmesh,
                          MeshDomain* domain,
                          Settings* settings,
                          QualityAssessor* qa,
                          MsqError& err );
      
      private:

        int innerIter;
        double maxTime, gradNorm;
        // constants
        const double untBeta;
        const double successiveEps;
        int parallelIterations;
      public:
        bool m_do_untangle_only;


      };

    private:
      int innerIter;
      double gradNorm;
      int parallelIterations;
    public:

      PMMShapeImprover(int innerIter=100, double gradNorm = 1.e-8, int parallelIterations=20) : 
        innerIter(innerIter), gradNorm(gradNorm), parallelIterations(parallelIterations)
      {}

      static int count_invalid_elements(Mesh &mesh, ParallelMesh* pmesh, MeshDomain &domain);

      void run(Mesquite::Mesh &mesh, Mesquite::MeshDomain &domain, bool always_smooth=true, int debug=0);

      static void save_or_restore_debug_state(bool save)
      {
        static bool debug[3] = {false,false,false};
        if (save) 
          {
            debug[0] = MsqDebug::get(1);
            debug[1] = MsqDebug::get(2);
            debug[2] = MsqDebug::get(3);
          }
        else
          {
            if (debug[0]) MsqDebug::enable(1);
            if (debug[1]) MsqDebug::enable(2);
            if (debug[2]) MsqDebug::enable(3);
          }
      }

    };

  }
}

#endif
#endif
