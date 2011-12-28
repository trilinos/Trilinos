/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMShapeImprover_hpp
#define PMMShapeImprover_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <MsqDebug.hpp>
#include <InstructionQueue.hpp>
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

    class PMMShapeImprover
    {
    public:
      PMMShapeImprover() {}

      void run(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain, bool debug=false)
      {
        if (debug)
          {
            Mesquite::MsqDebug::enable(1);
            Mesquite::MsqDebug::enable(2);
            Mesquite::MsqDebug::enable(3);
          }
        Mesquite::MsqError mErr;
        Mesquite::ShapeImprovementWrapper siw(mErr);
        //siw.set_iteration_limit(1);

        siw.run_instructions(&mesh, &domain, mErr);
      }
    };

  }
}

#endif
#endif
#endif
