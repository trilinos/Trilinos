/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMLaplaceSmoother_hpp
#define PMMLaplaceSmoother_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <Mesquite.hpp>
#include <MsqError.hpp>
#include <InstructionQueue.hpp>
#include <SmartLaplacianSmoother.hpp>
#include <LaplaceWrapper.hpp>

#include <IdealWeightInverseMeanRatio.hpp> 
#include <LaplacianSmoother.hpp>
#include <QualityAssessor.hpp>
#include <InstructionQueue.hpp>
#include <TerminationCriterion.hpp>
#include <MsqError.hpp>

#include <stk_percept/mesh/mod/mesquite-interface/PMMShapeImprover.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMeshDomain.hpp>

/// work derived from:
//-------------------------------------------------------------------------
// Filename      : SCVFracMesquiteLaplaceSmoother.hpp
//
// Purpose       : sculptor interface to mesquite Laplacian smoother  
//
// Description   : implements the Laplacian smoother
//
// Creator       : Steve Owen
//
// Creation Date : April 2011
//
// Owner         : Steve Owen
//-------------------------------------------------------------------------

namespace stk {
  namespace percept {

    class PMMLaplaceSmoother
    {
    public:
      PMMLaplaceSmoother() {}

      void run(Mesquite::Mesh &mesh, Mesquite::MeshDomain &domain, bool always_smooth=true, int debug=0)
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
            std::cout << "tmp srk PMMLaplaceSmoother num_invalid before= " << num_invalid 
                      << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : " ")
                      << std::endl;
          }

        if (num_invalid || always_smooth)
          {
            Mesquite::LaplaceWrapper lw;
            lw.set_iteration_limit(1);
            if (pmesh)
              lw.run_instructions(pmesh, &domain, mErr);
            else
              lw.run_instructions(&mesh, &domain, mErr);
            if (check_quality)
              {
                num_invalid = PMMShapeImprover::count_invalid_elements(mesh, pmesh, domain);
                std::cout << "tmp srk PMMLaplaceSmoother num_invalid after= " << num_invalid << " " 
                          << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
                              " SUCCESS: smoothed and removed invalid elements ")
                          << std::endl;
              }
          }
      }
    };

  }
}

#endif
#endif
