/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMLaplaceSmoother1_hpp
#define PMMLaplaceSmoother1_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

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
#include <Wrapper.hpp>

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

    class PMMLaplaceSmoother1 : public Mesquite::LaplaceWrapper
    {
      Mesquite::LaplacianSmoother m_smoother;
      int m_numIterMax;
    public:
  
      PMMLaplaceSmoother1(int numIterMax=1) : Mesquite::LaplaceWrapper(), m_numIterMax(numIterMax) {}

      virtual ~PMMLaplaceSmoother1() {}
      Mesquite::LaplacianSmoother& get_smoother() { return m_smoother; }

      void run(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain, bool always_smooth=true, int debug=0)
      {
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
            num_invalid = PMMShapeImprover::count_invalid_elements(mesh, domain);
            std::cout << "tmp srk PMMLaplaceSmoother1 num_invalid before= " << num_invalid 
                      << (num_invalid ? " WARNING: invalid elements exist before Mesquite smoothing" : " ")
                      << std::endl;
          }

        if (num_invalid || always_smooth)
          {

            this->set_iteration_limit(m_numIterMax);
            this->run_instructions(&mesh, &domain, mErr);
            if (check_quality)
              {
                num_invalid = PMMShapeImprover::count_invalid_elements(mesh, domain);
                std::cout << "tmp srk PMMLaplaceSmoother1 num_invalid after= " << num_invalid << " " 
                          << (num_invalid ? " ERROR still have invalid elements after Mesquite smoothing" : 
                              " SUCCESS: smoothed and removed invalid elements ")
                          << std::endl;
              }
          }
      }
      
    protected:

      virtual void run_wrapper( Mesquite::Mesh* mesh,
                                Mesquite::ParallelMesh* pmesh,
                                Mesquite::MeshDomain* geom,
                                Mesquite::Settings* settings,
                                Mesquite::QualityAssessor* qa,
                                Mesquite::MsqError& err );
  
    private:
  
    };


  }
}

#endif
#endif
#endif
