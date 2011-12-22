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
    public:
  
      PMMLaplaceSmoother1() : Mesquite::LaplaceWrapper() {}

      virtual ~PMMLaplaceSmoother1() {}
      Mesquite::LaplacianSmoother& get_smoother() { return m_smoother; }

      void run(PerceptMesquiteMesh &mesh, PerceptMesquiteMeshDomain &domain)
      {
        Mesquite::MsqError mErr;
        this->set_iteration_limit(1);
        this->run_instructions(&mesh, &domain, mErr);
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
