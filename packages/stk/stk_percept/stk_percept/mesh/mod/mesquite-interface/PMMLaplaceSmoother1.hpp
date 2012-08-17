/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMLaplaceSmoother1_hpp
#define PMMLaplaceSmoother1_hpp

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
    public:
  
      PMMLaplaceSmoother1(double max_cpu_time=0.0, double max_vertex_movement=1.e-4, int numIterMax=100, bool doCulling=false) : Mesquite::LaplaceWrapper() {
        this->set_cpu_time_limit(max_cpu_time);
        this->set_vertex_movement_limit_factor(max_vertex_movement);
        this->set_iteration_limit(numIterMax);
        this->enable_culling(doCulling);
      }

      virtual ~PMMLaplaceSmoother1() {}
      Mesquite::LaplacianSmoother& get_smoother() { return m_smoother; }

      void run(Mesquite::Mesh &mesh, Mesquite::MeshDomain &domain, bool always_smooth=true, int debug=0);
      
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
