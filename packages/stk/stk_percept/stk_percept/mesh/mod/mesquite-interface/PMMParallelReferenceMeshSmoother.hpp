/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother_hpp
#define PMMParallelReferenceMeshSmoother_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelShapeImprover.hpp>
#include <boost/unordered_map.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    class PMMParallelReferenceMeshSmoother : public PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper {
     
    public:  

      typedef std::vector<double> Vector;
      typedef boost::unordered_map<stk::mesh::Entity *, Vector  > NodeMap;
      
        
      PMMParallelReferenceMeshSmoother(int inner_iterations = 100,
                                       double cpu_time = 0.0, 
                                       double grad_norm =1.e-8,
                                       int parallel_iterations = 20)
        : PMMParallelShapeImprover::PMMParallelShapeImprovementWrapper(inner_iterations, cpu_time, grad_norm, parallel_iterations)
      {}


    protected:

      void run_wrapper( Mesh* mesh,
                        ParallelMesh* pmesh,
                        MeshDomain* domain,
                        Settings* settings,
                        QualityAssessor* qa,
                        MsqError& err );

      void run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                              MsqError& err );

      void sync_fields(int iter=0);
      bool check_convergence();
      
    private:
      NodeMap m_current_position;
      NodeMap m_delta;
      NodeMap m_weight;
      double m_dmax;
      double m_alpha;
      double m_alpha_prev;

      PerceptMesquiteMesh *m_pmm;
      PerceptMesh *m_eMesh;

    public:


    };


  }
}

#endif
#endif
#endif
