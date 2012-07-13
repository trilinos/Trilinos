/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother1_hpp
#define PMMParallelReferenceMeshSmoother1_hpp

#if !defined(__IBMCPP__)
#ifdef STK_BUILT_IN_SIERRA

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    /// A Jacobian based optimization smoother - 1/A - 1/W (A = local current Jacobian, W is for original mesh)
    class PMMParallelReferenceMeshSmoother1 : public PMMParallelReferenceMeshSmoother {
     
    public:  
        
      PMMParallelReferenceMeshSmoother1(double max_edge_length_factor=0.05,
                                        int inner_iterations = 100,
                                        double cpu_time = 0.0, 
                                        double grad_norm =1.e-8,
                                        int parallel_iterations = 20)
        : PMMParallelReferenceMeshSmoother(inner_iterations, cpu_time, grad_norm, parallel_iterations),
          m_max_edge_length_factor(max_edge_length_factor)
      {}


    protected:
      double m_max_edge_length_factor;

      virtual void run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                                      MsqError& err );

      virtual double total_metric(Mesh *mesh, double alpha, double multiplicative_edge_scaling=1.0);
      virtual double metric(stk::mesh::Entity& entity);
      

    };


  }
}

#endif
#endif
#endif
