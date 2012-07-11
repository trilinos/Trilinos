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

    class PMMParallelReferenceMeshSmoother1 : public PMMParallelReferenceMeshSmoother {
     
    public:  
        
      PMMParallelReferenceMeshSmoother1(int inner_iterations = 100,
                                       double cpu_time = 0.0, 
                                       double grad_norm =1.e-8,
                                       int parallel_iterations = 20)
        : PMMParallelReferenceMeshSmoother(inner_iterations, cpu_time, grad_norm, parallel_iterations)
      {}


    protected:

      virtual void run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                                      MsqError& err );

      virtual double total_metric(Mesh *mesh, double alpha, double edge_scaling);
      virtual double metric(stk::mesh::Entity& entity);
      



    };


  }
}

#endif
#endif
#endif
