/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother2_hpp
#define PMMParallelReferenceMeshSmoother2_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother1.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    /// A Jacobian based optimization smoother, eg. 1/A - 1/W (A = local current Jacobian, W is for original mesh)
    /// LOCAL patch version
    class PMMParallelReferenceMeshSmoother2 : public PMMParallelReferenceMeshSmoother1 {
     
    public:  
        
      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length 
      PMMParallelReferenceMeshSmoother2(double max_edge_length_factor=1.0,
                                        int inner_iterations = 100,
                                        double cpu_time = 0.0, 
                                        double grad_norm =1.e-8,
                                        int parallel_iterations = 20)
        : PMMParallelReferenceMeshSmoother1(max_edge_length_factor, inner_iterations, cpu_time, grad_norm, parallel_iterations)
      {}


    protected:

      virtual double run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                                      MsqError& err );

    };


  }
}

#endif
#endif

