/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef ReferenceMeshSmoother1_hpp
#define ReferenceMeshSmoother1_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/smoother/ReferenceMeshSmoother.hpp>

namespace stk {
  namespace percept {

    /// A Jacobian based optimization smoother, e.g. 1/A - 1/W, W/A - I, etc. (A = local current Jacobian, W is for original mesh)
    /// Conjugate-gradient version, element-based metrics
    class ReferenceMeshSmoother1 : public ReferenceMeshSmoother {
     
    public:  
        
      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length 
      ReferenceMeshSmoother1(PerceptMesh *eMesh,
                            stk::mesh::Selector *boundary_selector=0,
                            MeshGeometry *meshGeometry=0, 
                            int inner_iterations = 100,
                            double grad_norm =1.e-8,
                            int parallel_iterations = 20)
        : ReferenceMeshSmoother(eMesh, boundary_selector, meshGeometry, inner_iterations, grad_norm, parallel_iterations)
        ,m_max_edge_length_factor(1.0)
      {}

    protected:
      double m_max_edge_length_factor;

      virtual void get_gradient();
      virtual void get_scale();
      void debug_print(double alpha);

      virtual double run_one_iteration();


      virtual double total_metric( double alpha, double multiplicative_edge_scaling, bool& valid, int *num_invalid=0);
      virtual double metric(stk::mesh::Entity& entity, bool& valid);
      virtual void update_node_positions( double alpha);
      virtual bool check_convergence();

      double nodal_metric(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d,  bool& valid );
      void nodal_gradient(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d,  bool& valid, double *ng);
      double nodal_edge_length_ave(stk::mesh::Entity& node);
      void get_edge_lengths(PerceptMesh * eMesh);
      
    };


  }
}

#endif
#endif

