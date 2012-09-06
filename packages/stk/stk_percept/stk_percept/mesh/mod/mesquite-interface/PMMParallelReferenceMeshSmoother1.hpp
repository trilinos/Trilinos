/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMParallelReferenceMeshSmoother1_hpp
#define PMMParallelReferenceMeshSmoother1_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PMMParallelReferenceMeshSmoother.hpp>

namespace stk {
  namespace percept {

    using namespace Mesquite;

    /// A Jacobian based optimization smoother, e.g. 1/A - 1/W, W/A - I, etc. (A = local current Jacobian, W is for original mesh)
    /// Conjugate-gradient version, element-based metrics
    class PMMParallelReferenceMeshSmoother1 : public PMMParallelReferenceMeshSmoother {
     
    public:  
        
      /// max_edge_length_factor: used for scaling gradients to approximately this value times local edge length 
      PMMParallelReferenceMeshSmoother1(double max_edge_length_factor=1.0,
                                        int inner_iterations = 100,
                                        double cpu_time = 0.0, 
                                        double grad_norm =1.e-8,
                                        int parallel_iterations = 20)
        : PMMParallelReferenceMeshSmoother(inner_iterations, cpu_time, grad_norm, parallel_iterations),
          m_max_edge_length_factor(max_edge_length_factor), m_scale(1.), m_pmm(0)
      {}


    protected:
      double m_max_edge_length_factor;
      double m_scale;

      virtual void get_gradient( Mesh* mesh, MeshDomain *domain);
      virtual void get_scale( Mesh* mesh, MeshDomain *domain);
      void debug_print(double alpha);

      virtual double run_one_iteration( Mesh* mesh,  MeshDomain *domain,
                                      MsqError& err );

      virtual double total_metric(Mesh *mesh, double alpha, double multiplicative_edge_scaling, bool& valid);
      virtual double metric(stk::mesh::Entity& entity, bool& valid);
      virtual void update_node_positions(Mesh* mesh, double alpha);
      virtual bool check_convergence();

      double nodal_metric(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d,  bool& valid );
      void nodal_gradient(stk::mesh::Entity& node, double alpha, double *coord_current, double *cg_d,  bool& valid, double *ng);
      double nodal_edge_length_ave(stk::mesh::Entity& node);
      
      PerceptMesquiteMesh *m_pmm;

    };


  }
}

#endif
#endif

