/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMSmootherMetric_hpp
#define PMMSmootherMetric_hpp

#if !defined(__IBMCPP__) && defined(STK_BUILT_IN_SIERRA)

#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/JacobianUtil.hpp>
#include <stk_percept/mesh/mod/mesquite-interface/PMMMsqMatrix.hpp>
#include <stk_percept/math/Math.hpp>

namespace stk {
  namespace percept {

    enum CombineOp {
      COP_SUM,
      COP_MIN,
      COP_MAX
    };

    class PMMSmootherMetric
    {
    public:
      PMMSmootherMetric(PerceptMesh *eMesh) : m_topology_data(0), m_eMesh(eMesh), m_node(0), m_is_nodal(false), m_combine(COP_SUM)
      {
        m_coord_field_current   = eMesh->get_coordinates_field();
        m_coord_field_original  = eMesh->get_field("coordinates_NM1");
      }

      virtual double metric(stk::mesh::Entity& element, bool& valid)=0;

      const CellTopologyData * m_topology_data ;
      void set_node(stk::mesh::Entity *node) { m_node=node; }
      stk::mesh::Entity *get_node() { return m_node; }
      void set_combine_op(CombineOp combine) { m_combine= combine; }
      CombineOp get_combine_op() { return m_combine; }
      bool is_nodal() { return m_is_nodal; }

    protected:
      PerceptMesh *m_eMesh;
      stk::mesh::Entity *m_node; // for metrics that are node-based
      bool m_is_nodal;
      CombineOp m_combine;
      stk::mesh::FieldBase *m_coord_field_current;
      stk::mesh::FieldBase *m_coord_field_original;
      
    };


    class PMMSmootherMetricUntangle : public PMMSmootherMetric
    {
      double m_beta_mult;
    public:
      PMMSmootherMetricUntangle(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {
        //int spatialDim= eMesh->get_spatial_dim();
        m_beta_mult = 0.05;
      }
      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA, jacW;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_untangle=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            if (detAi <= 0.)
              {
                valid = false;
              }
#if 0
            double untangle_metric = 0.0;
            double beta = m_beta_mult*detWi;
            double temp_var = detAi - beta;

            double fval=0.0;
            if(temp_var<0.0){
              fval = -temp_var;
            }
            else
              {
                //fval = -0.001*temp_var;
              }

#endif
            //fval = Math::my_max_hi(-temp_var,0.0,beta*0.001);
            val_untangle += std::max(-(detAi - m_beta_mult*detWi),0.0);
          }
        val = val_untangle;
        return val;
      }
    };

    class PMMSmootherMetricShapeSizeOrient : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricShapeSizeOrient(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}
      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA, jacW;
        //jacA.m_scale_to_unit = true;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        MsqMatrix<3,3> Ident; 
        //ident.identity();
        identity(Ident);

        MsqMatrix<3,3> AI, Atmp, WAI;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            //double detWi = jacW.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            if (std::fabs(detAi) > 1.e-10)
              {
                //shape_metric = sqr_Frobenius(jacW.m_J[i]*inverse(jacA.m_J[i]) - Ident);
                MsqMatrix<3,3>& W = jacW.m_J[i];
                MsqMatrix<3,3>& A = jacA.m_J[i];
                inverse(A, AI);
                product(W, AI, WAI);
                difference(WAI, Ident, Atmp);
                shape_metric = my_sqr_Frobenius(Atmp);
                //VERIFY_OP_ON(std::fabs(shape_metric_new - shape_metric), <, 1.e-5, "hmm");
              }
            //val_shape += std::sqrt(shape_metric);
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }
    };

    class PMMSmootherMetricScaledJacobian0 : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricScaledJacobian0(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) 
      { m_is_nodal=true; m_combine=COP_MAX; }

      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        VERIFY_OP_ON(m_node, !=, 0, "must set a node");
        valid = true;
        JacobianUtil jacA, jacSA, jacW;
        jacSA.m_scale_to_unit = true;

        double SA_ = 0.0, A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacSA(SA_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        stk::mesh::PairIterRelation elem_nodes = element.relations(m_eMesh->node_rank());
        VERIFY_OP_ON((int)elem_nodes.size(), ==, jacA.m_num_nodes, "node num mismatch");
        val_shape = 0.0;
        bool found = false;
        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            if (elem_nodes[i].entity() == m_node)
              {
                double detAi = jacA.m_detJ[i];
                double detSAi = jacSA.m_detJ[i];
                double detWi = jacW.m_detJ[i];
                if (detAi < 0)
                  {
                    valid = false;
                  }
                double shape_metric = 0.0;
                //MsqMatrix<3,3>& A = jacA.m_J[i];
                double scale_factor = detWi;
                scale_factor = 1.0;
                double fac = 0.2;
                shape_metric = scale_factor* (detSAi > fac ? fac : detSAi);
                val_shape = shape_metric;
                //std::cout << "tmp srk i= " << i << " detAi = " << detAi << " detSAi= " << detSAi << " shape_metric= " << shape_metric << " val_shape= " << val_shape << std::endl;
                found = true;
                break;
              }
          }
        VERIFY_OP_ON(found, ==, true, "logic err");
        val = -val_shape;
        //std::cout << "tmp srk val = " << val << std::endl;
        return val;
      }
    };

    class PMMSmootherMetricScaledJacobian : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricScaledJacobian(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}

      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA, jacSA, jacW;
        jacSA.m_scale_to_unit = true;

        double SA_ = 0.0, A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacSA(SA_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        val_shape = 0.0;
        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detSAi = jacSA.m_detJ[i];
            double detWi = jacW.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            //MsqMatrix<3,3>& A = jacA.m_J[i];
            double scale_factor = detWi;
            scale_factor = 1.0;
            double sign_SA = (detSAi > 0.0 ? 1.0 : -1.0);
            double fac = 0.2;
            shape_metric = scale_factor* (detSAi > fac ? fac*fac : sign_SA*detSAi*detSAi);
            val_shape += shape_metric;
            //std::cout << "tmp srk i= " << i << " detAi = " << detAi << " detSAi= " << detSAi << " shape_metric= " << shape_metric << " val_shape= " << val_shape << std::endl;
          }

        val = -val_shape;
        //std::cout << "tmp srk val = " << val << std::endl;
        return val;
      }
    };

    class PMMSmootherMetricShapeB1 : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricShapeB1(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}

      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA, jacW;
        //jacA.m_scale_to_unit = true;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        MsqMatrix<3,3> Ident; 
        //ident.identity();
        identity(Ident);

        MsqMatrix<3,3> WI, T;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            //double detWi = jacW.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            if (std::fabs(detAi) > 1.e-10)
              {
                //shape_metric = sqr_Frobenius(jacW.m_J[i]*inverse(jacA.m_J[i]) - Ident);
                MsqMatrix<3,3>& W = jacW.m_J[i];
                MsqMatrix<3,3>& A = jacA.m_J[i];
                inverse(W, WI);
                product(A, WI, T);
                double f = Frobenius(T);
                double d = det(T);
                double den = 3 * MSQ_SQRT_THREE * d;
                shape_metric = (f*f*f)/den - 1.0;
              }
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }

    };


  }
}

#endif
#endif
