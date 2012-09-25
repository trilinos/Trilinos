/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMSmootherMetric_hpp
#define PMMSmootherMetric_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

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
      virtual double length_scaling_power() { return 1.0; }
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
      virtual double length_scaling_power() { return 3.0; }

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
      virtual double length_scaling_power() { return 1.0; }
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
        identity(Ident);

        MsqMatrix<3,3> AI, Atmp, WAI;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            if (std::fabs(detAi) > 1.e-10)
              {
                MsqMatrix<3,3>& W = jacW.m_J[i];
                MsqMatrix<3,3>& A = jacA.m_J[i];
                inverse(A, AI);
                product(W, AI, WAI);
                difference(WAI, Ident, Atmp);
                shape_metric = my_sqr_Frobenius(Atmp);
              }
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }
    };

    class PMMSmootherMetricScaledJacobianNodal : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricScaledJacobianNodal(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) 
      { m_is_nodal=true; m_combine=COP_MAX; }
      virtual double length_scaling_power() { return 1.0; }

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

    class PMMSmootherMetricScaledJacobianElemental : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricScaledJacobianElemental(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}
      virtual double length_scaling_power() { return 1.0; }

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
        //val_shape = std::numeric_limits<double>::max();
        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            double detSAi = jacSA.m_detJ[i];
            //double detWi = jacW.m_detJ[i];
            double FAi = Frobenius(jacA.m_detJ[i]);
            double FWi = Frobenius(jacW.m_detJ[i]);
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            //MsqMatrix<3,3>& A = jacA.m_J[i];
            double scale_factor = (FAi < FWi ? FAi/FWi : (FAi > 1.e-6? FWi / FAi : 1.0));
            scale_factor = FAi/FWi;
            //scale_factor = 1.0;
            //double sign_SA = (detSAi > 0.0 ? 1.0 : -1.0);
            //double fac = 0.2;
            //shape_metric = scale_factor* (detSAi > fac ? fac : detSAi);
            shape_metric = scale_factor*detSAi;
            val_shape += shape_metric;
            //val_shape = std::min(val_shape, shape_metric);
            //std::cout << "tmp srk i= " << i << " detAi = " << detAi << " detSAi= " << detSAi << " shape_metric= " << shape_metric << " val_shape= " << val_shape << " scale_factor= " << scale_factor << " FAi= " << FAi << " FWi= " << FWi << std::endl;
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

      virtual double length_scaling_power() { return 1.0; }
      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA, jacW;

        int spatialDim = m_eMesh->get_spatial_dim();
        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;
        MsqMatrix<3,3> Ident; 
        identity(Ident);

        MsqMatrix<3,3> WI, T;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            MsqMatrix<3,3>& W = jacW.m_J[i];
            MsqMatrix<3,3>& A = jacA.m_J[i];

            // frob2 = h^2 + h^2 + 1
            // frob21 = 2 h^2
            // f = h sqrt(2)
            // det = h*h
            // met = f*f / (det*2) - 1
            // frob3 = 3 h^2
            // f = h sqrt(3)
            // det = h*h*h
            // met = f*f*f/(3^3/2 *det) - 1 = f*f*f/(3*sqrt(3)*det) - 1
            double shape_metric = 0.0;
            if (std::fabs(detAi) > 1.e-15)
              {
                inverse(W, WI);
                product(A, WI, T);
                double d = det(T);
                double f = my_sqr_Frobenius(T);
                if (spatialDim==2)
                  {
                    // all our jacobians are 3D, with a 1 in the 3,3 slot for 2d, so we subtract it here
                    f = f - 1.0;
                    f = std::sqrt(f);
                    double fac = 2.0;
                    double den = fac * d;
                    shape_metric = (f*f)/den - 1.0;
                  }
                else
                  {
                    f = std::sqrt(f);
                    double fac = 3.0*std::sqrt(3.0);
                    double den = fac * d;
                    shape_metric = (f*f*f)/den - 1.0;
                  }
                //shape_metric = std::fabs(shape_metric);
                //shape_metric = f/std::pow(den,1./3.) - 1.0;
              }
            val_shape += shape_metric;
            //val_shape += std::fabs(shape_metric);
            //val_shape += shape_metric*shape_metric;
          }
        val = val_shape;
        //val = val_shape*val_shape;
        return val;
      }

    };

    class PMMSmootherMetricLaplace : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricLaplace(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}

      virtual double length_scaling_power() { return 2.0; }
      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA;
        //JacobianUtil jacW;

        double A_ = 0.0;
        //double W_ = 0.0;
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            MsqMatrix<3,3>& A = jacA.m_J[i];
            double shape_metric = 0.0;
            shape_metric = sqr_Frobenius(A);
            val_shape += shape_metric;
          }
        val = val_shape;
        return val;
      }

    };

    class PMMSmootherMetricVolumetricEnergy : public PMMSmootherMetric
    {
    public:
      PMMSmootherMetricVolumetricEnergy(PerceptMesh *eMesh) : PMMSmootherMetric(eMesh) {}

      virtual double length_scaling_power() { return 6.0; }
      virtual double metric(stk::mesh::Entity& element, bool& valid)
      {
        valid = true;
        JacobianUtil jacA;
        //JacobianUtil jacW;

        double A_ = 0.0;
        //double W_ = 0.0;
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        //jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_shape=0.0;

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double detAi = jacA.m_detJ[i];
            if (detAi < 0)
              {
                valid = false;
              }
            double shape_metric = detAi*detAi;
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
