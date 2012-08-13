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
#include <stk_percept/math/Math.hpp>

namespace stk {
  namespace percept {

    class PMMSmootherMetric
    {
    public:
      PMMSmootherMetric(PerceptMesh *eMesh) : m_topology_data(0), m_eMesh(eMesh)
      {
        m_coord_field_current   = eMesh->get_coordinates_field();
        m_coord_field_original  = eMesh->get_field("coordinates_NM1");
      }

      virtual double metric(stk::mesh::Entity& element, bool& valid)=0;

      const CellTopologyData * m_topology_data ;

    protected:
      PerceptMesh *m_eMesh;
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
        //jacA.m_scale_to_unit = true;

        double A_ = 0.0, W_ = 0.0; // current and reference detJ
        jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
        jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
        double val=0.0, val_untangle=0.0;
        //double A_tot=0, W_tot=0;
        //MsqMatrix<3,3> ident; 
        //ident.identity();

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double Ai = jacA.mMetrics[i];
            double Wi = jacW.mMetrics[i];
            if (Ai <= 0.)
              {
                valid = false;
              }
            //A_tot += Ai;
            //W_tot += Wi;
            double untangle_metric = 0.0;
            double beta = m_beta_mult*Wi;
            double temp_var = Ai - beta;
            double fval=0.0;

            if(temp_var<0.0){
              //fval=std::fabs(temp_var)-temp_var;
              fval = -temp_var;
            }
            else
              {
                //fval = -0.001*temp_var;
              }

            //fval = Math::my_max_hi(-temp_var,0.0,beta*0.001);

            //fval = Ai;
            //untangle_metric = fval*fval;  // ensures that gradient is not zero without metric being zero
            untangle_metric = fval;
            val_untangle += untangle_metric;
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
        MsqMatrix<3,3> ident; 
        ident.identity();

        for (int i=0; i < jacA.m_num_nodes; i++)
          {
            double Ai = jacA.mMetrics[i];
            //double Wi = jacW.mMetrics[i];
            if (Ai < 0)
              {
                valid = false;
              }
            double shape_metric = 0.0;
            if (std::fabs(Ai) > 1.e-10)
              {
                shape_metric = sqr_Frobenius(jacW.mJ[i]*inverse(jacA.mJ[i]) - ident);
                //shape_metric = sqr_Frobenius(jacA.mJ[i]*inverse(jacW.mJ[i]) - ident);
              }
            //val_shape += std::sqrt(shape_metric);
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
