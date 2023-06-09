// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.



/* ========================================================================================== */
/*  This is generated code, do not edit - see SmootherMetricGen.mhpp and *.m, *.nb files      */
/* ========================================================================================== */


#ifndef SmootherMetricUntangleGen_hpp
#define SmootherMetricUntangleGen_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif // GCC_VERSION
#endif // __GNUC__

#ifdef __INTEL_COMPILER
#pragma warning disable 1599
#pragma warning disable 1478
#endif // __INTEL_COMPILER

#if defined (__clang__) && !defined(__INTEL_COMPILER)
#pragma clang diagnostic ignored "-Wunused-variable"
#endif

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/ScalingMatricesGen.hpp>
#include <percept/mesh/mod/smoother/CodeGenHelper.hpp>  // for MyHeaviside, MyPow, etc.

namespace percept {

  class SmootherMetricUntangleGen : public SmootherMetricUntangle
  {
    MDArray A, W, WI;
    stk::mesh::FieldBase *m_cg_lambda_field;
    stk::mesh::FieldBase *m_cg_normal_field;
    stk::mesh::FieldBase *m_coord_field_0;
  public:
    SmootherMetricUntangleGen(PerceptMesh *eMesh) : SmootherMetricUntangle(eMesh) {
      m_cg_lambda_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_lambda");
      m_cg_normal_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_normal");
      m_coord_field_0      = eMesh->get_field(stk::topology::NODE_RANK, "coordinates_0");
      A.resize(3,3);
      W.resize(3,3);
      WI.resize(3,3);
    }

    virtual bool has_gradient() { return true; }
    virtual bool has_gradient_and_hessian() { return true; }

    virtual double metric(stk::mesh::Entity element, bool& valid)
    {
      static double grad[8][4];
      static double hess[8][4][8][4];
      return grad_and_hessian(element, valid, grad, hess, false, false);
    }

    /// computes metric and its gradient - see Mesquite::TShapeB1, TQualityMetric, TargetMetricUtil
    virtual double grad_metric(stk::mesh::Entity element, bool& valid, double grad[8][4])
    {
      double hess[8][4][8][4];
      return grad_and_hessian(element, valid, grad, hess, true, false);
    }

    virtual double grad_and_hessian(stk::mesh::Entity element, bool& valid, double grad[8][4], double hess[8][4][8][4], const bool getGrad=false, const bool getHess=false)
    {
      const double beta = m_beta_mult;
      valid = true;

      stk::mesh::Entity const *v_i = m_eMesh->get_bulk_data()->begin_nodes(element);
      size_t num_nodes = m_eMesh->get_bulk_data()->num_nodes(element);

#undef AVERTEX
#define AVERTEX(vi)  static_cast<double*>(stk::mesh::field_data( *m_coord_field_current, vi ))
#undef WVERTEX
#define WVERTEX(vi)  static_cast<double*>(stk::mesh::field_data( *m_coord_field_original, vi ))

      // double A_ = 0.0, W_ = 0.0; // current and reference detJ
      // jacA(A_, *m_eMesh, element, m_coord_field_current, m_topology_data);
      // jacW(W_, *m_eMesh, element, m_coord_field_original, m_topology_data);
      double val=0.0, val_metric=0.0;

      if (getGrad) memset(&grad[0][0],0,8*4*sizeof(double));
      if (getHess) memset(&hess[0][0][0][0],0,8*4*8*4*sizeof(double));

      stk::topology::topology_t stk_topo = m_eMesh->topology(element);
      const MDArray& Sc = ScalingMatrices::s_scalingMatrices.get(stk_topo);
      
       const double Sc00 = Sc(0,0) ;

       const double Sc01 = Sc(0,1) ;

       const double Sc02 = Sc(0,2) ;

       const double Sc10 = Sc(1,0) ;

       const double Sc11 = Sc(1,1) ;

       const double Sc12 = Sc(1,2) ;

       const double Sc20 = Sc(2,0) ;

       const double Sc21 = Sc(2,1) ;

       const double Sc22 = Sc(2,2) ;

      double *x0, *x1, *x2, *x3, *y0;
      double *wx0, *wx1, *wx2, *wx3;
      double *X[8];
      double *WX[8];
      double L[8];
      double *Y[8];
      double  Lambda0,  Lambda1,  Lambda2,  Lambda3;
      double onB[8];
      double *Norm[8];

      int spatialDim = m_eMesh->get_spatial_dim();

      static std::vector<double> normals(3, 0.0);
      for (int i=0; i < num_nodes; i++)
        {
          Norm[i] = (m_cg_normal_field? static_cast<double*>(stk::mesh::field_data( *m_cg_normal_field, v_i[i] )) : 0);

          X[i] = AVERTEX(v_i[i]);
          WX[i] = WVERTEX(v_i[i]);
          L[i] = (m_cg_lambda_field? *static_cast<double*>(stk::mesh::field_data( *m_cg_lambda_field, v_i[i] )) : 0.0);
          Y[i] = (m_coord_field_0 ? static_cast<double*>(stk::mesh::field_data( *m_coord_field_0, v_i[i] )) : 0);
          std::pair<bool,int> fixed;
          onB[i] = 0.0;
        }



#define X0(i) (spatialDim == 3 ? x0[i] : i == 2 ? 0 : x0[i])
#define X1(i) (spatialDim == 3 ? x1[i] : i == 2 ? 0 : x1[i])
#define X2(i) (spatialDim == 3 ? x2[i] : i == 2 ? 0 : x2[i])
#define X3(i) (spatialDim == 3 ? x3[i] : i == 2 ? 1 : 0)

#define Y0(i) (spatialDim == 3 ? y0[i] : i == 2 ? 0 : y0[i])

#define WX0(i) (spatialDim == 3 ? wx0[i] : i == 2 ? 0 : wx0[i])
#define WX1(i) (spatialDim == 3 ? wx1[i] : i == 2 ? 0 : wx1[i])
#define WX2(i) (spatialDim == 3 ? wx2[i] : i == 2 ? 0 : wx2[i])
#define WX3(i) (spatialDim == 3 ? wx3[i] : i == 2 ? 1 : 0)

#define II(i) indices[i]
#define GRAD(i,j) grad[i][j]
#define HESS(i,j,k,l) hess[i][j][k][l]

#define normal(j) (spatialDim == 3 ? Norm[i][j] : j == 2 ? 0 : Norm[i][j])

      for (int i=0; i < num_nodes; i++)
        {

          const int *indices = Indices::s_indices.get_indices(stk_topo, i);
          x0 = X[indices[0]];
          x1 = X[indices[1]];
          x2 = X[indices[2]];
          x3 = X[indices[3]];
          wx0 = WX[indices[0]];
          wx1 = WX[indices[1]];
          wx2 = WX[indices[2]];
          wx3 = WX[indices[3]];

          Lambda0 = L[indices[0]];
          Lambda1 = L[indices[1]];
          Lambda2 = L[indices[2]];
          Lambda3 = L[indices[3]];

          y0 = Y[indices[0]];

          double onBoundary = onB[indices[0]];

          double vv = 0.0, sdetA = 0.0, sdetW = 0.0;
          if (!getGrad && !getHess)
            {
              
       const double A00 = -X0(0) + X1(0) ;

       const double A01 = -X0(0) + X2(0) ;

       const double A02 = -X0(0) + X3(0) ;

       const double A10 = -X0(1) + X1(1) ;

       const double A11 = -X0(1) + X2(1) ;

       const double A12 = -X0(1) + X3(1) ;

       const double A20 = -X0(2) + X1(2) ;

       const double A21 = -X0(2) + X2(2) ;

       const double A22 = -X0(2) + X3(2) ;

       const double W00 = -WX0(0) + WX1(0) ;

       const double W01 = -WX0(0) + WX2(0) ;

       const double W02 = -WX0(0) + WX3(0) ;

       const double W10 = -WX0(1) + WX1(1) ;

       const double W11 = -WX0(1) + WX2(1) ;

       const double W12 = -WX0(1) + WX3(1) ;

       const double W20 = -WX0(2) + WX1(2) ;

       const double W21 = -WX0(2) + WX2(2) ;

       const double W22 = -WX0(2) + WX3(2) ;

       const double detA = A02*(-(A11*A20) + A10*A21) + A01*(A12*A20 - A10*A22) + A00*(-(A12*A21) + A11*A22) ;

       const double detW = -(W02*W11*W20) + W01*W12*W20 + W02*W10*W21 - W00*W12*W21 - W01*W10*W22 + W00*W11*W22 ;

       //const double detWI = MyInverse(detW) ;

       //const double WI00 = detWI*(-(W12*W21) + W11*W22) ;

       //const double WI01 = detWI*(W02*W21 - W01*W22) ;

       //const double WI02 = detWI*(-(W02*W11) + W01*W12) ;

       //const double WI10 = detWI*(W12*W20 - W10*W22) ;

       //const double WI11 = detWI*(-(W02*W20) + W00*W22) ;

       //const double WI12 = detWI*(W02*W10 - W00*W12) ;

       //const double WI20 = detWI*(-(W11*W20) + W10*W21) ;

       //const double WI21 = detWI*(W01*W20 - W00*W21) ;

       //const double WI22 = detWI*(-(W01*W10) + W00*W11) ;

       //const double T00 = A00*WI00 + A01*WI10 + A02*WI20 ;

       //const double T01 = A00*WI01 + A01*WI11 + A02*WI21 ;

       //const double T02 = A00*WI02 + A01*WI12 + A02*WI22 ;

       //const double T10 = A10*WI00 + A11*WI10 + A12*WI20 ;

       //const double T11 = A10*WI01 + A11*WI11 + A12*WI21 ;

       //const double T12 = A10*WI02 + A11*WI12 + A12*WI22 ;

       //const double T20 = A20*WI00 + A21*WI10 + A22*WI20 ;

       //const double T21 = A20*WI01 + A21*WI11 + A22*WI21 ;

       //const double T22 = A20*WI02 + A21*WI12 + A22*WI22 ;

       const double vU = -detA + beta*detW ;

       const double hv = MyHeaviside(vU,0) ;

       const double met = hv*MyPow2(vU) ;

#ifndef NDEBUG
              if (detW <= 0.0)
                {
                  std::cout << "detW= 0, i = " << i << " detA= " << detA << " m_topology_data= " << m_topology_data << std::endl;
                  m_eMesh->print_entity(element);
                  m_eMesh->print_entity(element, m_coord_field_original);
                  if (m_topology_data)
                    {
                      shards::CellTopology topology(m_topology_data);
                      std::cout << "topology = " << topology.getName() << std::endl;
                    }
                }
#endif
              vv = met; sdetA = detA; sdetW = detW;
            }
          else
            {
              if (getGrad && !getHess)
                {
                  
       const double A00 = -X0(0) + X1(0) ;

       const double A01 = -X0(0) + X2(0) ;

       const double A02 = -X0(0) + X3(0) ;

       const double A10 = -X0(1) + X1(1) ;

       const double A11 = -X0(1) + X2(1) ;

       const double A12 = -X0(1) + X3(1) ;

       const double A20 = -X0(2) + X1(2) ;

       const double A21 = -X0(2) + X2(2) ;

       const double A22 = -X0(2) + X3(2) ;

       const double W00 = -WX0(0) + WX1(0) ;

       const double W01 = -WX0(0) + WX2(0) ;

       const double W02 = -WX0(0) + WX3(0) ;

       const double W10 = -WX0(1) + WX1(1) ;

       const double W11 = -WX0(1) + WX2(1) ;

       const double W12 = -WX0(1) + WX3(1) ;

       const double W20 = -WX0(2) + WX1(2) ;

       const double W21 = -WX0(2) + WX2(2) ;

       const double W22 = -WX0(2) + WX3(2) ;

       const double detA = A02*(-(A11*A20) + A10*A21) + A01*(A12*A20 - A10*A22) + A00*(-(A12*A21) + A11*A22) ;

       const double detW = -(W02*W11*W20) + W01*W12*W20 + W02*W10*W21 - W00*W12*W21 - W01*W10*W22 + W00*W11*W22 ;

       const double detWI = MyInverse(detW) ;

       const double WI00 = detWI*(-(W12*W21) + W11*W22) ;

       const double WI01 = detWI*(W02*W21 - W01*W22) ;

       const double WI02 = detWI*(-(W02*W11) + W01*W12) ;

       const double WI10 = detWI*(W12*W20 - W10*W22) ;

       const double WI11 = detWI*(-(W02*W20) + W00*W22) ;

       const double WI12 = detWI*(W02*W10 - W00*W12) ;

       const double WI20 = detWI*(-(W11*W20) + W10*W21) ;

       const double WI21 = detWI*(W01*W20 - W00*W21) ;

       const double WI22 = detWI*(-(W01*W10) + W00*W11) ;

       const double T00 = A00*WI00 + A01*WI10 + A02*WI20 ;

       const double T01 = A00*WI01 + A01*WI11 + A02*WI21 ;

       const double T02 = A00*WI02 + A01*WI12 + A02*WI22 ;

       const double T10 = A10*WI00 + A11*WI10 + A12*WI20 ;

       const double T11 = A10*WI01 + A11*WI11 + A12*WI21 ;

       const double T12 = A10*WI02 + A11*WI12 + A12*WI22 ;

       const double T20 = A20*WI00 + A21*WI10 + A22*WI20 ;

       const double T21 = A20*WI01 + A21*WI11 + A22*WI21 ;

       const double T22 = A20*WI02 + A21*WI12 + A22*WI22 ;

       const double vU = -detA + beta*detW ;

       const double hv = MyHeaviside(vU,0) ;

       const double met = hv*MyPow2(vU) ;

       const double D_A00_100000000000 = -1 ;

       const double D_A00_000100000000 = 1 ;

       const double D_A01_100000000000 = -1 ;

       const double D_A01_000000100000 = 1 ;

       const double D_A02_100000000000 = -1 ;

       const double D_A02_000000000100 = 1 ;

       const double D_A10_010000000000 = -1 ;

       const double D_A10_000010000000 = 1 ;

       const double D_A11_010000000000 = -1 ;

       const double D_A11_000000010000 = 1 ;

       const double D_A12_010000000000 = -1 ;

       const double D_A12_000000000010 = 1 ;

       const double D_A20_001000000000 = -1 ;

       const double D_A20_000001000000 = 1 ;

       const double D_A21_001000000000 = -1 ;

       const double D_A21_000000001000 = 1 ;

       const double D_A22_001000000000 = -1 ;

       const double D_A22_000000000001 = 1 ;

       const double D_detA_100000000000 = (-(A12*A21) + A11*A22)*D_A00_100000000000 + (A12*A20 -\
 
  A10*A22)*D_A01_100000000000 + (-(A11*A20) + A10*A21)*D_A02_100000000000 ;

       const double D_detA_010000000000 = A02*(A21*D_A10_010000000000 - A20*D_A11_010000000000) +\
 
  A01*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) + A00*(A22*D_A11_010000000000 - A21*D_A12_010000000000) ;

       const double D_detA_001000000000 = A02*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) +\
 
  A01*(A12*D_A20_001000000000 - A10*D_A22_001000000000) + A00*(-(A12*D_A21_001000000000) + A11*D_A22_001000000000) ;

       const double D_detA_000100000000 = (-(A12*A21) + A11*A22)*D_A00_000100000000 ;

       const double D_detA_000010000000 = A02*A21*D_A10_000010000000 - A01*A22*D_A10_000010000000 ;

       const double D_detA_000001000000 = -(A02*A11*D_A20_000001000000) + A01*A12*D_A20_000001000000 ;

       const double D_detA_000000100000 = (A12*A20 - A10*A22)*D_A01_000000100000 ;

       const double D_detA_000000010000 = -(A02*A20*D_A11_000000010000) + A00*A22*D_A11_000000010000 ;

       const double D_detA_000000001000 = A02*A10*D_A21_000000001000 - A00*A12*D_A21_000000001000 ;

       const double D_detA_000000000100 = (-(A11*A20) + A10*A21)*D_A02_000000000100 ;

       const double D_detA_000000000010 = A01*A20*D_A12_000000000010 - A00*A21*D_A12_000000000010 ;

       const double D_detA_000000000001 = -(A01*A10*D_A22_000000000001) + A00*A11*D_A22_000000000001 ;

       const double D_T00_100000000000 = D_A00_100000000000*WI00 + D_A01_100000000000*WI10 + D_A02_100000000000*WI20 ;

       const double D_T00_000100000000 = D_A00_000100000000*WI00 ;

       const double D_T00_000000100000 = D_A01_000000100000*WI10 ;

       const double D_T00_000000000100 = D_A02_000000000100*WI20 ;

       const double D_T01_100000000000 = D_A00_100000000000*WI01 + D_A01_100000000000*WI11 + D_A02_100000000000*WI21 ;

       const double D_T01_000100000000 = D_A00_000100000000*WI01 ;

       const double D_T01_000000100000 = D_A01_000000100000*WI11 ;

       const double D_T01_000000000100 = D_A02_000000000100*WI21 ;

       const double D_T02_100000000000 = D_A00_100000000000*WI02 + D_A01_100000000000*WI12 + D_A02_100000000000*WI22 ;

       const double D_T02_000100000000 = D_A00_000100000000*WI02 ;

       const double D_T02_000000100000 = D_A01_000000100000*WI12 ;

       const double D_T02_000000000100 = D_A02_000000000100*WI22 ;

       const double D_T10_010000000000 = D_A10_010000000000*WI00 + D_A11_010000000000*WI10 + D_A12_010000000000*WI20 ;

       const double D_T10_000010000000 = D_A10_000010000000*WI00 ;

       const double D_T10_000000010000 = D_A11_000000010000*WI10 ;

       const double D_T10_000000000010 = D_A12_000000000010*WI20 ;

       const double D_T11_010000000000 = D_A10_010000000000*WI01 + D_A11_010000000000*WI11 + D_A12_010000000000*WI21 ;

       const double D_T11_000010000000 = D_A10_000010000000*WI01 ;

       const double D_T11_000000010000 = D_A11_000000010000*WI11 ;

       const double D_T11_000000000010 = D_A12_000000000010*WI21 ;

       const double D_T12_010000000000 = D_A10_010000000000*WI02 + D_A11_010000000000*WI12 + D_A12_010000000000*WI22 ;

       const double D_T12_000010000000 = D_A10_000010000000*WI02 ;

       const double D_T12_000000010000 = D_A11_000000010000*WI12 ;

       const double D_T12_000000000010 = D_A12_000000000010*WI22 ;

       const double D_T20_001000000000 = D_A20_001000000000*WI00 + D_A21_001000000000*WI10 + D_A22_001000000000*WI20 ;

       const double D_T20_000001000000 = D_A20_000001000000*WI00 ;

       const double D_T20_000000001000 = D_A21_000000001000*WI10 ;

       const double D_T20_000000000001 = D_A22_000000000001*WI20 ;

       const double D_T21_001000000000 = D_A20_001000000000*WI01 + D_A21_001000000000*WI11 + D_A22_001000000000*WI21 ;

       const double D_T21_000001000000 = D_A20_000001000000*WI01 ;

       const double D_T21_000000001000 = D_A21_000000001000*WI11 ;

       const double D_T21_000000000001 = D_A22_000000000001*WI21 ;

       const double D_T22_001000000000 = D_A20_001000000000*WI02 + D_A21_001000000000*WI12 + D_A22_001000000000*WI22 ;

       const double D_T22_000001000000 = D_A20_000001000000*WI02 ;

       const double D_T22_000000001000 = D_A21_000000001000*WI12 ;

       const double D_T22_000000000001 = D_A22_000000000001*WI22 ;

       const double D_vU_100000000000 = -D_detA_100000000000 ;

       const double D_vU_010000000000 = -D_detA_010000000000 ;

       const double D_vU_001000000000 = -D_detA_001000000000 ;

       const double D_vU_000100000000 = -D_detA_000100000000 ;

       const double D_vU_000010000000 = -D_detA_000010000000 ;

       const double D_vU_000001000000 = -D_detA_000001000000 ;

       const double D_vU_000000100000 = -D_detA_000000100000 ;

       const double D_vU_000000010000 = -D_detA_000000010000 ;

       const double D_vU_000000001000 = -D_detA_000000001000 ;

       const double D_vU_000000000100 = -D_detA_000000000100 ;

       const double D_vU_000000000010 = -D_detA_000000000010 ;

       const double D_vU_000000000001 = -D_detA_000000000001 ;

       const double D_met_100000000000 = 2*D_vU_100000000000*hv*vU ;

       const double D_met_010000000000 = 2*D_vU_010000000000*hv*vU ;

       const double D_met_001000000000 = 2*D_vU_001000000000*hv*vU ;

       const double D_met_000100000000 = 2*D_vU_000100000000*hv*vU ;

       const double D_met_000010000000 = 2*D_vU_000010000000*hv*vU ;

       const double D_met_000001000000 = 2*D_vU_000001000000*hv*vU ;

       const double D_met_000000100000 = 2*D_vU_000000100000*hv*vU ;

       const double D_met_000000010000 = 2*D_vU_000000010000*hv*vU ;

       const double D_met_000000001000 = 2*D_vU_000000001000*hv*vU ;

       const double D_met_000000000100 = 2*D_vU_000000000100*hv*vU ;

       const double D_met_000000000010 = 2*D_vU_000000000010*hv*vU ;

       const double D_met_000000000001 = 2*D_vU_000000000001*hv*vU ;
                  
       GRAD(II(0),0) += D_met_100000000000 ;

       GRAD(II(0),1) += D_met_010000000000 ;

       GRAD(II(0),2) += D_met_001000000000 ;

       GRAD(II(1),0) += D_met_000100000000 ;

       GRAD(II(1),1) += D_met_000010000000 ;

       GRAD(II(1),2) += D_met_000001000000 ;

       GRAD(II(2),0) += D_met_000000100000 ;

       GRAD(II(2),1) += D_met_000000010000 ;

       GRAD(II(2),2) += D_met_000000001000 ;

       GRAD(II(3),0) += D_met_000000000100 ;

       GRAD(II(3),1) += D_met_000000000010 ;

       GRAD(II(3),2) += D_met_000000000001 ;
                  vv = met; sdetA = detA; sdetW = detW;
                }
              else if (!getGrad && getHess)
                {
                  
       const double A00 = -X0(0) + X1(0) ;

       const double A01 = -X0(0) + X2(0) ;

       const double A02 = -X0(0) + X3(0) ;

       const double A10 = -X0(1) + X1(1) ;

       const double A11 = -X0(1) + X2(1) ;

       const double A12 = -X0(1) + X3(1) ;

       const double A20 = -X0(2) + X1(2) ;

       const double A21 = -X0(2) + X2(2) ;

       const double A22 = -X0(2) + X3(2) ;

       const double W00 = -WX0(0) + WX1(0) ;

       const double W01 = -WX0(0) + WX2(0) ;

       const double W02 = -WX0(0) + WX3(0) ;

       const double W10 = -WX0(1) + WX1(1) ;

       const double W11 = -WX0(1) + WX2(1) ;

       const double W12 = -WX0(1) + WX3(1) ;

       const double W20 = -WX0(2) + WX1(2) ;

       const double W21 = -WX0(2) + WX2(2) ;

       const double W22 = -WX0(2) + WX3(2) ;

       const double detA = A02*(-(A11*A20) + A10*A21) + A01*(A12*A20 - A10*A22) + A00*(-(A12*A21) + A11*A22) ;

       const double detW = -(W02*W11*W20) + W01*W12*W20 + W02*W10*W21 - W00*W12*W21 - W01*W10*W22 + W00*W11*W22 ;

       const double detWI = MyInverse(detW) ;

       const double WI00 = detWI*(-(W12*W21) + W11*W22) ;

       const double WI01 = detWI*(W02*W21 - W01*W22) ;

       const double WI02 = detWI*(-(W02*W11) + W01*W12) ;

       const double WI10 = detWI*(W12*W20 - W10*W22) ;

       const double WI11 = detWI*(-(W02*W20) + W00*W22) ;

       const double WI12 = detWI*(W02*W10 - W00*W12) ;

       const double WI20 = detWI*(-(W11*W20) + W10*W21) ;

       const double WI21 = detWI*(W01*W20 - W00*W21) ;

       const double WI22 = detWI*(-(W01*W10) + W00*W11) ;

       const double T00 = A00*WI00 + A01*WI10 + A02*WI20 ;

       const double T01 = A00*WI01 + A01*WI11 + A02*WI21 ;

       const double T02 = A00*WI02 + A01*WI12 + A02*WI22 ;

       const double T10 = A10*WI00 + A11*WI10 + A12*WI20 ;

       const double T11 = A10*WI01 + A11*WI11 + A12*WI21 ;

       const double T12 = A10*WI02 + A11*WI12 + A12*WI22 ;

       const double T20 = A20*WI00 + A21*WI10 + A22*WI20 ;

       const double T21 = A20*WI01 + A21*WI11 + A22*WI21 ;

       const double T22 = A20*WI02 + A21*WI12 + A22*WI22 ;

       const double vU = -detA + beta*detW ;

       const double hv = MyHeaviside(vU,0) ;

       const double met = hv*MyPow2(vU) ;

       const double D_A00_100000000000 = -1 ;

       const double D_A00_000100000000 = 1 ;

       const double D_A01_100000000000 = -1 ;

       const double D_A01_000000100000 = 1 ;

       const double D_A02_100000000000 = -1 ;

       const double D_A02_000000000100 = 1 ;

       const double D_A10_010000000000 = -1 ;

       const double D_A10_000010000000 = 1 ;

       const double D_A11_010000000000 = -1 ;

       const double D_A11_000000010000 = 1 ;

       const double D_A12_010000000000 = -1 ;

       const double D_A12_000000000010 = 1 ;

       const double D_A20_001000000000 = -1 ;

       const double D_A20_000001000000 = 1 ;

       const double D_A21_001000000000 = -1 ;

       const double D_A21_000000001000 = 1 ;

       const double D_A22_001000000000 = -1 ;

       const double D_A22_000000000001 = 1 ;

       const double D_detA_100000000000 = (-(A12*A21) + A11*A22)*D_A00_100000000000 + (A12*A20 -\
 
  A10*A22)*D_A01_100000000000 + (-(A11*A20) + A10*A21)*D_A02_100000000000 ;

       const double D_detA_010000000000 = A02*(A21*D_A10_010000000000 - A20*D_A11_010000000000) +\
 
  A01*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) + A00*(A22*D_A11_010000000000 - A21*D_A12_010000000000) ;

       const double D_detA_001000000000 = A02*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) +\
 
  A01*(A12*D_A20_001000000000 - A10*D_A22_001000000000) + A00*(-(A12*D_A21_001000000000) + A11*D_A22_001000000000) ;

       const double D_detA_000100000000 = (-(A12*A21) + A11*A22)*D_A00_000100000000 ;

       const double D_detA_000010000000 = A02*A21*D_A10_000010000000 - A01*A22*D_A10_000010000000 ;

       const double D_detA_000001000000 = -(A02*A11*D_A20_000001000000) + A01*A12*D_A20_000001000000 ;

       const double D_detA_000000100000 = (A12*A20 - A10*A22)*D_A01_000000100000 ;

       const double D_detA_000000010000 = -(A02*A20*D_A11_000000010000) + A00*A22*D_A11_000000010000 ;

       const double D_detA_000000001000 = A02*A10*D_A21_000000001000 - A00*A12*D_A21_000000001000 ;

       const double D_detA_000000000100 = (-(A11*A20) + A10*A21)*D_A02_000000000100 ;

       const double D_detA_000000000010 = A01*A20*D_A12_000000000010 - A00*A21*D_A12_000000000010 ;

       const double D_detA_000000000001 = -(A01*A10*D_A22_000000000001) + A00*A11*D_A22_000000000001 ;

       const double D_T00_100000000000 = D_A00_100000000000*WI00 + D_A01_100000000000*WI10 + D_A02_100000000000*WI20 ;

       const double D_T00_000100000000 = D_A00_000100000000*WI00 ;

       const double D_T00_000000100000 = D_A01_000000100000*WI10 ;

       const double D_T00_000000000100 = D_A02_000000000100*WI20 ;

       const double D_T01_100000000000 = D_A00_100000000000*WI01 + D_A01_100000000000*WI11 + D_A02_100000000000*WI21 ;

       const double D_T01_000100000000 = D_A00_000100000000*WI01 ;

       const double D_T01_000000100000 = D_A01_000000100000*WI11 ;

       const double D_T01_000000000100 = D_A02_000000000100*WI21 ;

       const double D_T02_100000000000 = D_A00_100000000000*WI02 + D_A01_100000000000*WI12 + D_A02_100000000000*WI22 ;

       const double D_T02_000100000000 = D_A00_000100000000*WI02 ;

       const double D_T02_000000100000 = D_A01_000000100000*WI12 ;

       const double D_T02_000000000100 = D_A02_000000000100*WI22 ;

       const double D_T10_010000000000 = D_A10_010000000000*WI00 + D_A11_010000000000*WI10 + D_A12_010000000000*WI20 ;

       const double D_T10_000010000000 = D_A10_000010000000*WI00 ;

       const double D_T10_000000010000 = D_A11_000000010000*WI10 ;

       const double D_T10_000000000010 = D_A12_000000000010*WI20 ;

       const double D_T11_010000000000 = D_A10_010000000000*WI01 + D_A11_010000000000*WI11 + D_A12_010000000000*WI21 ;

       const double D_T11_000010000000 = D_A10_000010000000*WI01 ;

       const double D_T11_000000010000 = D_A11_000000010000*WI11 ;

       const double D_T11_000000000010 = D_A12_000000000010*WI21 ;

       const double D_T12_010000000000 = D_A10_010000000000*WI02 + D_A11_010000000000*WI12 + D_A12_010000000000*WI22 ;

       const double D_T12_000010000000 = D_A10_000010000000*WI02 ;

       const double D_T12_000000010000 = D_A11_000000010000*WI12 ;

       const double D_T12_000000000010 = D_A12_000000000010*WI22 ;

       const double D_T20_001000000000 = D_A20_001000000000*WI00 + D_A21_001000000000*WI10 + D_A22_001000000000*WI20 ;

       const double D_T20_000001000000 = D_A20_000001000000*WI00 ;

       const double D_T20_000000001000 = D_A21_000000001000*WI10 ;

       const double D_T20_000000000001 = D_A22_000000000001*WI20 ;

       const double D_T21_001000000000 = D_A20_001000000000*WI01 + D_A21_001000000000*WI11 + D_A22_001000000000*WI21 ;

       const double D_T21_000001000000 = D_A20_000001000000*WI01 ;

       const double D_T21_000000001000 = D_A21_000000001000*WI11 ;

       const double D_T21_000000000001 = D_A22_000000000001*WI21 ;

       const double D_T22_001000000000 = D_A20_001000000000*WI02 + D_A21_001000000000*WI12 + D_A22_001000000000*WI22 ;

       const double D_T22_000001000000 = D_A20_000001000000*WI02 ;

       const double D_T22_000000001000 = D_A21_000000001000*WI12 ;

       const double D_T22_000000000001 = D_A22_000000000001*WI22 ;

       const double D_vU_100000000000 = -D_detA_100000000000 ;

       const double D_vU_010000000000 = -D_detA_010000000000 ;

       const double D_vU_001000000000 = -D_detA_001000000000 ;

       const double D_vU_000100000000 = -D_detA_000100000000 ;

       const double D_vU_000010000000 = -D_detA_000010000000 ;

       const double D_vU_000001000000 = -D_detA_000001000000 ;

       const double D_vU_000000100000 = -D_detA_000000100000 ;

       const double D_vU_000000010000 = -D_detA_000000010000 ;

       const double D_vU_000000001000 = -D_detA_000000001000 ;

       const double D_vU_000000000100 = -D_detA_000000000100 ;

       const double D_vU_000000000010 = -D_detA_000000000010 ;

       const double D_vU_000000000001 = -D_detA_000000000001 ;

       const double D_met_100000000000 = 2*D_vU_100000000000*hv*vU ;

       const double D_met_010000000000 = 2*D_vU_010000000000*hv*vU ;

       const double D_met_001000000000 = 2*D_vU_001000000000*hv*vU ;

       const double D_met_000100000000 = 2*D_vU_000100000000*hv*vU ;

       const double D_met_000010000000 = 2*D_vU_000010000000*hv*vU ;

       const double D_met_000001000000 = 2*D_vU_000001000000*hv*vU ;

       const double D_met_000000100000 = 2*D_vU_000000100000*hv*vU ;

       const double D_met_000000010000 = 2*D_vU_000000010000*hv*vU ;

       const double D_met_000000001000 = 2*D_vU_000000001000*hv*vU ;

       const double D_met_000000000100 = 2*D_vU_000000000100*hv*vU ;

       const double D_met_000000000010 = 2*D_vU_000000000010*hv*vU ;

       const double D_met_000000000001 = 2*D_vU_000000000001*hv*vU ;

       const double D_detA_110000000000 = D_A02_100000000000*(A21*D_A10_010000000000 - A20*D_A11_010000000000) +\
 
  D_A01_100000000000*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) + D_A00_100000000000*(A22*D_A11_010000000000\
 
  - A21*D_A12_010000000000) ;

       const double D_detA_101000000000 = D_A02_100000000000*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) +\
 
  D_A01_100000000000*(A12*D_A20_001000000000 - A10*D_A22_001000000000) + D_A00_100000000000*(-(A12*D_A21_001000000000)\
 
  + A11*D_A22_001000000000) ;

       const double D_detA_100010000000 = -(A22*D_A01_100000000000*D_A10_000010000000) +\
 
  A21*D_A02_100000000000*D_A10_000010000000 ;

       const double D_detA_100001000000 = A12*D_A01_100000000000*D_A20_000001000000 -\
 
  A11*D_A02_100000000000*D_A20_000001000000 ;

       const double D_detA_100000010000 = A22*D_A00_100000000000*D_A11_000000010000 -\
 
  A20*D_A02_100000000000*D_A11_000000010000 ;

       const double D_detA_100000001000 = -(A12*D_A00_100000000000*D_A21_000000001000) +\
 
  A10*D_A02_100000000000*D_A21_000000001000 ;

       const double D_detA_100000000010 = -(A21*D_A00_100000000000*D_A12_000000000010) +\
 
  A20*D_A01_100000000000*D_A12_000000000010 ;

       const double D_detA_100000000001 = A11*D_A00_100000000000*D_A22_000000000001 -\
 
  A10*D_A01_100000000000*D_A22_000000000001 ;

       const double D_detA_011000000000 = A02*(-(D_A11_010000000000*D_A20_001000000000) +\
 
  D_A10_010000000000*D_A21_001000000000) + A01*(D_A12_010000000000*D_A20_001000000000 -\
 
  D_A10_010000000000*D_A22_001000000000) + A00*(-(D_A12_010000000000*D_A21_001000000000) +\
 
  D_A11_010000000000*D_A22_001000000000) ;

       const double D_detA_010100000000 = D_A00_000100000000*(A22*D_A11_010000000000 - A21*D_A12_010000000000) ;

       const double D_detA_010001000000 = -(A02*D_A11_010000000000*D_A20_000001000000) +\
 
  A01*D_A12_010000000000*D_A20_000001000000 ;

       const double D_detA_010000100000 = D_A01_000000100000*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) ;

       const double D_detA_010000001000 = A02*D_A10_010000000000*D_A21_000000001000 -\
 
  A00*D_A12_010000000000*D_A21_000000001000 ;

       const double D_detA_010000000100 = D_A02_000000000100*(A21*D_A10_010000000000 - A20*D_A11_010000000000) ;

       const double D_detA_010000000001 = -(A01*D_A10_010000000000*D_A22_000000000001) +\
 
  A00*D_A11_010000000000*D_A22_000000000001 ;

       const double D_detA_001100000000 = D_A00_000100000000*(-(A12*D_A21_001000000000) + A11*D_A22_001000000000) ;

       const double D_detA_001010000000 = A02*D_A10_000010000000*D_A21_001000000000 -\
 
  A01*D_A10_000010000000*D_A22_001000000000 ;

       const double D_detA_001000100000 = D_A01_000000100000*(A12*D_A20_001000000000 - A10*D_A22_001000000000) ;

       const double D_detA_001000010000 = -(A02*D_A11_000000010000*D_A20_001000000000) +\
 
  A00*D_A11_000000010000*D_A22_001000000000 ;

       const double D_detA_001000000100 = D_A02_000000000100*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) ;

       const double D_detA_001000000010 = A01*D_A12_000000000010*D_A20_001000000000 -\
 
  A00*D_A12_000000000010*D_A21_001000000000 ;

       const double D_detA_000100010000 = A22*D_A00_000100000000*D_A11_000000010000 ;

       const double D_detA_000100001000 = -(A12*D_A00_000100000000*D_A21_000000001000) ;

       const double D_detA_000100000010 = -(A21*D_A00_000100000000*D_A12_000000000010) ;

       const double D_detA_000100000001 = A11*D_A00_000100000000*D_A22_000000000001 ;

       const double D_detA_000010100000 = -(A22*D_A01_000000100000*D_A10_000010000000) ;

       const double D_detA_000010001000 = A02*D_A10_000010000000*D_A21_000000001000 ;

       const double D_detA_000010000100 = A21*D_A02_000000000100*D_A10_000010000000 ;

       const double D_detA_000010000001 = -(A01*D_A10_000010000000*D_A22_000000000001) ;

       const double D_detA_000001100000 = A12*D_A01_000000100000*D_A20_000001000000 ;

       const double D_detA_000001010000 = -(A02*D_A11_000000010000*D_A20_000001000000) ;

       const double D_detA_000001000100 = -(A11*D_A02_000000000100*D_A20_000001000000) ;

       const double D_detA_000001000010 = A01*D_A12_000000000010*D_A20_000001000000 ;

       const double D_detA_000000100010 = A20*D_A01_000000100000*D_A12_000000000010 ;

       const double D_detA_000000100001 = -(A10*D_A01_000000100000*D_A22_000000000001) ;

       const double D_detA_000000010100 = -(A20*D_A02_000000000100*D_A11_000000010000) ;

       const double D_detA_000000010001 = A00*D_A11_000000010000*D_A22_000000000001 ;

       const double D_detA_000000001100 = A10*D_A02_000000000100*D_A21_000000001000 ;

       const double D_detA_000000001010 = -(A00*D_A12_000000000010*D_A21_000000001000) ;

       const double D_vU_110000000000 = -D_detA_110000000000 ;

       const double D_vU_101000000000 = -D_detA_101000000000 ;

       const double D_vU_100010000000 = -D_detA_100010000000 ;

       const double D_vU_100001000000 = -D_detA_100001000000 ;

       const double D_vU_100000010000 = -D_detA_100000010000 ;

       const double D_vU_100000001000 = -D_detA_100000001000 ;

       const double D_vU_100000000010 = -D_detA_100000000010 ;

       const double D_vU_100000000001 = -D_detA_100000000001 ;

       const double D_vU_011000000000 = -D_detA_011000000000 ;

       const double D_vU_010100000000 = -D_detA_010100000000 ;

       const double D_vU_010001000000 = -D_detA_010001000000 ;

       const double D_vU_010000100000 = -D_detA_010000100000 ;

       const double D_vU_010000001000 = -D_detA_010000001000 ;

       const double D_vU_010000000100 = -D_detA_010000000100 ;

       const double D_vU_010000000001 = -D_detA_010000000001 ;

       const double D_vU_001100000000 = -D_detA_001100000000 ;

       const double D_vU_001010000000 = -D_detA_001010000000 ;

       const double D_vU_001000100000 = -D_detA_001000100000 ;

       const double D_vU_001000010000 = -D_detA_001000010000 ;

       const double D_vU_001000000100 = -D_detA_001000000100 ;

       const double D_vU_001000000010 = -D_detA_001000000010 ;

       const double D_vU_000100010000 = -D_detA_000100010000 ;

       const double D_vU_000100001000 = -D_detA_000100001000 ;

       const double D_vU_000100000010 = -D_detA_000100000010 ;

       const double D_vU_000100000001 = -D_detA_000100000001 ;

       const double D_vU_000010100000 = -D_detA_000010100000 ;

       const double D_vU_000010001000 = -D_detA_000010001000 ;

       const double D_vU_000010000100 = -D_detA_000010000100 ;

       const double D_vU_000010000001 = -D_detA_000010000001 ;

       const double D_vU_000001100000 = -D_detA_000001100000 ;

       const double D_vU_000001010000 = -D_detA_000001010000 ;

       const double D_vU_000001000100 = -D_detA_000001000100 ;

       const double D_vU_000001000010 = -D_detA_000001000010 ;

       const double D_vU_000000100010 = -D_detA_000000100010 ;

       const double D_vU_000000100001 = -D_detA_000000100001 ;

       const double D_vU_000000010100 = -D_detA_000000010100 ;

       const double D_vU_000000010001 = -D_detA_000000010001 ;

       const double D_vU_000000001100 = -D_detA_000000001100 ;

       const double D_vU_000000001010 = -D_detA_000000001010 ;

       const double D_met_200000000000 = 2*hv*MyPow2(D_vU_100000000000) ;

       const double D_met_110000000000 = 2*D_vU_010000000000*D_vU_100000000000*hv + 2*D_vU_110000000000*hv*vU ;

       const double D_met_101000000000 = 2*D_vU_001000000000*D_vU_100000000000*hv + 2*D_vU_101000000000*hv*vU ;

       const double D_met_100100000000 = 2*D_vU_000100000000*D_vU_100000000000*hv ;

       const double D_met_100010000000 = 2*D_vU_000010000000*D_vU_100000000000*hv + 2*D_vU_100010000000*hv*vU ;

       const double D_met_100001000000 = 2*D_vU_000001000000*D_vU_100000000000*hv + 2*D_vU_100001000000*hv*vU ;

       const double D_met_100000100000 = 2*D_vU_000000100000*D_vU_100000000000*hv ;

       const double D_met_100000010000 = 2*D_vU_000000010000*D_vU_100000000000*hv + 2*D_vU_100000010000*hv*vU ;

       const double D_met_100000001000 = 2*D_vU_000000001000*D_vU_100000000000*hv + 2*D_vU_100000001000*hv*vU ;

       const double D_met_100000000100 = 2*D_vU_000000000100*D_vU_100000000000*hv ;

       const double D_met_100000000010 = 2*D_vU_000000000010*D_vU_100000000000*hv + 2*D_vU_100000000010*hv*vU ;

       const double D_met_100000000001 = 2*D_vU_000000000001*D_vU_100000000000*hv + 2*D_vU_100000000001*hv*vU ;

       const double D_met_020000000000 = 2*hv*MyPow2(D_vU_010000000000) ;

       const double D_met_011000000000 = 2*D_vU_001000000000*D_vU_010000000000*hv + 2*D_vU_011000000000*hv*vU ;

       const double D_met_010100000000 = 2*D_vU_000100000000*D_vU_010000000000*hv + 2*D_vU_010100000000*hv*vU ;

       const double D_met_010010000000 = 2*D_vU_000010000000*D_vU_010000000000*hv ;

       const double D_met_010001000000 = 2*D_vU_000001000000*D_vU_010000000000*hv + 2*D_vU_010001000000*hv*vU ;

       const double D_met_010000100000 = 2*D_vU_000000100000*D_vU_010000000000*hv + 2*D_vU_010000100000*hv*vU ;

       const double D_met_010000010000 = 2*D_vU_000000010000*D_vU_010000000000*hv ;

       const double D_met_010000001000 = 2*D_vU_000000001000*D_vU_010000000000*hv + 2*D_vU_010000001000*hv*vU ;

       const double D_met_010000000100 = 2*D_vU_000000000100*D_vU_010000000000*hv + 2*D_vU_010000000100*hv*vU ;

       const double D_met_010000000010 = 2*D_vU_000000000010*D_vU_010000000000*hv ;

       const double D_met_010000000001 = 2*D_vU_000000000001*D_vU_010000000000*hv + 2*D_vU_010000000001*hv*vU ;

       const double D_met_002000000000 = 2*hv*MyPow2(D_vU_001000000000) ;

       const double D_met_001100000000 = 2*D_vU_000100000000*D_vU_001000000000*hv + 2*D_vU_001100000000*hv*vU ;

       const double D_met_001010000000 = 2*D_vU_000010000000*D_vU_001000000000*hv + 2*D_vU_001010000000*hv*vU ;

       const double D_met_001001000000 = 2*D_vU_000001000000*D_vU_001000000000*hv ;

       const double D_met_001000100000 = 2*D_vU_000000100000*D_vU_001000000000*hv + 2*D_vU_001000100000*hv*vU ;

       const double D_met_001000010000 = 2*D_vU_000000010000*D_vU_001000000000*hv + 2*D_vU_001000010000*hv*vU ;

       const double D_met_001000001000 = 2*D_vU_000000001000*D_vU_001000000000*hv ;

       const double D_met_001000000100 = 2*D_vU_000000000100*D_vU_001000000000*hv + 2*D_vU_001000000100*hv*vU ;

       const double D_met_001000000010 = 2*D_vU_000000000010*D_vU_001000000000*hv + 2*D_vU_001000000010*hv*vU ;

       const double D_met_001000000001 = 2*D_vU_000000000001*D_vU_001000000000*hv ;

       const double D_met_000200000000 = 2*hv*MyPow2(D_vU_000100000000) ;

       const double D_met_000110000000 = 2*D_vU_000010000000*D_vU_000100000000*hv ;

       const double D_met_000101000000 = 2*D_vU_000001000000*D_vU_000100000000*hv ;

       const double D_met_000100100000 = 2*D_vU_000000100000*D_vU_000100000000*hv ;

       const double D_met_000100010000 = 2*D_vU_000000010000*D_vU_000100000000*hv + 2*D_vU_000100010000*hv*vU ;

       const double D_met_000100001000 = 2*D_vU_000000001000*D_vU_000100000000*hv + 2*D_vU_000100001000*hv*vU ;

       const double D_met_000100000100 = 2*D_vU_000000000100*D_vU_000100000000*hv ;

       const double D_met_000100000010 = 2*D_vU_000000000010*D_vU_000100000000*hv + 2*D_vU_000100000010*hv*vU ;

       const double D_met_000100000001 = 2*D_vU_000000000001*D_vU_000100000000*hv + 2*D_vU_000100000001*hv*vU ;

       const double D_met_000020000000 = 2*hv*MyPow2(D_vU_000010000000) ;

       const double D_met_000011000000 = 2*D_vU_000001000000*D_vU_000010000000*hv ;

       const double D_met_000010100000 = 2*D_vU_000000100000*D_vU_000010000000*hv + 2*D_vU_000010100000*hv*vU ;

       const double D_met_000010010000 = 2*D_vU_000000010000*D_vU_000010000000*hv ;

       const double D_met_000010001000 = 2*D_vU_000000001000*D_vU_000010000000*hv + 2*D_vU_000010001000*hv*vU ;

       const double D_met_000010000100 = 2*D_vU_000000000100*D_vU_000010000000*hv + 2*D_vU_000010000100*hv*vU ;

       const double D_met_000010000010 = 2*D_vU_000000000010*D_vU_000010000000*hv ;

       const double D_met_000010000001 = 2*D_vU_000000000001*D_vU_000010000000*hv + 2*D_vU_000010000001*hv*vU ;

       const double D_met_000002000000 = 2*hv*MyPow2(D_vU_000001000000) ;

       const double D_met_000001100000 = 2*D_vU_000000100000*D_vU_000001000000*hv + 2*D_vU_000001100000*hv*vU ;

       const double D_met_000001010000 = 2*D_vU_000000010000*D_vU_000001000000*hv + 2*D_vU_000001010000*hv*vU ;

       const double D_met_000001001000 = 2*D_vU_000000001000*D_vU_000001000000*hv ;

       const double D_met_000001000100 = 2*D_vU_000000000100*D_vU_000001000000*hv + 2*D_vU_000001000100*hv*vU ;

       const double D_met_000001000010 = 2*D_vU_000000000010*D_vU_000001000000*hv + 2*D_vU_000001000010*hv*vU ;

       const double D_met_000001000001 = 2*D_vU_000000000001*D_vU_000001000000*hv ;

       const double D_met_000000200000 = 2*hv*MyPow2(D_vU_000000100000) ;

       const double D_met_000000110000 = 2*D_vU_000000010000*D_vU_000000100000*hv ;

       const double D_met_000000101000 = 2*D_vU_000000001000*D_vU_000000100000*hv ;

       const double D_met_000000100100 = 2*D_vU_000000000100*D_vU_000000100000*hv ;

       const double D_met_000000100010 = 2*D_vU_000000000010*D_vU_000000100000*hv + 2*D_vU_000000100010*hv*vU ;

       const double D_met_000000100001 = 2*D_vU_000000000001*D_vU_000000100000*hv + 2*D_vU_000000100001*hv*vU ;

       const double D_met_000000020000 = 2*hv*MyPow2(D_vU_000000010000) ;

       const double D_met_000000011000 = 2*D_vU_000000001000*D_vU_000000010000*hv ;

       const double D_met_000000010100 = 2*D_vU_000000000100*D_vU_000000010000*hv + 2*D_vU_000000010100*hv*vU ;

       const double D_met_000000010010 = 2*D_vU_000000000010*D_vU_000000010000*hv ;

       const double D_met_000000010001 = 2*D_vU_000000000001*D_vU_000000010000*hv + 2*D_vU_000000010001*hv*vU ;

       const double D_met_000000002000 = 2*hv*MyPow2(D_vU_000000001000) ;

       const double D_met_000000001100 = 2*D_vU_000000000100*D_vU_000000001000*hv + 2*D_vU_000000001100*hv*vU ;

       const double D_met_000000001010 = 2*D_vU_000000000010*D_vU_000000001000*hv + 2*D_vU_000000001010*hv*vU ;

       const double D_met_000000001001 = 2*D_vU_000000000001*D_vU_000000001000*hv ;

       const double D_met_000000000200 = 2*hv*MyPow2(D_vU_000000000100) ;

       const double D_met_000000000110 = 2*D_vU_000000000010*D_vU_000000000100*hv ;

       const double D_met_000000000101 = 2*D_vU_000000000001*D_vU_000000000100*hv ;

       const double D_met_000000000020 = 2*hv*MyPow2(D_vU_000000000010) ;

       const double D_met_000000000011 = 2*D_vU_000000000001*D_vU_000000000010*hv ;

       const double D_met_000000000002 = 2*hv*MyPow2(D_vU_000000000001) ;
                  
       HESS(II(0),0,II(0),0) += D_met_200000000000 ;

       HESS(II(0),0,II(0),1) += D_met_110000000000 ;

       HESS(II(0),0,II(0),2) += D_met_101000000000 ;

       HESS(II(0),0,II(1),0) += D_met_100100000000 ;

       HESS(II(0),0,II(1),1) += D_met_100010000000 ;

       HESS(II(0),0,II(1),2) += D_met_100001000000 ;

       HESS(II(0),0,II(2),0) += D_met_100000100000 ;

       HESS(II(0),0,II(2),1) += D_met_100000010000 ;

       HESS(II(0),0,II(2),2) += D_met_100000001000 ;

       HESS(II(0),0,II(3),0) += D_met_100000000100 ;

       HESS(II(0),0,II(3),1) += D_met_100000000010 ;

       HESS(II(0),0,II(3),2) += D_met_100000000001 ;

       HESS(II(0),1,II(0),0) += D_met_110000000000 ;

       HESS(II(0),1,II(0),1) += D_met_020000000000 ;

       HESS(II(0),1,II(0),2) += D_met_011000000000 ;

       HESS(II(0),1,II(1),0) += D_met_010100000000 ;

       HESS(II(0),1,II(1),1) += D_met_010010000000 ;

       HESS(II(0),1,II(1),2) += D_met_010001000000 ;

       HESS(II(0),1,II(2),0) += D_met_010000100000 ;

       HESS(II(0),1,II(2),1) += D_met_010000010000 ;

       HESS(II(0),1,II(2),2) += D_met_010000001000 ;

       HESS(II(0),1,II(3),0) += D_met_010000000100 ;

       HESS(II(0),1,II(3),1) += D_met_010000000010 ;

       HESS(II(0),1,II(3),2) += D_met_010000000001 ;

       HESS(II(0),2,II(0),0) += D_met_101000000000 ;

       HESS(II(0),2,II(0),1) += D_met_011000000000 ;

       HESS(II(0),2,II(0),2) += D_met_002000000000 ;

       HESS(II(0),2,II(1),0) += D_met_001100000000 ;

       HESS(II(0),2,II(1),1) += D_met_001010000000 ;

       HESS(II(0),2,II(1),2) += D_met_001001000000 ;

       HESS(II(0),2,II(2),0) += D_met_001000100000 ;

       HESS(II(0),2,II(2),1) += D_met_001000010000 ;

       HESS(II(0),2,II(2),2) += D_met_001000001000 ;

       HESS(II(0),2,II(3),0) += D_met_001000000100 ;

       HESS(II(0),2,II(3),1) += D_met_001000000010 ;

       HESS(II(0),2,II(3),2) += D_met_001000000001 ;

       HESS(II(1),0,II(0),0) += D_met_100100000000 ;

       HESS(II(1),0,II(0),1) += D_met_010100000000 ;

       HESS(II(1),0,II(0),2) += D_met_001100000000 ;

       HESS(II(1),0,II(1),0) += D_met_000200000000 ;

       HESS(II(1),0,II(1),1) += D_met_000110000000 ;

       HESS(II(1),0,II(1),2) += D_met_000101000000 ;

       HESS(II(1),0,II(2),0) += D_met_000100100000 ;

       HESS(II(1),0,II(2),1) += D_met_000100010000 ;

       HESS(II(1),0,II(2),2) += D_met_000100001000 ;

       HESS(II(1),0,II(3),0) += D_met_000100000100 ;

       HESS(II(1),0,II(3),1) += D_met_000100000010 ;

       HESS(II(1),0,II(3),2) += D_met_000100000001 ;

       HESS(II(1),1,II(0),0) += D_met_100010000000 ;

       HESS(II(1),1,II(0),1) += D_met_010010000000 ;

       HESS(II(1),1,II(0),2) += D_met_001010000000 ;

       HESS(II(1),1,II(1),0) += D_met_000110000000 ;

       HESS(II(1),1,II(1),1) += D_met_000020000000 ;

       HESS(II(1),1,II(1),2) += D_met_000011000000 ;

       HESS(II(1),1,II(2),0) += D_met_000010100000 ;

       HESS(II(1),1,II(2),1) += D_met_000010010000 ;

       HESS(II(1),1,II(2),2) += D_met_000010001000 ;

       HESS(II(1),1,II(3),0) += D_met_000010000100 ;

       HESS(II(1),1,II(3),1) += D_met_000010000010 ;

       HESS(II(1),1,II(3),2) += D_met_000010000001 ;

       HESS(II(1),2,II(0),0) += D_met_100001000000 ;

       HESS(II(1),2,II(0),1) += D_met_010001000000 ;

       HESS(II(1),2,II(0),2) += D_met_001001000000 ;

       HESS(II(1),2,II(1),0) += D_met_000101000000 ;

       HESS(II(1),2,II(1),1) += D_met_000011000000 ;

       HESS(II(1),2,II(1),2) += D_met_000002000000 ;

       HESS(II(1),2,II(2),0) += D_met_000001100000 ;

       HESS(II(1),2,II(2),1) += D_met_000001010000 ;

       HESS(II(1),2,II(2),2) += D_met_000001001000 ;

       HESS(II(1),2,II(3),0) += D_met_000001000100 ;

       HESS(II(1),2,II(3),1) += D_met_000001000010 ;

       HESS(II(1),2,II(3),2) += D_met_000001000001 ;

       HESS(II(2),0,II(0),0) += D_met_100000100000 ;

       HESS(II(2),0,II(0),1) += D_met_010000100000 ;

       HESS(II(2),0,II(0),2) += D_met_001000100000 ;

       HESS(II(2),0,II(1),0) += D_met_000100100000 ;

       HESS(II(2),0,II(1),1) += D_met_000010100000 ;

       HESS(II(2),0,II(1),2) += D_met_000001100000 ;

       HESS(II(2),0,II(2),0) += D_met_000000200000 ;

       HESS(II(2),0,II(2),1) += D_met_000000110000 ;

       HESS(II(2),0,II(2),2) += D_met_000000101000 ;

       HESS(II(2),0,II(3),0) += D_met_000000100100 ;

       HESS(II(2),0,II(3),1) += D_met_000000100010 ;

       HESS(II(2),0,II(3),2) += D_met_000000100001 ;

       HESS(II(2),1,II(0),0) += D_met_100000010000 ;

       HESS(II(2),1,II(0),1) += D_met_010000010000 ;

       HESS(II(2),1,II(0),2) += D_met_001000010000 ;

       HESS(II(2),1,II(1),0) += D_met_000100010000 ;

       HESS(II(2),1,II(1),1) += D_met_000010010000 ;

       HESS(II(2),1,II(1),2) += D_met_000001010000 ;

       HESS(II(2),1,II(2),0) += D_met_000000110000 ;

       HESS(II(2),1,II(2),1) += D_met_000000020000 ;

       HESS(II(2),1,II(2),2) += D_met_000000011000 ;

       HESS(II(2),1,II(3),0) += D_met_000000010100 ;

       HESS(II(2),1,II(3),1) += D_met_000000010010 ;

       HESS(II(2),1,II(3),2) += D_met_000000010001 ;

       HESS(II(2),2,II(0),0) += D_met_100000001000 ;

       HESS(II(2),2,II(0),1) += D_met_010000001000 ;

       HESS(II(2),2,II(0),2) += D_met_001000001000 ;

       HESS(II(2),2,II(1),0) += D_met_000100001000 ;

       HESS(II(2),2,II(1),1) += D_met_000010001000 ;

       HESS(II(2),2,II(1),2) += D_met_000001001000 ;

       HESS(II(2),2,II(2),0) += D_met_000000101000 ;

       HESS(II(2),2,II(2),1) += D_met_000000011000 ;

       HESS(II(2),2,II(2),2) += D_met_000000002000 ;

       HESS(II(2),2,II(3),0) += D_met_000000001100 ;

       HESS(II(2),2,II(3),1) += D_met_000000001010 ;

       HESS(II(2),2,II(3),2) += D_met_000000001001 ;

       HESS(II(3),0,II(0),0) += D_met_100000000100 ;

       HESS(II(3),0,II(0),1) += D_met_010000000100 ;

       HESS(II(3),0,II(0),2) += D_met_001000000100 ;

       HESS(II(3),0,II(1),0) += D_met_000100000100 ;

       HESS(II(3),0,II(1),1) += D_met_000010000100 ;

       HESS(II(3),0,II(1),2) += D_met_000001000100 ;

       HESS(II(3),0,II(2),0) += D_met_000000100100 ;

       HESS(II(3),0,II(2),1) += D_met_000000010100 ;

       HESS(II(3),0,II(2),2) += D_met_000000001100 ;

       HESS(II(3),0,II(3),0) += D_met_000000000200 ;

       HESS(II(3),0,II(3),1) += D_met_000000000110 ;

       HESS(II(3),0,II(3),2) += D_met_000000000101 ;

       HESS(II(3),1,II(0),0) += D_met_100000000010 ;

       HESS(II(3),1,II(0),1) += D_met_010000000010 ;

       HESS(II(3),1,II(0),2) += D_met_001000000010 ;

       HESS(II(3),1,II(1),0) += D_met_000100000010 ;

       HESS(II(3),1,II(1),1) += D_met_000010000010 ;

       HESS(II(3),1,II(1),2) += D_met_000001000010 ;

       HESS(II(3),1,II(2),0) += D_met_000000100010 ;

       HESS(II(3),1,II(2),1) += D_met_000000010010 ;

       HESS(II(3),1,II(2),2) += D_met_000000001010 ;

       HESS(II(3),1,II(3),0) += D_met_000000000110 ;

       HESS(II(3),1,II(3),1) += D_met_000000000020 ;

       HESS(II(3),1,II(3),2) += D_met_000000000011 ;

       HESS(II(3),2,II(0),0) += D_met_100000000001 ;

       HESS(II(3),2,II(0),1) += D_met_010000000001 ;

       HESS(II(3),2,II(0),2) += D_met_001000000001 ;

       HESS(II(3),2,II(1),0) += D_met_000100000001 ;

       HESS(II(3),2,II(1),1) += D_met_000010000001 ;

       HESS(II(3),2,II(1),2) += D_met_000001000001 ;

       HESS(II(3),2,II(2),0) += D_met_000000100001 ;

       HESS(II(3),2,II(2),1) += D_met_000000010001 ;

       HESS(II(3),2,II(2),2) += D_met_000000001001 ;

       HESS(II(3),2,II(3),0) += D_met_000000000101 ;

       HESS(II(3),2,II(3),1) += D_met_000000000011 ;

       HESS(II(3),2,II(3),2) += D_met_000000000002 ;
                  vv = met; sdetA = detA; sdetW = detW;
                }
              else
                {
                  
       const double A00 = -X0(0) + X1(0) ;

       const double A01 = -X0(0) + X2(0) ;

       const double A02 = -X0(0) + X3(0) ;

       const double A10 = -X0(1) + X1(1) ;

       const double A11 = -X0(1) + X2(1) ;

       const double A12 = -X0(1) + X3(1) ;

       const double A20 = -X0(2) + X1(2) ;

       const double A21 = -X0(2) + X2(2) ;

       const double A22 = -X0(2) + X3(2) ;

       const double W00 = -WX0(0) + WX1(0) ;

       const double W01 = -WX0(0) + WX2(0) ;

       const double W02 = -WX0(0) + WX3(0) ;

       const double W10 = -WX0(1) + WX1(1) ;

       const double W11 = -WX0(1) + WX2(1) ;

       const double W12 = -WX0(1) + WX3(1) ;

       const double W20 = -WX0(2) + WX1(2) ;

       const double W21 = -WX0(2) + WX2(2) ;

       const double W22 = -WX0(2) + WX3(2) ;

       const double detA = A02*(-(A11*A20) + A10*A21) + A01*(A12*A20 - A10*A22) + A00*(-(A12*A21) + A11*A22) ;

       const double detW = -(W02*W11*W20) + W01*W12*W20 + W02*W10*W21 - W00*W12*W21 - W01*W10*W22 + W00*W11*W22 ;

       const double detWI = MyInverse(detW) ;

       const double WI00 = detWI*(-(W12*W21) + W11*W22) ;

       const double WI01 = detWI*(W02*W21 - W01*W22) ;

       const double WI02 = detWI*(-(W02*W11) + W01*W12) ;

       const double WI10 = detWI*(W12*W20 - W10*W22) ;

       const double WI11 = detWI*(-(W02*W20) + W00*W22) ;

       const double WI12 = detWI*(W02*W10 - W00*W12) ;

       const double WI20 = detWI*(-(W11*W20) + W10*W21) ;

       const double WI21 = detWI*(W01*W20 - W00*W21) ;

       const double WI22 = detWI*(-(W01*W10) + W00*W11) ;

       const double T00 = A00*WI00 + A01*WI10 + A02*WI20 ;

       const double T01 = A00*WI01 + A01*WI11 + A02*WI21 ;

       const double T02 = A00*WI02 + A01*WI12 + A02*WI22 ;

       const double T10 = A10*WI00 + A11*WI10 + A12*WI20 ;

       const double T11 = A10*WI01 + A11*WI11 + A12*WI21 ;

       const double T12 = A10*WI02 + A11*WI12 + A12*WI22 ;

       const double T20 = A20*WI00 + A21*WI10 + A22*WI20 ;

       const double T21 = A20*WI01 + A21*WI11 + A22*WI21 ;

       const double T22 = A20*WI02 + A21*WI12 + A22*WI22 ;

       const double vU = -detA + beta*detW ;

       const double hv = MyHeaviside(vU,0) ;

       const double met = hv*MyPow2(vU) ;

       const double D_A00_100000000000 = -1 ;

       const double D_A00_000100000000 = 1 ;

       const double D_A01_100000000000 = -1 ;

       const double D_A01_000000100000 = 1 ;

       const double D_A02_100000000000 = -1 ;

       const double D_A02_000000000100 = 1 ;

       const double D_A10_010000000000 = -1 ;

       const double D_A10_000010000000 = 1 ;

       const double D_A11_010000000000 = -1 ;

       const double D_A11_000000010000 = 1 ;

       const double D_A12_010000000000 = -1 ;

       const double D_A12_000000000010 = 1 ;

       const double D_A20_001000000000 = -1 ;

       const double D_A20_000001000000 = 1 ;

       const double D_A21_001000000000 = -1 ;

       const double D_A21_000000001000 = 1 ;

       const double D_A22_001000000000 = -1 ;

       const double D_A22_000000000001 = 1 ;

       const double D_detA_100000000000 = (-(A12*A21) + A11*A22)*D_A00_100000000000 + (A12*A20 -\
 
  A10*A22)*D_A01_100000000000 + (-(A11*A20) + A10*A21)*D_A02_100000000000 ;

       const double D_detA_010000000000 = A02*(A21*D_A10_010000000000 - A20*D_A11_010000000000) +\
 
  A01*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) + A00*(A22*D_A11_010000000000 - A21*D_A12_010000000000) ;

       const double D_detA_001000000000 = A02*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) +\
 
  A01*(A12*D_A20_001000000000 - A10*D_A22_001000000000) + A00*(-(A12*D_A21_001000000000) + A11*D_A22_001000000000) ;

       const double D_detA_000100000000 = (-(A12*A21) + A11*A22)*D_A00_000100000000 ;

       const double D_detA_000010000000 = A02*A21*D_A10_000010000000 - A01*A22*D_A10_000010000000 ;

       const double D_detA_000001000000 = -(A02*A11*D_A20_000001000000) + A01*A12*D_A20_000001000000 ;

       const double D_detA_000000100000 = (A12*A20 - A10*A22)*D_A01_000000100000 ;

       const double D_detA_000000010000 = -(A02*A20*D_A11_000000010000) + A00*A22*D_A11_000000010000 ;

       const double D_detA_000000001000 = A02*A10*D_A21_000000001000 - A00*A12*D_A21_000000001000 ;

       const double D_detA_000000000100 = (-(A11*A20) + A10*A21)*D_A02_000000000100 ;

       const double D_detA_000000000010 = A01*A20*D_A12_000000000010 - A00*A21*D_A12_000000000010 ;

       const double D_detA_000000000001 = -(A01*A10*D_A22_000000000001) + A00*A11*D_A22_000000000001 ;

       const double D_T00_100000000000 = D_A00_100000000000*WI00 + D_A01_100000000000*WI10 + D_A02_100000000000*WI20 ;

       const double D_T00_000100000000 = D_A00_000100000000*WI00 ;

       const double D_T00_000000100000 = D_A01_000000100000*WI10 ;

       const double D_T00_000000000100 = D_A02_000000000100*WI20 ;

       const double D_T01_100000000000 = D_A00_100000000000*WI01 + D_A01_100000000000*WI11 + D_A02_100000000000*WI21 ;

       const double D_T01_000100000000 = D_A00_000100000000*WI01 ;

       const double D_T01_000000100000 = D_A01_000000100000*WI11 ;

       const double D_T01_000000000100 = D_A02_000000000100*WI21 ;

       const double D_T02_100000000000 = D_A00_100000000000*WI02 + D_A01_100000000000*WI12 + D_A02_100000000000*WI22 ;

       const double D_T02_000100000000 = D_A00_000100000000*WI02 ;

       const double D_T02_000000100000 = D_A01_000000100000*WI12 ;

       const double D_T02_000000000100 = D_A02_000000000100*WI22 ;

       const double D_T10_010000000000 = D_A10_010000000000*WI00 + D_A11_010000000000*WI10 + D_A12_010000000000*WI20 ;

       const double D_T10_000010000000 = D_A10_000010000000*WI00 ;

       const double D_T10_000000010000 = D_A11_000000010000*WI10 ;

       const double D_T10_000000000010 = D_A12_000000000010*WI20 ;

       const double D_T11_010000000000 = D_A10_010000000000*WI01 + D_A11_010000000000*WI11 + D_A12_010000000000*WI21 ;

       const double D_T11_000010000000 = D_A10_000010000000*WI01 ;

       const double D_T11_000000010000 = D_A11_000000010000*WI11 ;

       const double D_T11_000000000010 = D_A12_000000000010*WI21 ;

       const double D_T12_010000000000 = D_A10_010000000000*WI02 + D_A11_010000000000*WI12 + D_A12_010000000000*WI22 ;

       const double D_T12_000010000000 = D_A10_000010000000*WI02 ;

       const double D_T12_000000010000 = D_A11_000000010000*WI12 ;

       const double D_T12_000000000010 = D_A12_000000000010*WI22 ;

       const double D_T20_001000000000 = D_A20_001000000000*WI00 + D_A21_001000000000*WI10 + D_A22_001000000000*WI20 ;

       const double D_T20_000001000000 = D_A20_000001000000*WI00 ;

       const double D_T20_000000001000 = D_A21_000000001000*WI10 ;

       const double D_T20_000000000001 = D_A22_000000000001*WI20 ;

       const double D_T21_001000000000 = D_A20_001000000000*WI01 + D_A21_001000000000*WI11 + D_A22_001000000000*WI21 ;

       const double D_T21_000001000000 = D_A20_000001000000*WI01 ;

       const double D_T21_000000001000 = D_A21_000000001000*WI11 ;

       const double D_T21_000000000001 = D_A22_000000000001*WI21 ;

       const double D_T22_001000000000 = D_A20_001000000000*WI02 + D_A21_001000000000*WI12 + D_A22_001000000000*WI22 ;

       const double D_T22_000001000000 = D_A20_000001000000*WI02 ;

       const double D_T22_000000001000 = D_A21_000000001000*WI12 ;

       const double D_T22_000000000001 = D_A22_000000000001*WI22 ;

       const double D_vU_100000000000 = -D_detA_100000000000 ;

       const double D_vU_010000000000 = -D_detA_010000000000 ;

       const double D_vU_001000000000 = -D_detA_001000000000 ;

       const double D_vU_000100000000 = -D_detA_000100000000 ;

       const double D_vU_000010000000 = -D_detA_000010000000 ;

       const double D_vU_000001000000 = -D_detA_000001000000 ;

       const double D_vU_000000100000 = -D_detA_000000100000 ;

       const double D_vU_000000010000 = -D_detA_000000010000 ;

       const double D_vU_000000001000 = -D_detA_000000001000 ;

       const double D_vU_000000000100 = -D_detA_000000000100 ;

       const double D_vU_000000000010 = -D_detA_000000000010 ;

       const double D_vU_000000000001 = -D_detA_000000000001 ;

       const double D_met_100000000000 = 2*D_vU_100000000000*hv*vU ;

       const double D_met_010000000000 = 2*D_vU_010000000000*hv*vU ;

       const double D_met_001000000000 = 2*D_vU_001000000000*hv*vU ;

       const double D_met_000100000000 = 2*D_vU_000100000000*hv*vU ;

       const double D_met_000010000000 = 2*D_vU_000010000000*hv*vU ;

       const double D_met_000001000000 = 2*D_vU_000001000000*hv*vU ;

       const double D_met_000000100000 = 2*D_vU_000000100000*hv*vU ;

       const double D_met_000000010000 = 2*D_vU_000000010000*hv*vU ;

       const double D_met_000000001000 = 2*D_vU_000000001000*hv*vU ;

       const double D_met_000000000100 = 2*D_vU_000000000100*hv*vU ;

       const double D_met_000000000010 = 2*D_vU_000000000010*hv*vU ;

       const double D_met_000000000001 = 2*D_vU_000000000001*hv*vU ;

       const double D_detA_110000000000 = D_A02_100000000000*(A21*D_A10_010000000000 - A20*D_A11_010000000000) +\
 
  D_A01_100000000000*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) + D_A00_100000000000*(A22*D_A11_010000000000\
 
  - A21*D_A12_010000000000) ;

       const double D_detA_101000000000 = D_A02_100000000000*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) +\
 
  D_A01_100000000000*(A12*D_A20_001000000000 - A10*D_A22_001000000000) + D_A00_100000000000*(-(A12*D_A21_001000000000)\
 
  + A11*D_A22_001000000000) ;

       const double D_detA_100010000000 = -(A22*D_A01_100000000000*D_A10_000010000000) +\
 
  A21*D_A02_100000000000*D_A10_000010000000 ;

       const double D_detA_100001000000 = A12*D_A01_100000000000*D_A20_000001000000 -\
 
  A11*D_A02_100000000000*D_A20_000001000000 ;

       const double D_detA_100000010000 = A22*D_A00_100000000000*D_A11_000000010000 -\
 
  A20*D_A02_100000000000*D_A11_000000010000 ;

       const double D_detA_100000001000 = -(A12*D_A00_100000000000*D_A21_000000001000) +\
 
  A10*D_A02_100000000000*D_A21_000000001000 ;

       const double D_detA_100000000010 = -(A21*D_A00_100000000000*D_A12_000000000010) +\
 
  A20*D_A01_100000000000*D_A12_000000000010 ;

       const double D_detA_100000000001 = A11*D_A00_100000000000*D_A22_000000000001 -\
 
  A10*D_A01_100000000000*D_A22_000000000001 ;

       const double D_detA_011000000000 = A02*(-(D_A11_010000000000*D_A20_001000000000) +\
 
  D_A10_010000000000*D_A21_001000000000) + A01*(D_A12_010000000000*D_A20_001000000000 -\
 
  D_A10_010000000000*D_A22_001000000000) + A00*(-(D_A12_010000000000*D_A21_001000000000) +\
 
  D_A11_010000000000*D_A22_001000000000) ;

       const double D_detA_010100000000 = D_A00_000100000000*(A22*D_A11_010000000000 - A21*D_A12_010000000000) ;

       const double D_detA_010001000000 = -(A02*D_A11_010000000000*D_A20_000001000000) +\
 
  A01*D_A12_010000000000*D_A20_000001000000 ;

       const double D_detA_010000100000 = D_A01_000000100000*(-(A22*D_A10_010000000000) + A20*D_A12_010000000000) ;

       const double D_detA_010000001000 = A02*D_A10_010000000000*D_A21_000000001000 -\
 
  A00*D_A12_010000000000*D_A21_000000001000 ;

       const double D_detA_010000000100 = D_A02_000000000100*(A21*D_A10_010000000000 - A20*D_A11_010000000000) ;

       const double D_detA_010000000001 = -(A01*D_A10_010000000000*D_A22_000000000001) +\
 
  A00*D_A11_010000000000*D_A22_000000000001 ;

       const double D_detA_001100000000 = D_A00_000100000000*(-(A12*D_A21_001000000000) + A11*D_A22_001000000000) ;

       const double D_detA_001010000000 = A02*D_A10_000010000000*D_A21_001000000000 -\
 
  A01*D_A10_000010000000*D_A22_001000000000 ;

       const double D_detA_001000100000 = D_A01_000000100000*(A12*D_A20_001000000000 - A10*D_A22_001000000000) ;

       const double D_detA_001000010000 = -(A02*D_A11_000000010000*D_A20_001000000000) +\
 
  A00*D_A11_000000010000*D_A22_001000000000 ;

       const double D_detA_001000000100 = D_A02_000000000100*(-(A11*D_A20_001000000000) + A10*D_A21_001000000000) ;

       const double D_detA_001000000010 = A01*D_A12_000000000010*D_A20_001000000000 -\
 
  A00*D_A12_000000000010*D_A21_001000000000 ;

       const double D_detA_000100010000 = A22*D_A00_000100000000*D_A11_000000010000 ;

       const double D_detA_000100001000 = -(A12*D_A00_000100000000*D_A21_000000001000) ;

       const double D_detA_000100000010 = -(A21*D_A00_000100000000*D_A12_000000000010) ;

       const double D_detA_000100000001 = A11*D_A00_000100000000*D_A22_000000000001 ;

       const double D_detA_000010100000 = -(A22*D_A01_000000100000*D_A10_000010000000) ;

       const double D_detA_000010001000 = A02*D_A10_000010000000*D_A21_000000001000 ;

       const double D_detA_000010000100 = A21*D_A02_000000000100*D_A10_000010000000 ;

       const double D_detA_000010000001 = -(A01*D_A10_000010000000*D_A22_000000000001) ;

       const double D_detA_000001100000 = A12*D_A01_000000100000*D_A20_000001000000 ;

       const double D_detA_000001010000 = -(A02*D_A11_000000010000*D_A20_000001000000) ;

       const double D_detA_000001000100 = -(A11*D_A02_000000000100*D_A20_000001000000) ;

       const double D_detA_000001000010 = A01*D_A12_000000000010*D_A20_000001000000 ;

       const double D_detA_000000100010 = A20*D_A01_000000100000*D_A12_000000000010 ;

       const double D_detA_000000100001 = -(A10*D_A01_000000100000*D_A22_000000000001) ;

       const double D_detA_000000010100 = -(A20*D_A02_000000000100*D_A11_000000010000) ;

       const double D_detA_000000010001 = A00*D_A11_000000010000*D_A22_000000000001 ;

       const double D_detA_000000001100 = A10*D_A02_000000000100*D_A21_000000001000 ;

       const double D_detA_000000001010 = -(A00*D_A12_000000000010*D_A21_000000001000) ;

       const double D_vU_110000000000 = -D_detA_110000000000 ;

       const double D_vU_101000000000 = -D_detA_101000000000 ;

       const double D_vU_100010000000 = -D_detA_100010000000 ;

       const double D_vU_100001000000 = -D_detA_100001000000 ;

       const double D_vU_100000010000 = -D_detA_100000010000 ;

       const double D_vU_100000001000 = -D_detA_100000001000 ;

       const double D_vU_100000000010 = -D_detA_100000000010 ;

       const double D_vU_100000000001 = -D_detA_100000000001 ;

       const double D_vU_011000000000 = -D_detA_011000000000 ;

       const double D_vU_010100000000 = -D_detA_010100000000 ;

       const double D_vU_010001000000 = -D_detA_010001000000 ;

       const double D_vU_010000100000 = -D_detA_010000100000 ;

       const double D_vU_010000001000 = -D_detA_010000001000 ;

       const double D_vU_010000000100 = -D_detA_010000000100 ;

       const double D_vU_010000000001 = -D_detA_010000000001 ;

       const double D_vU_001100000000 = -D_detA_001100000000 ;

       const double D_vU_001010000000 = -D_detA_001010000000 ;

       const double D_vU_001000100000 = -D_detA_001000100000 ;

       const double D_vU_001000010000 = -D_detA_001000010000 ;

       const double D_vU_001000000100 = -D_detA_001000000100 ;

       const double D_vU_001000000010 = -D_detA_001000000010 ;

       const double D_vU_000100010000 = -D_detA_000100010000 ;

       const double D_vU_000100001000 = -D_detA_000100001000 ;

       const double D_vU_000100000010 = -D_detA_000100000010 ;

       const double D_vU_000100000001 = -D_detA_000100000001 ;

       const double D_vU_000010100000 = -D_detA_000010100000 ;

       const double D_vU_000010001000 = -D_detA_000010001000 ;

       const double D_vU_000010000100 = -D_detA_000010000100 ;

       const double D_vU_000010000001 = -D_detA_000010000001 ;

       const double D_vU_000001100000 = -D_detA_000001100000 ;

       const double D_vU_000001010000 = -D_detA_000001010000 ;

       const double D_vU_000001000100 = -D_detA_000001000100 ;

       const double D_vU_000001000010 = -D_detA_000001000010 ;

       const double D_vU_000000100010 = -D_detA_000000100010 ;

       const double D_vU_000000100001 = -D_detA_000000100001 ;

       const double D_vU_000000010100 = -D_detA_000000010100 ;

       const double D_vU_000000010001 = -D_detA_000000010001 ;

       const double D_vU_000000001100 = -D_detA_000000001100 ;

       const double D_vU_000000001010 = -D_detA_000000001010 ;

       const double D_met_200000000000 = 2*hv*MyPow2(D_vU_100000000000) ;

       const double D_met_110000000000 = 2*D_vU_010000000000*D_vU_100000000000*hv + 2*D_vU_110000000000*hv*vU ;

       const double D_met_101000000000 = 2*D_vU_001000000000*D_vU_100000000000*hv + 2*D_vU_101000000000*hv*vU ;

       const double D_met_100100000000 = 2*D_vU_000100000000*D_vU_100000000000*hv ;

       const double D_met_100010000000 = 2*D_vU_000010000000*D_vU_100000000000*hv + 2*D_vU_100010000000*hv*vU ;

       const double D_met_100001000000 = 2*D_vU_000001000000*D_vU_100000000000*hv + 2*D_vU_100001000000*hv*vU ;

       const double D_met_100000100000 = 2*D_vU_000000100000*D_vU_100000000000*hv ;

       const double D_met_100000010000 = 2*D_vU_000000010000*D_vU_100000000000*hv + 2*D_vU_100000010000*hv*vU ;

       const double D_met_100000001000 = 2*D_vU_000000001000*D_vU_100000000000*hv + 2*D_vU_100000001000*hv*vU ;

       const double D_met_100000000100 = 2*D_vU_000000000100*D_vU_100000000000*hv ;

       const double D_met_100000000010 = 2*D_vU_000000000010*D_vU_100000000000*hv + 2*D_vU_100000000010*hv*vU ;

       const double D_met_100000000001 = 2*D_vU_000000000001*D_vU_100000000000*hv + 2*D_vU_100000000001*hv*vU ;

       const double D_met_020000000000 = 2*hv*MyPow2(D_vU_010000000000) ;

       const double D_met_011000000000 = 2*D_vU_001000000000*D_vU_010000000000*hv + 2*D_vU_011000000000*hv*vU ;

       const double D_met_010100000000 = 2*D_vU_000100000000*D_vU_010000000000*hv + 2*D_vU_010100000000*hv*vU ;

       const double D_met_010010000000 = 2*D_vU_000010000000*D_vU_010000000000*hv ;

       const double D_met_010001000000 = 2*D_vU_000001000000*D_vU_010000000000*hv + 2*D_vU_010001000000*hv*vU ;

       const double D_met_010000100000 = 2*D_vU_000000100000*D_vU_010000000000*hv + 2*D_vU_010000100000*hv*vU ;

       const double D_met_010000010000 = 2*D_vU_000000010000*D_vU_010000000000*hv ;

       const double D_met_010000001000 = 2*D_vU_000000001000*D_vU_010000000000*hv + 2*D_vU_010000001000*hv*vU ;

       const double D_met_010000000100 = 2*D_vU_000000000100*D_vU_010000000000*hv + 2*D_vU_010000000100*hv*vU ;

       const double D_met_010000000010 = 2*D_vU_000000000010*D_vU_010000000000*hv ;

       const double D_met_010000000001 = 2*D_vU_000000000001*D_vU_010000000000*hv + 2*D_vU_010000000001*hv*vU ;

       const double D_met_002000000000 = 2*hv*MyPow2(D_vU_001000000000) ;

       const double D_met_001100000000 = 2*D_vU_000100000000*D_vU_001000000000*hv + 2*D_vU_001100000000*hv*vU ;

       const double D_met_001010000000 = 2*D_vU_000010000000*D_vU_001000000000*hv + 2*D_vU_001010000000*hv*vU ;

       const double D_met_001001000000 = 2*D_vU_000001000000*D_vU_001000000000*hv ;

       const double D_met_001000100000 = 2*D_vU_000000100000*D_vU_001000000000*hv + 2*D_vU_001000100000*hv*vU ;

       const double D_met_001000010000 = 2*D_vU_000000010000*D_vU_001000000000*hv + 2*D_vU_001000010000*hv*vU ;

       const double D_met_001000001000 = 2*D_vU_000000001000*D_vU_001000000000*hv ;

       const double D_met_001000000100 = 2*D_vU_000000000100*D_vU_001000000000*hv + 2*D_vU_001000000100*hv*vU ;

       const double D_met_001000000010 = 2*D_vU_000000000010*D_vU_001000000000*hv + 2*D_vU_001000000010*hv*vU ;

       const double D_met_001000000001 = 2*D_vU_000000000001*D_vU_001000000000*hv ;

       const double D_met_000200000000 = 2*hv*MyPow2(D_vU_000100000000) ;

       const double D_met_000110000000 = 2*D_vU_000010000000*D_vU_000100000000*hv ;

       const double D_met_000101000000 = 2*D_vU_000001000000*D_vU_000100000000*hv ;

       const double D_met_000100100000 = 2*D_vU_000000100000*D_vU_000100000000*hv ;

       const double D_met_000100010000 = 2*D_vU_000000010000*D_vU_000100000000*hv + 2*D_vU_000100010000*hv*vU ;

       const double D_met_000100001000 = 2*D_vU_000000001000*D_vU_000100000000*hv + 2*D_vU_000100001000*hv*vU ;

       const double D_met_000100000100 = 2*D_vU_000000000100*D_vU_000100000000*hv ;

       const double D_met_000100000010 = 2*D_vU_000000000010*D_vU_000100000000*hv + 2*D_vU_000100000010*hv*vU ;

       const double D_met_000100000001 = 2*D_vU_000000000001*D_vU_000100000000*hv + 2*D_vU_000100000001*hv*vU ;

       const double D_met_000020000000 = 2*hv*MyPow2(D_vU_000010000000) ;

       const double D_met_000011000000 = 2*D_vU_000001000000*D_vU_000010000000*hv ;

       const double D_met_000010100000 = 2*D_vU_000000100000*D_vU_000010000000*hv + 2*D_vU_000010100000*hv*vU ;

       const double D_met_000010010000 = 2*D_vU_000000010000*D_vU_000010000000*hv ;

       const double D_met_000010001000 = 2*D_vU_000000001000*D_vU_000010000000*hv + 2*D_vU_000010001000*hv*vU ;

       const double D_met_000010000100 = 2*D_vU_000000000100*D_vU_000010000000*hv + 2*D_vU_000010000100*hv*vU ;

       const double D_met_000010000010 = 2*D_vU_000000000010*D_vU_000010000000*hv ;

       const double D_met_000010000001 = 2*D_vU_000000000001*D_vU_000010000000*hv + 2*D_vU_000010000001*hv*vU ;

       const double D_met_000002000000 = 2*hv*MyPow2(D_vU_000001000000) ;

       const double D_met_000001100000 = 2*D_vU_000000100000*D_vU_000001000000*hv + 2*D_vU_000001100000*hv*vU ;

       const double D_met_000001010000 = 2*D_vU_000000010000*D_vU_000001000000*hv + 2*D_vU_000001010000*hv*vU ;

       const double D_met_000001001000 = 2*D_vU_000000001000*D_vU_000001000000*hv ;

       const double D_met_000001000100 = 2*D_vU_000000000100*D_vU_000001000000*hv + 2*D_vU_000001000100*hv*vU ;

       const double D_met_000001000010 = 2*D_vU_000000000010*D_vU_000001000000*hv + 2*D_vU_000001000010*hv*vU ;

       const double D_met_000001000001 = 2*D_vU_000000000001*D_vU_000001000000*hv ;

       const double D_met_000000200000 = 2*hv*MyPow2(D_vU_000000100000) ;

       const double D_met_000000110000 = 2*D_vU_000000010000*D_vU_000000100000*hv ;

       const double D_met_000000101000 = 2*D_vU_000000001000*D_vU_000000100000*hv ;

       const double D_met_000000100100 = 2*D_vU_000000000100*D_vU_000000100000*hv ;

       const double D_met_000000100010 = 2*D_vU_000000000010*D_vU_000000100000*hv + 2*D_vU_000000100010*hv*vU ;

       const double D_met_000000100001 = 2*D_vU_000000000001*D_vU_000000100000*hv + 2*D_vU_000000100001*hv*vU ;

       const double D_met_000000020000 = 2*hv*MyPow2(D_vU_000000010000) ;

       const double D_met_000000011000 = 2*D_vU_000000001000*D_vU_000000010000*hv ;

       const double D_met_000000010100 = 2*D_vU_000000000100*D_vU_000000010000*hv + 2*D_vU_000000010100*hv*vU ;

       const double D_met_000000010010 = 2*D_vU_000000000010*D_vU_000000010000*hv ;

       const double D_met_000000010001 = 2*D_vU_000000000001*D_vU_000000010000*hv + 2*D_vU_000000010001*hv*vU ;

       const double D_met_000000002000 = 2*hv*MyPow2(D_vU_000000001000) ;

       const double D_met_000000001100 = 2*D_vU_000000000100*D_vU_000000001000*hv + 2*D_vU_000000001100*hv*vU ;

       const double D_met_000000001010 = 2*D_vU_000000000010*D_vU_000000001000*hv + 2*D_vU_000000001010*hv*vU ;

       const double D_met_000000001001 = 2*D_vU_000000000001*D_vU_000000001000*hv ;

       const double D_met_000000000200 = 2*hv*MyPow2(D_vU_000000000100) ;

       const double D_met_000000000110 = 2*D_vU_000000000010*D_vU_000000000100*hv ;

       const double D_met_000000000101 = 2*D_vU_000000000001*D_vU_000000000100*hv ;

       const double D_met_000000000020 = 2*hv*MyPow2(D_vU_000000000010) ;

       const double D_met_000000000011 = 2*D_vU_000000000001*D_vU_000000000010*hv ;

       const double D_met_000000000002 = 2*hv*MyPow2(D_vU_000000000001) ;
                  
       GRAD(II(0),0) += D_met_100000000000 ;

       GRAD(II(0),1) += D_met_010000000000 ;

       GRAD(II(0),2) += D_met_001000000000 ;

       GRAD(II(1),0) += D_met_000100000000 ;

       GRAD(II(1),1) += D_met_000010000000 ;

       GRAD(II(1),2) += D_met_000001000000 ;

       GRAD(II(2),0) += D_met_000000100000 ;

       GRAD(II(2),1) += D_met_000000010000 ;

       GRAD(II(2),2) += D_met_000000001000 ;

       GRAD(II(3),0) += D_met_000000000100 ;

       GRAD(II(3),1) += D_met_000000000010 ;

       GRAD(II(3),2) += D_met_000000000001 ;
                  
       HESS(II(0),0,II(0),0) += D_met_200000000000 ;

       HESS(II(0),0,II(0),1) += D_met_110000000000 ;

       HESS(II(0),0,II(0),2) += D_met_101000000000 ;

       HESS(II(0),0,II(1),0) += D_met_100100000000 ;

       HESS(II(0),0,II(1),1) += D_met_100010000000 ;

       HESS(II(0),0,II(1),2) += D_met_100001000000 ;

       HESS(II(0),0,II(2),0) += D_met_100000100000 ;

       HESS(II(0),0,II(2),1) += D_met_100000010000 ;

       HESS(II(0),0,II(2),2) += D_met_100000001000 ;

       HESS(II(0),0,II(3),0) += D_met_100000000100 ;

       HESS(II(0),0,II(3),1) += D_met_100000000010 ;

       HESS(II(0),0,II(3),2) += D_met_100000000001 ;

       HESS(II(0),1,II(0),0) += D_met_110000000000 ;

       HESS(II(0),1,II(0),1) += D_met_020000000000 ;

       HESS(II(0),1,II(0),2) += D_met_011000000000 ;

       HESS(II(0),1,II(1),0) += D_met_010100000000 ;

       HESS(II(0),1,II(1),1) += D_met_010010000000 ;

       HESS(II(0),1,II(1),2) += D_met_010001000000 ;

       HESS(II(0),1,II(2),0) += D_met_010000100000 ;

       HESS(II(0),1,II(2),1) += D_met_010000010000 ;

       HESS(II(0),1,II(2),2) += D_met_010000001000 ;

       HESS(II(0),1,II(3),0) += D_met_010000000100 ;

       HESS(II(0),1,II(3),1) += D_met_010000000010 ;

       HESS(II(0),1,II(3),2) += D_met_010000000001 ;

       HESS(II(0),2,II(0),0) += D_met_101000000000 ;

       HESS(II(0),2,II(0),1) += D_met_011000000000 ;

       HESS(II(0),2,II(0),2) += D_met_002000000000 ;

       HESS(II(0),2,II(1),0) += D_met_001100000000 ;

       HESS(II(0),2,II(1),1) += D_met_001010000000 ;

       HESS(II(0),2,II(1),2) += D_met_001001000000 ;

       HESS(II(0),2,II(2),0) += D_met_001000100000 ;

       HESS(II(0),2,II(2),1) += D_met_001000010000 ;

       HESS(II(0),2,II(2),2) += D_met_001000001000 ;

       HESS(II(0),2,II(3),0) += D_met_001000000100 ;

       HESS(II(0),2,II(3),1) += D_met_001000000010 ;

       HESS(II(0),2,II(3),2) += D_met_001000000001 ;

       HESS(II(1),0,II(0),0) += D_met_100100000000 ;

       HESS(II(1),0,II(0),1) += D_met_010100000000 ;

       HESS(II(1),0,II(0),2) += D_met_001100000000 ;

       HESS(II(1),0,II(1),0) += D_met_000200000000 ;

       HESS(II(1),0,II(1),1) += D_met_000110000000 ;

       HESS(II(1),0,II(1),2) += D_met_000101000000 ;

       HESS(II(1),0,II(2),0) += D_met_000100100000 ;

       HESS(II(1),0,II(2),1) += D_met_000100010000 ;

       HESS(II(1),0,II(2),2) += D_met_000100001000 ;

       HESS(II(1),0,II(3),0) += D_met_000100000100 ;

       HESS(II(1),0,II(3),1) += D_met_000100000010 ;

       HESS(II(1),0,II(3),2) += D_met_000100000001 ;

       HESS(II(1),1,II(0),0) += D_met_100010000000 ;

       HESS(II(1),1,II(0),1) += D_met_010010000000 ;

       HESS(II(1),1,II(0),2) += D_met_001010000000 ;

       HESS(II(1),1,II(1),0) += D_met_000110000000 ;

       HESS(II(1),1,II(1),1) += D_met_000020000000 ;

       HESS(II(1),1,II(1),2) += D_met_000011000000 ;

       HESS(II(1),1,II(2),0) += D_met_000010100000 ;

       HESS(II(1),1,II(2),1) += D_met_000010010000 ;

       HESS(II(1),1,II(2),2) += D_met_000010001000 ;

       HESS(II(1),1,II(3),0) += D_met_000010000100 ;

       HESS(II(1),1,II(3),1) += D_met_000010000010 ;

       HESS(II(1),1,II(3),2) += D_met_000010000001 ;

       HESS(II(1),2,II(0),0) += D_met_100001000000 ;

       HESS(II(1),2,II(0),1) += D_met_010001000000 ;

       HESS(II(1),2,II(0),2) += D_met_001001000000 ;

       HESS(II(1),2,II(1),0) += D_met_000101000000 ;

       HESS(II(1),2,II(1),1) += D_met_000011000000 ;

       HESS(II(1),2,II(1),2) += D_met_000002000000 ;

       HESS(II(1),2,II(2),0) += D_met_000001100000 ;

       HESS(II(1),2,II(2),1) += D_met_000001010000 ;

       HESS(II(1),2,II(2),2) += D_met_000001001000 ;

       HESS(II(1),2,II(3),0) += D_met_000001000100 ;

       HESS(II(1),2,II(3),1) += D_met_000001000010 ;

       HESS(II(1),2,II(3),2) += D_met_000001000001 ;

       HESS(II(2),0,II(0),0) += D_met_100000100000 ;

       HESS(II(2),0,II(0),1) += D_met_010000100000 ;

       HESS(II(2),0,II(0),2) += D_met_001000100000 ;

       HESS(II(2),0,II(1),0) += D_met_000100100000 ;

       HESS(II(2),0,II(1),1) += D_met_000010100000 ;

       HESS(II(2),0,II(1),2) += D_met_000001100000 ;

       HESS(II(2),0,II(2),0) += D_met_000000200000 ;

       HESS(II(2),0,II(2),1) += D_met_000000110000 ;

       HESS(II(2),0,II(2),2) += D_met_000000101000 ;

       HESS(II(2),0,II(3),0) += D_met_000000100100 ;

       HESS(II(2),0,II(3),1) += D_met_000000100010 ;

       HESS(II(2),0,II(3),2) += D_met_000000100001 ;

       HESS(II(2),1,II(0),0) += D_met_100000010000 ;

       HESS(II(2),1,II(0),1) += D_met_010000010000 ;

       HESS(II(2),1,II(0),2) += D_met_001000010000 ;

       HESS(II(2),1,II(1),0) += D_met_000100010000 ;

       HESS(II(2),1,II(1),1) += D_met_000010010000 ;

       HESS(II(2),1,II(1),2) += D_met_000001010000 ;

       HESS(II(2),1,II(2),0) += D_met_000000110000 ;

       HESS(II(2),1,II(2),1) += D_met_000000020000 ;

       HESS(II(2),1,II(2),2) += D_met_000000011000 ;

       HESS(II(2),1,II(3),0) += D_met_000000010100 ;

       HESS(II(2),1,II(3),1) += D_met_000000010010 ;

       HESS(II(2),1,II(3),2) += D_met_000000010001 ;

       HESS(II(2),2,II(0),0) += D_met_100000001000 ;

       HESS(II(2),2,II(0),1) += D_met_010000001000 ;

       HESS(II(2),2,II(0),2) += D_met_001000001000 ;

       HESS(II(2),2,II(1),0) += D_met_000100001000 ;

       HESS(II(2),2,II(1),1) += D_met_000010001000 ;

       HESS(II(2),2,II(1),2) += D_met_000001001000 ;

       HESS(II(2),2,II(2),0) += D_met_000000101000 ;

       HESS(II(2),2,II(2),1) += D_met_000000011000 ;

       HESS(II(2),2,II(2),2) += D_met_000000002000 ;

       HESS(II(2),2,II(3),0) += D_met_000000001100 ;

       HESS(II(2),2,II(3),1) += D_met_000000001010 ;

       HESS(II(2),2,II(3),2) += D_met_000000001001 ;

       HESS(II(3),0,II(0),0) += D_met_100000000100 ;

       HESS(II(3),0,II(0),1) += D_met_010000000100 ;

       HESS(II(3),0,II(0),2) += D_met_001000000100 ;

       HESS(II(3),0,II(1),0) += D_met_000100000100 ;

       HESS(II(3),0,II(1),1) += D_met_000010000100 ;

       HESS(II(3),0,II(1),2) += D_met_000001000100 ;

       HESS(II(3),0,II(2),0) += D_met_000000100100 ;

       HESS(II(3),0,II(2),1) += D_met_000000010100 ;

       HESS(II(3),0,II(2),2) += D_met_000000001100 ;

       HESS(II(3),0,II(3),0) += D_met_000000000200 ;

       HESS(II(3),0,II(3),1) += D_met_000000000110 ;

       HESS(II(3),0,II(3),2) += D_met_000000000101 ;

       HESS(II(3),1,II(0),0) += D_met_100000000010 ;

       HESS(II(3),1,II(0),1) += D_met_010000000010 ;

       HESS(II(3),1,II(0),2) += D_met_001000000010 ;

       HESS(II(3),1,II(1),0) += D_met_000100000010 ;

       HESS(II(3),1,II(1),1) += D_met_000010000010 ;

       HESS(II(3),1,II(1),2) += D_met_000001000010 ;

       HESS(II(3),1,II(2),0) += D_met_000000100010 ;

       HESS(II(3),1,II(2),1) += D_met_000000010010 ;

       HESS(II(3),1,II(2),2) += D_met_000000001010 ;

       HESS(II(3),1,II(3),0) += D_met_000000000110 ;

       HESS(II(3),1,II(3),1) += D_met_000000000020 ;

       HESS(II(3),1,II(3),2) += D_met_000000000011 ;

       HESS(II(3),2,II(0),0) += D_met_100000000001 ;

       HESS(II(3),2,II(0),1) += D_met_010000000001 ;

       HESS(II(3),2,II(0),2) += D_met_001000000001 ;

       HESS(II(3),2,II(1),0) += D_met_000100000001 ;

       HESS(II(3),2,II(1),1) += D_met_000010000001 ;

       HESS(II(3),2,II(1),2) += D_met_000001000001 ;

       HESS(II(3),2,II(2),0) += D_met_000000100001 ;

       HESS(II(3),2,II(2),1) += D_met_000000010001 ;

       HESS(II(3),2,II(2),2) += D_met_000000001001 ;

       HESS(II(3),2,II(3),0) += D_met_000000000101 ;

       HESS(II(3),2,II(3),1) += D_met_000000000011 ;

       HESS(II(3),2,II(3),2) += D_met_000000000002 ;
                  vv = met; sdetA = detA; sdetW = detW;
                }
            }

          VERIFY_OP_ON(sdetW, >, 0.0, "bad reference mesh");
          if (sdetA <= 0.)
            {
              valid = false;
            }
          val_metric += vv;

        }
      val = val_metric;
      return val;
    }



  };

}

#undef X0
#undef X1
#undef X2
#undef X3
#undef WX0
#undef WX1
#undef WX2
#undef WX3

#undef GRAD
#undef HESS
#undef II
#undef normal

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic pop
#endif // GCC_VERSION
#endif // __GNUC__

#endif
#endif
