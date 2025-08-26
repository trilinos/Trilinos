// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>

#if !defined(NO_GEOM_SUPPORT)

#include "JacobianUtil.hpp"

#include "mpi.h"

  namespace percept {

    // Prism and Hex element descriptions
    static const int locs_prism[6][4] = {{0, 1, 2, 3}, {1, 2, 0, 4},
                                         {2, 0, 1, 5}, {3, 5, 4, 0},
                                         {4, 3, 5, 1}, {5, 4, 3, 2}};
    static const int locs_hex[8][4] = {{0, 1, 3, 4}, {1, 2, 0, 5},
                                       {2, 3, 1, 6}, {3, 0, 2, 7},
                                       {4, 7, 5, 0}, {5, 4, 6, 1},
                                       {6, 5, 7, 2}, {7, 6, 4, 3}};

    template<>    bool  JacobianUtilImpl<STKMesh>::jacobian_matrix_1D_in_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[2])
    {
      A.zero();
      A(0,0) = std::sqrt(SQR(x[1][0] - x[0][0]) + SQR(x[1][1] - x[0][1]) );
      A(1,1) = 1.0;
      A(2,2) = 1.0;
      detJ = A(0,0);
      return detJ == 0.0;
    }

    template<>    bool  JacobianUtilImpl<STKMesh>::jacobian_matrix_1D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[2])
    {
      A.zero();
      A(0,0) = std::sqrt(SQR(x[1][0] - x[0][0]) + SQR(x[1][1] - x[0][1]) + SQR(x[1][2] - x[0][2]) );
      A(1,1) = 1.0;
      A(2,2) = 1.0;
      detJ = A(0,0);
      return detJ == 0.0;
    }

    template<>    bool  JacobianUtilImpl<STKMesh>::jacobian_matrix_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
    {
      /* Calculate A */
      // x_xi, x_eta, x_zeta => A(ixyz, ixietazeta) = dx_i/dxi_j
      A(0,0) = (x[1][0] - x[0][0]);
      A(0,1) = (x[2][0] - x[0][0]);
      A(0,2) = 0;

      A(1,0) = (x[1][1] - x[0][1]);
      A(1,1) = (x[2][1] - x[0][1]);
      A(1,2) = 0;

      A(2,0) = 0; // (x[1][2] - x[0][2]);
      A(2,1) = 0; // (x[2][2] - x[0][2]);
      A(2,2) = 1.0;

      detJ = det(A);
      return detJ < 0.0;
    }

    template<>    bool  JacobianUtilImpl<STKMesh>::jacobian_matrix_2D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3], double *area_vector)
    {
      double x0[3] = {(x[1][0] - x[0][0]), (x[1][1] - x[0][1]), (x[1][2] - x[0][2])};
      double x1[3] = {(x[2][0] - x[0][0]), (x[2][1] - x[0][1]), (x[2][2] - x[0][2])};
      double v0[3] = {0,0,0};
      double *v = (area_vector ? area_vector : &v0[0]);
      Math::cross_3d(x0, x1, v);

      double vl = Math::norm_3d(v);
      if (vl == 0.0) vl = 1.0;
      A(0,0) = v[0]/vl;
      A(0,1) = v[1]/vl;
      A(0,2) = v[2]/vl;

      A(1,0) = x0[0];
      A(1,1) = x0[1];
      A(1,2) = x0[2];

      A(2,0) = x1[0];
      A(2,1) = x1[1];
      A(2,2) = x1[2];

      detJ =  det(A);
      return detJ < 0.0;
    }

    /// calculates dMetric_dA (@param dMdA) times dA/dx_n_i to get the gradient of the metric dMetric_dx_n_i, where
    /// x_n_i is the i'th coordinate of the n'th node in the element.
    /// @param indices are the 4 indices associated with this corner of the element, passed in from the grad_metric function
    /// Note: this is the dMetric_dx_n_i term associated with a particular corner, so there are up to nnode of these passed
    /// in from the grad_metric function
    /// @see jacobian_matrix_3D
    template<>    void  JacobianUtilImpl<STKMesh>::grad_util_2d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int /*spd*/, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < 3; j++)
          grad[i][j]=0.0;

      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      grad[indices[2]][0] += dMdA(0,1)*(+1); grad[indices[0]][0] += dMdA(0,1)*(-1);
      //grad[indices[3]][0] += dMdA(0,2)*(+1); grad[indices[0]][0] += dMdA(0,2)*(-1);

      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      grad[indices[2]][1] += dMdA(1,1)*(+1); grad[indices[0]][1] += dMdA(1,1)*(-1);
      //grad[indices[3]][1] += dMdA(1,2)*(+1); grad[indices[0]][1] += dMdA(1,2)*(-1);

      //grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      //grad[indices[2]][2] += dMdA(2,1)*(+1); grad[indices[0]][2] += dMdA(2,1)*(-1);
      //grad[indices[3]][2] += dMdA(2,2)*(+1); grad[indices[0]][2] += dMdA(2,2)*(-1);

    }

    template<>    void  JacobianUtilImpl<STKMesh>::grad_util(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int /*spd*/, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < 3; j++)
          grad[i][j]=0.0;

      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      grad[indices[2]][0] += dMdA(0,1)*(+1); grad[indices[0]][0] += dMdA(0,1)*(-1);
      grad[indices[3]][0] += dMdA(0,2)*(+1); grad[indices[0]][0] += dMdA(0,2)*(-1);

      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      grad[indices[2]][1] += dMdA(1,1)*(+1); grad[indices[0]][1] += dMdA(1,1)*(-1);
      grad[indices[3]][1] += dMdA(1,2)*(+1); grad[indices[0]][1] += dMdA(1,2)*(-1);

      grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      grad[indices[2]][2] += dMdA(2,1)*(+1); grad[indices[0]][2] += dMdA(2,1)*(-1);
      grad[indices[3]][2] += dMdA(2,2)*(+1); grad[indices[0]][2] += dMdA(2,2)*(-1);

    }

    template<>    void  JacobianUtilImpl<STKMesh>::grad_util_pyramid_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < spd; j++)
          grad[i][j]=0.0;

      const double h = 0.5;		/* h = 1 / (2*height) */

      // grad[i][j] += dMdA(j,m)*D[A(j,m),x[i,j]]

      /* Calculate M = A*inv(W). */
      //matr[0] = x1[0] - x0[0];  A(0,0)
      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      //matr[1] = x2[0] - x0[0];  A(0,1)
      grad[indices[2]][0] += dMdA(0,1)*(+1); grad[indices[0]][0] += dMdA(0,1)*(-1);
      //matr[2] = (2.0*x3[0] - x1[0] - x2[0])*h;  A(0,2)
      grad[indices[3]][0] += dMdA(0,2)*(+2.0*h); grad[indices[1]][0] += dMdA(0,2)*(-1.0*h); grad[indices[2]][0] += dMdA(0,2)*(-1.0*h);

      //matr[3] = x1[1] - x0[1];  A(1,0)
      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      //matr[4] = x2[1] - x0[1];  A(1,1)
      grad[indices[2]][1] += dMdA(1,1)*(+1); grad[indices[0]][1] += dMdA(1,1)*(-1);
      //matr[5] = (2.0*x3[1] - x1[1] - x2[1])*h;  A(1,2)
      grad[indices[3]][1] += dMdA(1,2)*(+2.0*h); grad[indices[1]][1] += dMdA(1,2)*(-1.0*h); grad[indices[2]][1] += dMdA(1,2)*(-1.0*h);

      //matr[6] = x1[2] - x0[2]; A(2,0)
      grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      //matr[7] = x2[2] - x0[2];  A(2,1)
      grad[indices[2]][2] += dMdA(2,1)*(+1); grad[indices[0]][2] += dMdA(2,1)*(-1);
      //matr[8] = (2.0*x3[2] - x1[2] - x2[2])*h;  A(2,2)
      grad[indices[3]][2] += dMdA(2,2)*(+2.0*h); grad[indices[1]][2] += dMdA(2,2)*(-1.0*h); grad[indices[2]][2] += dMdA(2,2)*(-1.0*h);
    }

    template<>    void  JacobianUtilImpl<STKMesh>::grad_util_tet_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < spd; j++)
          grad[i][j]=0.0;

      /* Calculate M = A*inv(W). */
      //matr[0] = x1[0] - x0[0];
      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      //matr[1] = (2.0*x2[0] - (x1[0] + x0[0]))*isqrt3;
      grad[indices[2]][0] += dMdA(0,1)*(+2.0*isqrt3); grad[indices[1]][0] += dMdA(0,1)*(-1.0*isqrt3); grad[indices[0]][0] += dMdA(0,1)*(-1.0*isqrt3);
      //matr[2] = (3.0*x3[0] - x2[0] - (x1[0] + x0[0]))*isqrt6;
      grad[indices[3]][0] += dMdA(0,2)*(+3.0*isqrt6); grad[indices[2]][0] += dMdA(0,2)*(-1.0*isqrt6); grad[indices[1]][0] += dMdA(0,2)*(-1.0*isqrt6); grad[indices[0]][0] += dMdA(0,2)*(-1.0*isqrt6);

      //matr[3] = x1[1] - x0[1];
      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      //matr[4] = (2.0*x2[1] - (x1[1] + x0[1]))*isqrt3;
      grad[indices[2]][1] += dMdA(1,1)*(+2.0*isqrt3); grad[indices[1]][1] += dMdA(1,1)*(-1.0*isqrt3); grad[indices[0]][1] += dMdA(1,1)*(-1.0*isqrt3);
      //matr[5] = (3.0*x3[1] - x2[1] - (x1[1] + x0[1]))*isqrt6;
      grad[indices[3]][1] += dMdA(1,2)*(+3.0*isqrt6); grad[indices[2]][1] += dMdA(1,2)*(-1.0*isqrt6); grad[indices[1]][1] += dMdA(1,2)*(-1.0*isqrt6); grad[indices[0]][1] += dMdA(1,2)*(-1.0*isqrt6);

      //matr[6] = x1[2] - x0[2];
      grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      //matr[7] = (2.0*x2[2] - (x1[2] + x0[2]))*isqrt3;
      grad[indices[2]][2] += dMdA(2,1)*(+2.0*isqrt3); grad[indices[1]][2] += dMdA(2,1)*(-1.0*isqrt3); grad[indices[0]][2] += dMdA(2,1)*(-1.0*isqrt3);
      //matr[8] = (3.0*x3[2] - x2[2] - (x1[2] + x0[2]))*isqrt6;
      grad[indices[3]][2] += dMdA(2,2)*(+3.0*isqrt6); grad[indices[2]][2] += dMdA(2,2)*(-1.0*isqrt6); grad[indices[1]][2] += dMdA(2,2)*(-1.0*isqrt6); grad[indices[0]][2] += dMdA(2,2)*(-1.0*isqrt6);

    }

    template<>    void  JacobianUtilImpl<STKMesh>::grad_util_wedge_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < spd; j++)
          grad[i][j]=0.0;

      /* Calculate M = A*inv(W). */
      //matr[0] = x1[0] - x0[0];
      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      //matr[1] = (2.0*x2[0] - (x1[0] + x0[0]))*isqrt3;
      grad[indices[2]][0] += dMdA(0,1)*(+2.0*isqrt3); grad[indices[1]][0] += dMdA(0,1)*(-1.0*isqrt3); grad[indices[0]][0] += dMdA(0,1)*(-1.0*isqrt3);
      //matr[2] = x3[0] - x0[0];
      grad[indices[3]][0] += dMdA(0,2)*(+1); grad[indices[0]][0] += dMdA(0,2)*(-1);

      //matr[3] = x1[1] - x0[1];
      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      //matr[4] = (2.0*x2[1] - (x1[1] + x0[1]))*isqrt3;
      grad[indices[2]][1] += dMdA(1,1)*(+2.0*isqrt3); grad[indices[1]][1] += dMdA(1,1)*(-1.0*isqrt3); grad[indices[0]][1] += dMdA(1,1)*(-1.0*isqrt3);
      //matr[5] = x3[1] - x0[1];
      grad[indices[3]][1] += dMdA(1,2)*(+1); grad[indices[0]][1] += dMdA(1,2)*(-1);

      //matr[6] = x1[2] - x0[2];
      grad[indices[1]][2] += dMdA(2,0)*(+1); grad[indices[0]][2] += dMdA(2,0)*(-1);
      //matr[7] = (2.0*x2[2] - (x1[2] + x0[2]))*isqrt3;
      grad[indices[2]][2] += dMdA(2,1)*(+2.0*isqrt3); grad[indices[1]][2] += dMdA(2,1)*(-1.0*isqrt3); grad[indices[0]][2] += dMdA(2,1)*(-1.0*isqrt3);
      //matr[8] = x3[2] - x0[2];
      grad[indices[3]][2] += dMdA(2,2)*(+1); grad[indices[0]][2] += dMdA(2,2)*(-1);

    }

    template<>    void  JacobianUtilImpl<STKMesh>::grad_util_tri_2d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int /*spd*/, const int *indices, const int /*nind*/)
    {
      for (int i=0; i < nnode; i++)
        for (int j=0; j < 3; j++)
          grad[i][j]=0.0;


      /* Calculate M = [A*inv(W) n] */
      //matr[0] = x[1][0] - x[0][0];
      grad[indices[1]][0] += dMdA(0,0)*(+1); grad[indices[0]][0] += dMdA(0,0)*(-1);
      //matr[1] = (2.0*x[2][0] - x[1][0] - x[0][0])*isqrt3;
      grad[indices[2]][0] += dMdA(0,1)*(+2.0*isqrt3); grad[indices[1]][0] += dMdA(0,1)*(-1.0*isqrt3); grad[indices[0]][0] += dMdA(0,1)*(-1.0*isqrt3);

      //matr[3] = x[1][1] - x[0][1];
      grad[indices[1]][1] += dMdA(1,0)*(+1); grad[indices[0]][1] += dMdA(1,0)*(-1);
      //matr[4] = (2.0*x[2][1] - x[1][1] - x[0][1])*isqrt3;
      grad[indices[2]][1] += dMdA(1,1)*(+2.0*isqrt3); grad[indices[1]][1] += dMdA(1,1)*(-1.0*isqrt3); grad[indices[0]][1] += dMdA(1,1)*(-1.0*isqrt3);

    }

    template<>    void  JacobianUtilImpl<STKMesh>::check_unhandled_topo(PerceptMesh& eMesh, const CellTopologyData * topology_data)
    {
      if (!m_use_approximate_quadratic_jacobian)
        {
          shards::CellTopology topology(topology_data);
          std::cout << "topology = " << topology.getName() << std::endl;
          VERIFY_MSG("unknown/unhandled topology in JacobianUtil: "+std::string(topology.getName()));
        }
      else
        {
          static int count = 0;
          if (count < 3)
            {
              shards::CellTopology topology(topology_data);
              if (eMesh.get_rank() == 0) 
                std::cout << "WARNING[percept::JacobianUtil]: topology [" << topology.getName() << "] is unhandled, and will use linear volume estimate" << std::endl;
              ++count;
            }
        }
    }
    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate()
    template<>    bool  JacobianUtilImpl<STKMesh>::operator()(double& m,  PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::FieldBase *coord_field,
                                  const CellTopologyData * topology_data )
    {
      EXCEPTWATCH;
      static DenseMatrix<3,3> J;
      int spatialDimension = eMesh.get_spatial_dim();

      //static Vector3D n;			// Surface normal for 2D objects

      int i=0;

      //static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = eMesh.get_cell_topology(element);

      const MyPairIterRelation v_i(eMesh, element, eMesh.node_rank() );
      m_num_nodes = v_i.size();

      const double *x2d[3] = {0,0,0};
      //const double *x3d[4] = {0,0,0,0};

#define VERTEX(vi)  static_cast<double*>(stk::mesh::field_data( *coord_field, vi.entity() ))

      switch(topology_data->key)
        {
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<6>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::ShellTriangle<3>::key:
        case shards::Triangle<3>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          x2d[0] = VERTEX(v_i[0]);
          x2d[1] = VERTEX(v_i[1]);
          x2d[2] = VERTEX(v_i[2]);
          if (spatialDimension == 3)
            metric_valid = jacobian_matrix_2D_in_3D(m, J, x2d);
          else
            metric_valid = jacobian_matrix_tri_2D(m, J, x2d);
          for (i = 0; i < 3; i++) { m_detJ[i] = m; m_J[i] = J; }
          break;

        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::ShellQuadrilateral<4>::key:
        case shards::Quadrilateral<4>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            x2d[0] =  VERTEX(v_i[locs_hex[i][0]]);
            x2d[1] =  VERTEX(v_i[locs_hex[i][1]]);
            x2d[2] =  VERTEX(v_i[locs_hex[i][2]]);
            if (spatialDimension == 3)
              metric_valid = jacobian_matrix_2D_in_3D(m_detJ[i], m_J[i], x2d);
            else
              metric_valid = jacobian_matrix_2D(m_detJ[i], m_J[i], x2d);
          }
          m = average_metrics(m_detJ, 4);
          break;

        // FIXME - tmp
        case shards::Tetrahedron<10>::key:
          metric_valid = true;
          for (i = 0; i < 10; i++) { m_detJ[i] = 1.0; identity(m_J[i]); }
          break;

        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<11>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Tetrahedron<4>::key:
          metric_valid = jacobian_matrix_tet_3D(m, m_J[0],
                                                VERTEX(v_i[0]),
                                                VERTEX(v_i[1]),
                                                VERTEX(v_i[2]),
                                                VERTEX(v_i[3]) );
          m_detJ[0] = m;
          for (i = 1; i < 4; i++) { m_detJ[i] = m; m_J[i] = m_J[0]; }
          break;

        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Pyramid<5>::key:
          {
            for (i = 0; i < 5; ++i) {
              metric_valid = jacobian_matrix_pyramid_3D_new(i,
                                                            m_detJ[i], m_J[i],
                                                            VERTEX(v_i[0]),
                                                            VERTEX(v_i[1]),
                                                            VERTEX(v_i[2]),
                                                            VERTEX(v_i[3]),
                                                            VERTEX(v_i[4]));
            }

            m = average_metrics(m_detJ, 5);
          }
          break;


        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Wedge<6>::key:
          for (i = 0; i < 6; ++i) {
            metric_valid = jacobian_matrix_wedge_3D(m_detJ[i], m_J[i],
                                                    VERTEX(v_i[locs_prism[i][0]]),
                                                    VERTEX(v_i[locs_prism[i][1]]),
                                                    VERTEX(v_i[locs_prism[i][2]]),
                                                    VERTEX(v_i[locs_prism[i][3]]));
          }
          m = average_metrics(m_detJ, 6);
          break;

        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
          check_unhandled_topo(eMesh, topology_data);
        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            metric_valid = jacobian_matrix_3D(m_detJ[i], m_J[i],
                                              VERTEX(v_i[locs_hex[i][0]]),
                                              VERTEX(v_i[locs_hex[i][1]]),
                                              VERTEX(v_i[locs_hex[i][2]]),
                                              VERTEX(v_i[locs_hex[i][3]]));
          }
          m = average_metrics(m_detJ, 8);
          break;


        case shards::Line<2>::key:
        case shards::Line<3>::key:
        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:

          if (1)
            {
              VERIFY_OP_ON(m_num_nodes, <=, 2, "m_num_nodes");
              const double *xv[2] = { VERTEX(v_i[0]), VERTEX(v_i[1]) };

              if (spatialDimension == 3)
                {
                  metric_valid = jacobian_matrix_1D_in_3D(m_detJ[0], m_J[0], xv);
                }
              else
                {
                  metric_valid = jacobian_matrix_1D_in_2D(m_detJ[0], m_J[0], xv);
                }
              m_J[1] = m_J[0];
              m_detJ[1] = m_detJ[0];
            }

          break;

          // unimplemented
        case shards::Node::key:
        case shards::Particle::key:
        case shards::Pentagon<5>::key:
        case shards::Hexagon<6>::key:

        default:
          check_unhandled_topo(eMesh, topology_data);
          break;

        } // end switch over element type

      return metric_valid;
    }

    /// modeled after code from Mesquite::IdealWeightMeanRatio::evaluate(), and TargetMetricUtil
    /// fills the mGrad member variable given the array of (member variable) m_dMetric_dA terms
    template<>    bool  JacobianUtilImpl<STKMesh>::grad_metric_util( PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::FieldBase */*coord_field*/,
                                         const CellTopologyData * topology_data )
    {
      static DenseMatrix<3,3> J;

      //static Vector3D n;			// Surface normal for 2D objects

      int i=0;

      //static const Vector3D d_con(1.0, 1.0, 1.0);

      bool metric_valid = false;
      if (!topology_data) topology_data = eMesh.get_cell_topology(element);

      const MyPairIterRelation v_i(eMesh, element, eMesh.node_rank() );
      m_num_nodes = v_i.size();

      const int indices_tri[3] = {0,1,2};
      const int indices_tet[4] = {0,1,2,3};
      switch(topology_data->key)
        {
        case shards::Triangle<3>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 3; i++) {
            grad_util_tri_2d(m_dMetric_dA[i], m_grad[i], 3, 2, indices_tri, 3);
          }
          break;

        case shards::Quadrilateral<4>::key:
          //n[0] = 0; n[1] = 0; n[2] = 1;
          for (i = 0; i < 4; ++i) {
            const int *indices_quad = locs_hex[i];
            grad_util_2d(m_dMetric_dA[i], m_grad[i], 4, 2, indices_quad, 3);
          }
          break;

        case shards::Tetrahedron<4>::key:
          for (i = 0; i < 4; ++i) {
            grad_util_tet_3d(m_dMetric_dA[i], m_grad[i], 4, 3, indices_tet, 4);
          }
          break;

        case shards::Pyramid<5>::key:
          for (i = 0; i < 4; ++i) {
            const int indices_pyr[4] = {i, (i+1)%4, (i+3)%4, 4};
            grad_util_pyramid_3d(m_dMetric_dA[i], m_grad[i], 5, 3, indices_pyr, 4);
          }
          // FIXME
          for ( i=0; i < 5; i++)
            {
              for (int j=0; j < 3; j++)
                m_grad[4][i][j] = 0.0;
            }
          for (int k=0; k < 4; k++)
            {
              for (i=0; i < 5; i++)
                {
                  for (int j=0; j < 3; j++)
                    m_grad[4][i][j] += m_grad[k][i][j]*0.25;
                }
            }
          break;

        case shards::Wedge<6>::key:
          for (i = 0; i < 6; ++i) {
            const int *indices_prism = locs_prism[i];
            grad_util_wedge_3d(m_dMetric_dA[i], m_grad[i], 6, 3, indices_prism, 4);
          }
          break;

        case shards::Hexahedron<8>::key:
          for (i = 0; i < 8; ++i) {
            const int *indices_hex = locs_hex[i];
            grad_util(m_dMetric_dA[i], m_grad[i], 8, 3, indices_hex, 4);
          }
          break;

        // FIXME - tmp
        case shards::Tetrahedron<10>::key:
          {
            check_unhandled_topo(eMesh, topology_data);
            metric_valid=true;
            for (int inode=0; inode < NNODES_MAX; ++inode)
              for (int j=0; j < NNODES_MAX; ++j)
                for (int k=0; k < 3; k++)
                  m_grad[inode][j][k] = 0.0;
          }
          break;

          // unimplemented
        case shards::Node::key:
        case shards::Particle::key:
        case shards::Line<2>::key:
        case shards::Line<3>::key:
        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:

        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:

        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:

        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<11>::key:

        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:

        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:

        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:

        case shards::Pentagon<5>::key:
        case shards::Hexagon<6>::key:

        default:
          check_unhandled_topo(eMesh, topology_data);
          break;

        } // end switch over element type

      return metric_valid;
    }

    template<>    void  JacobianUtilImpl<STKMesh>::edge_lengths(PerceptMesh& eMesh, stk::mesh::Entity element,
                                    double& min_edge_length, double& max_edge_length, double& ave_edge_length,
                                    stk::mesh::FieldBase *coord_field ,
                                    const CellTopologyData * topology_data_in  )
    {
      if (!topology_data_in) topology_data_in = eMesh.get_cell_topology(element);
      if (!coord_field) coord_field = eMesh.get_coordinates_field();
      ave_edge_length = eMesh.edge_length_ave(element, coord_field,  &min_edge_length, &max_edge_length, topology_data_in);
    }

    template<>    void  JacobianUtilImpl<STKMesh>::stretch_eigens(PerceptMesh& eMesh, stk::mesh::Entity element,
                                      double stretch_eigens[3],
                                      stk::mesh::FieldBase *coord_field ,
                                      const CellTopologyData * topology_data_in  )
    {
      if (!topology_data_in) topology_data_in = eMesh.get_cell_topology(element);
      if (!coord_field) coord_field = eMesh.get_coordinates_field();
      double averageDetJ = 0.0;
      bool valid = this->operator()(averageDetJ, eMesh, element, coord_field, topology_data_in);
      (void)valid;
      static DenseMatrix<3,3> Jt, Ui;
      for (int j = 0; j < 3; j++) stretch_eigens[j] = 0.0;
      for (int i=0; i < m_num_nodes; i++)
        {
          double eigen[3];
          Jt = transpose(m_J[i]);
          product(Jt, m_J[i], Ui);
          eigen_3x3(Ui, eigen);
          for (int j = 0; j < 3; j++) {
            if (eigen[j] < -1.e-8)
              {
                std::cout << "eigens not positive" << eigen[j];
              }
            VERIFY_OP_ON(eigen[j], >=, -1.e-8, "eigens not positive");
            stretch_eigens[j] += std::sqrt(std::fabs(eigen[j]))/double(m_num_nodes);
          }
        }

      for (int i = 1; i < 3; i++)
        {
          VERIFY_OP_ON(stretch_eigens[i], <=,  stretch_eigens[i-1], "eigens not sorted");
        }
    }

    template<>    double  JacobianUtilImpl<STKMesh>::getJacobianToVolumeScale(shards::CellTopology& cell_topo)
    {
      double volEqui = 1.0;
      switch(cell_topo.getKey() )
        {

          // Tet cells
        case shards::Tetrahedron<4>::key:
        case shards::Tetrahedron<8>::key:
        case shards::Tetrahedron<10>::key:
          volEqui = 1./6.;
          break;

          // Hex cells
        case shards::Hexahedron<8>::key:
        case shards::Hexahedron<20>::key:
        case shards::Hexahedron<27>::key:
          volEqui = 1.0;
          break;

          // Pyramid cells
        case shards::Pyramid<5>::key:
        case shards::Pyramid<13>::key:
        case shards::Pyramid<14>::key:
          volEqui = 1./6.;
          break;

          // Wedge cells
        case shards::Wedge<6>::key:
        case shards::Wedge<15>::key:
        case shards::Wedge<18>::key:
          volEqui = 1./2.;
          break;

        case shards::Triangle<3>::key:
        case shards::Triangle<4>::key:
        case shards::Triangle<6>::key:
          volEqui = 1./2.;
          break;

        case shards::Quadrilateral<4>::key:
        case shards::Quadrilateral<8>::key:
        case shards::Quadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellTriangle<3>::key:
        case shards::ShellTriangle<6>::key:
          volEqui = 1./2.;
          break;

        case shards::ShellQuadrilateral<4>::key:
        case shards::ShellQuadrilateral<8>::key:
        case shards::ShellQuadrilateral<9>::key:
          volEqui = 1.0;
          break;

        case shards::ShellLine<2>::key:
        case shards::ShellLine<3>::key:
        case shards::Beam<2>::key:
        case shards::Beam<3>::key:
          volEqui = 1.0;
          break;

        default:
          break;
        }//cell key
      return volEqui;
    }

  }


#endif
