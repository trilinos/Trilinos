// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef JacobianUtil_hpp
#define JacobianUtil_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)


#include <percept/PerceptMesh.hpp>
#include <percept/MeshType.hpp>
#include <percept/math/DenseMatrix.hpp>

  namespace percept {

    template<typename MeshType>
    class JacobianUtilBase
    {

    public:
      typedef double Vec3D[3];
      Vec3D mCoords[4];

      //enum { NELEM_TYPES = 10, NNODES_MAX = 10 };
      enum { NELEM_TYPES =  MeshType::NELEM_TYPES, NNODES_MAX = 10 };

      double   m_detJ[NNODES_MAX];
      DenseMatrix<3,3> m_J[NNODES_MAX];
      DenseMatrix<3,3> m_dMetric_dA[NNODES_MAX];
      double m_grad[NNODES_MAX][NNODES_MAX][3];
      int m_num_nodes;
      bool m_use_approximate_quadratic_jacobian;

      JacobianUtilBase(bool use_approximate_quadratic_jacobian=true) :
        m_num_nodes(0),
        m_use_approximate_quadratic_jacobian(use_approximate_quadratic_jacobian)
      {
      }

      double average_metrics(const double detJ[NNODES_MAX], const unsigned n)
      {
        double sum=0.0;
        for (unsigned i=0; i < n; i++) sum += detJ[i];
        return sum/double(n);
      }
    };


    template<typename MeshType>
    class JacobianUtilImpl : public JacobianUtilBase<MeshType>
    {

    public:

      using Base =  JacobianUtilBase<MeshType>;

      using  typename Base::Vec3D;
      using  Base::NELEM_TYPES;
      using  Base::NNODES_MAX;
      using  Base::m_detJ;
      using  Base::m_J;
      using  Base::m_dMetric_dA;
      using  Base::m_grad;
      using  Base::m_num_nodes;
      using  Base::m_use_approximate_quadratic_jacobian;


      JacobianUtilImpl(bool use_approximate_quadratic_jacobian=true) : Base(use_approximate_quadratic_jacobian)
      {
      }

      bool operator()(double& averageJ, PerceptMesh& eMesh, typename MeshType::MTElement element, typename MeshType::MTField *coord_field,
                      const typename MeshType::MTCellTopology * topology_data_in = 0 );

      bool operator()(double& averageJ,  typename MeshType::MTElement element, typename MeshType::MTField *coord_field,
                      const typename MeshType::MTCellTopology * topology_data_in = 0 )
      {
        PerceptMesh eMesh(&coord_field->mesh_meta_data(), &coord_field->get_mesh(), true);
        return this->operator()(averageJ, eMesh, element, coord_field, topology_data_in);
      }

      bool grad_metric_util( PerceptMesh& eMesh, typename MeshType::MTElement element, typename MeshType::MTField *coord_field,
                             const typename MeshType::MTCellTopology * topology_data );

      double getJacobianToVolumeScale(shards::CellTopology& cell_topo);

      /// compute an approximate element diameter (diameter of circumscribing sphere) by computing
      ///   max of distance between all pairs of vertices
      /// The average of distance between vertex pairs is returned in @param ave_pair_length_returned
      double approx_diameter(PerceptMesh& eMesh, stk::mesh::Entity element, double& ave_pair_length_returned,
                             stk::mesh::FieldBase *coord_field = 0,
                             const CellTopologyData * topology_data_in = 0 );

      // element "h" based on gradient of a given field - "h" is computed by the "h" in the direction
      //   of the field's gradient - uses the formula h = || grad u || / || J^-1 grad u ||
      //   If || grad u || ~0, it returns approx_diameter()
      double grad_based_diameter(PerceptMesh& eMesh, stk::mesh::Entity element,
                                 stk::mesh::FieldBase *field_for_grad,
                                 stk::mesh::FieldBase *coord_field = 0,
                                 const CellTopologyData * topology_data_in = 0 );

      /// compute edge length min, max and average between pairs of vertices that form element edges
      void edge_lengths(PerceptMesh& eMesh, stk::mesh::Entity element,
                        double& min_edge_length, double& max_edge_length, double& ave_edge_length,
                        stk::mesh::FieldBase *coord_field = 0,
                        const CellTopologyData * topology_data_in = 0 );

      /// return sorted (largest first) eigenvalues of U (stretch matrix) of the polar decomposition
      ///   of Jacobian, J = R U, where R is a rotation.  These represent stretch w.r.t. principal
      ///   axes of the element, and are thus somewhat representative of mesh parameter.  Here, J is
      ///   the average J from the vertex-based (corner-based) Jacobians.
      /// Uses the jacobian w.r.t. "equilateral" reference elements, i.e. unit edge lengths
      void stretch_eigens(PerceptMesh& eMesh, stk::mesh::Entity element,
                          double stretch_eigens[3],
                          stk::mesh::FieldBase *coord_field = 0,
                          const CellTopologyData * topology_data_in = 0 );

      bool jacobian_matrix_2D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3], double *normal=0);

    private:
      void check_unhandled_topo(PerceptMesh& eMesh, const CellTopologyData * topology_data);

      inline bool jacobian_matrix_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        A(0,0) = (x1[0] - x0[0]);
        A(0,1) = (x2[0] - x0[0]);
        A(0,2) = (x3[0] - x0[0]);

        A(1,0) = (x1[1] - x0[1]);
        A(1,1) = (x2[1] - x0[1]);
        A(1,2) = (x3[1] - x0[1]);

        A(2,0) = (x1[2] - x0[2]);
        A(2,1) = (x2[2] - x0[2]);
        A(2,2) = (x3[2] - x0[2]);

        detJ = det(A);
        return detJ < 0.0;
      }
      bool jacobian_matrix_1D_in_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[2]);
      bool jacobian_matrix_1D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[2]);

      bool jacobian_matrix_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);

      void grad_util_2d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);

      void grad_util_tri_2d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util_tet_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util_wedge_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util_pyramid_3d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);

      // following is from Mesquite.MeanRatioFunctions

      /*  1.0/sqrt(3.0)*/
#define isqrt3  5.77350269189625797959429519858e-01
      /*  2.0/sqrt(3.0)*/
#define tisqrt3  1.15470053837925159591885903972e+00
      /*  1.0/sqrt(6.0)*/
#define isqrt6   4.08248290463863052509822647505e-01
      /*  3.0/sqrt(6.0)*/
#define tisqrt6  1.22474487139158915752946794252e+00

      inline bool jacobian_matrix_tet_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        /*****************************************************************************/
        /* The following set of functions reference tetrahedral elements to a        */
        /* regular tetrahedron.  They are used when assessing the quality of a       */
        /* tetrahedral element.  A zero return value indicates success, while        */
        /* a nonzero value indicates failure.                                        */
        /*****************************************************************************/
        double *matr = &A(0,0);
        double f;

        /* Calculate M = A*inv(W). */
        f       = x1[0] + x0[0];
        matr[0] = x1[0] - x0[0];
        matr[1] = (2.0*x2[0] - f)*isqrt3;
        matr[2] = (3.0*x3[0] - x2[0] - f)*isqrt6;

        f       = x1[1] + x0[1];
        matr[3] = x1[1] - x0[1];
        matr[4] = (2.0*x2[1] - f)*isqrt3;
        matr[5] = (3.0*x3[1] - x2[1] - f)*isqrt6;

        f       = x1[2] + x0[2];
        matr[6] = x1[2] - x0[2];
        matr[7] = (2.0*x2[2] - f)*isqrt3;
        matr[8] = (3.0*x3[2] - x2[2] - f)*isqrt6;

        /* Calculate det(M). */
        detJ = det(A);
        return detJ < 0.0;
      }

      inline bool jacobian_matrix_wedge_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {

        /*********************************************************************
         * Reference tetrahedral elements to corners of an ideal wedge element.
         * Vertices should be ordered such that the first three vertices form
         * the ideally-equaliteral face of the tetrahedron (the end of the
         * wedge) and the first and fourth vertices form the tetrahedral edge
         * orthogonal to the ideally-equaliteral face (the lateral edge of the
         * edge.)
         *      1  1/2        0                     1  -1/sqrt(3)  0
         * W =  0  sqrt(3)/2  0           inv(W) =  0   2/sqrt(3)  0
         *      0  0          1                     0   0          1
         *
         *********************************************************************/
        double *matr = &A(0,0);

        /* Calculate M = A*inv(W). */
        matr[0] = x1[0] - x0[0];
        matr[1] = isqrt3 * (2 * x2[0] - x1[0] - x0[0]);
        matr[2] = x3[0] - x0[0];

        matr[3] = x1[1] - x0[1];
        matr[4] = isqrt3 * (2 * x2[1] - x1[1] - x0[1]);
        matr[5] = x3[1] - x0[1];

        matr[6] = x1[2] - x0[2];
        matr[7] = isqrt3 * (2 * x2[2] - x1[2] - x0[2]);
        matr[8] = x3[2] - x0[2];
        /* Calculate det(M). */
        detJ = det(A);
        return detJ < 0.0;
      }

      inline bool jacobian_matrix_pyramid_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        /*********************************************************************
         * Reference tetrahedral elements to halves of an ideal pyramid.
         * Vertices should be ordered such that the edge between the 2nd and
         * 3rd vertices is the one longer edge in the reference tetrahedron
         * (the diagonal of the base of the pyramid of which the tetrahedron
         * is one half of).
         *      1  0  1/2                     1  0 -1/sqrt(2)
         * W =  0  1  1/2           inv(W) =  0  1 -1/sqrt(2)
         *      0  0  1/sqrt(2)               0  0  sqrt(2)
         *
         *********************************************************************/

        const double h = 0.5;		/* h = 1 / (2*height) */

        double *matr = &A(0,0);

        /* Calculate M = A*inv(W). */
        matr[0] = x1[0] - x0[0];
        matr[1] = x2[0] - x0[0];
        matr[2] = (2.0*x3[0] - x1[0] - x2[0])*h;

        matr[3] = x1[1] - x0[1];
        matr[4] = x2[1] - x0[1];
        matr[5] = (2.0*x3[1] - x1[1] - x2[1])*h;

        matr[6] = x1[2] - x0[2];
        matr[7] = x2[2] - x0[2];
        matr[8] = (2.0*x3[2] - x1[2] - x2[2])*h;
        /* Calculate det(M). */
        detJ = det(A);
        return detJ < 0.0;
      }

      inline double XI(int I, int J, const double *x0, const double *x1, const double *x2, const double *x3, const double *x4)
      {
        switch (I-1)
          {
          case 0: return x0[J-1];
          case 1: return x1[J-1];
          case 2: return x2[J-1];
          case 3: return x3[J-1];
          case 4: return x4[J-1];
          }
        return 0;
      }

      /** Corrected basis functions, derivatives, and Jacobian at the node @param ibasis of the pyramid.
          @see cs.brown.edu/cgc/cgc98/final/final19.ps and   www.global-sci.org/aamm/freedownload/32-131.pdf
          Reference element is [0,1]x[0,1]x[0,1/2]
          See pyr-correct.nb in this directory
      */
      inline bool jacobian_matrix_pyramid_3D_new(const int ibasis, double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3, const double *x4)
      {
        // parametric coordinates at vertices
        double rst[5][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{.5,.5,.49999}};  //note: avoid singularity at top vertex

        double r = rst[ibasis][0];
        double s = rst[ibasis][1];
        double t = rst[ibasis][2];

        // bases
        double bases[] = {-(((-1 + r + t)*(-1 + s + t))/(-1 + 2*t)),
                          ((r - t)*(-1 + s + t))/(-1 + 2*t),
                          -(((r - t)*(s - t))/(-1 + 2*t)),
                          ((s - t)*(-1 + r + t))/(-1 + 2*t),
                          2*t};
        (void)bases;

#define List(x,y,z) {x,y,z}
#define Power(x,y) (y==2?(x)*(x):std::pow(x,y))
#define xi(I,J) XI(I,J,x0,x1,x2,x3,x4)

        // Jacobian at (r,s,t)[ibasis] for top vertex (1/2,1/2,1/2)
        double dxidxij_top_vertex[3][3] =
          List(List((-xi(1,1) + xi(2,1) + xi(3,1) - xi(4,1))/2.,(-xi(1,1) - xi(2,1) + xi(3,1) + xi(4,1))/2.,(-xi(1,1) - xi(2,1) - xi(3,1) - xi(4,1) + 4*xi(5,1))/2.),
               List((-xi(1,2) + xi(2,2) + xi(3,2) - xi(4,2))/2.,(-xi(1,2) - xi(2,2) + xi(3,2) + xi(4,2))/2.,(-xi(1,2) - xi(2,2) - xi(3,2) - xi(4,2) + 4*xi(5,2))/2.),
               List((-xi(1,3) + xi(2,3) + xi(3,3) - xi(4,3))/2.,(-xi(1,3) - xi(2,3) + xi(3,3) + xi(4,3))/2.,(-xi(1,3) - xi(2,3) - xi(3,3) - xi(4,3) + 4*xi(5,3))/2.));

        // at base vertices [0,1]x[0,1] t=0
        double dxidxij_base_vertices[3][3] =
          List(List((-1 + s)*xi(1,1) - (-1 + s)*xi(2,1) + s*(xi(3,1) - xi(4,1)),(-1 + r)*xi(1,1) + xi(4,1) - r*(xi(2,1) - xi(3,1) + xi(4,1)),
                    -xi(2,1) + r*(-1 + 2*s)*(xi(1,1) - xi(2,1) + xi(3,1) - xi(4,1)) - xi(4,1) + s*(-xi(1,1) + xi(2,1) - xi(3,1) + xi(4,1)) + 2*xi(5,1)),
               List((-1 + s)*xi(1,2) - (-1 + s)*xi(2,2) + s*(xi(3,2) - xi(4,2)),(-1 + r)*xi(1,2) + xi(4,2) - r*(xi(2,2) - xi(3,2) + xi(4,2)),
                    -xi(2,2) + r*(-1 + 2*s)*(xi(1,2) - xi(2,2) + xi(3,2) - xi(4,2)) - xi(4,2) + s*(-xi(1,2) + xi(2,2) - xi(3,2) + xi(4,2)) + 2*xi(5,2)),
               List((-1 + s)*xi(1,3) - (-1 + s)*xi(2,3) + s*(xi(3,3) - xi(4,3)),(-1 + r)*xi(1,3) + xi(4,3) - r*(xi(2,3) - xi(3,3) + xi(4,3)),
                    -xi(2,3) + r*(-1 + 2*s)*(xi(1,3) - xi(2,3) + xi(3,3) - xi(4,3)) - xi(4,3) + s*(-xi(1,3) + xi(2,3) - xi(3,3) + xi(4,3)) + 2*xi(5,3))) ;

        // at general point
        double dxidxij[3][3] =
          List(List((-((-1 + s + t)*xi(1,1)) + (-1 + s + t)*xi(2,1) + (-s + t)*(xi(3,1) - xi(4,1)))/(-1 + 2*t),
                    (-((-1 + r + t)*xi(1,1)) + r*xi(2,1) - r*xi(3,1) - xi(4,1) + r*xi(4,1) + t*(-xi(2,1) + xi(3,1) + xi(4,1)))/(-1 + 2*t),
                    -((r*xi(1,1) + xi(2,1) - r*xi(2,1) + r*xi(3,1) - (-1 + 2*r)*s*(xi(1,1) - xi(2,1) + xi(3,1) - xi(4,1)) + xi(4,1) - r*xi(4,1) -
                       2*t*(xi(1,1) + xi(2,1) + xi(3,1) + xi(4,1) - 4*xi(5,1)) + 2*Power(t,2)*(xi(1,1) + xi(2,1) + xi(3,1) + xi(4,1) - 4*xi(5,1)) - 2*xi(5,1))/
                      Power(1 - 2*t,2))),List((-((-1 + s + t)*xi(1,2)) + (-1 + s + t)*xi(2,2) + (-s + t)*(xi(3,2) - xi(4,2)))/(-1 + 2*t),
                                              (-((-1 + r + t)*xi(1,2)) + r*xi(2,2) - r*xi(3,2) - xi(4,2) + r*xi(4,2) + t*(-xi(2,2) + xi(3,2) + xi(4,2)))/(-1 + 2*t),
                                              -((r*xi(1,2) + xi(2,2) - r*xi(2,2) + r*xi(3,2) - (-1 + 2*r)*s*(xi(1,2) - xi(2,2) + xi(3,2) - xi(4,2)) + xi(4,2) - r*xi(4,2) -
                                                 2*t*(xi(1,2) + xi(2,2) + xi(3,2) + xi(4,2) - 4*xi(5,2)) + 2*Power(t,2)*(xi(1,2) + xi(2,2) + xi(3,2) + xi(4,2) - 4*xi(5,2)) - 2*xi(5,2))/
                                                Power(1 - 2*t,2))),List((-((-1 + s + t)*xi(1,3)) + (-1 + s + t)*xi(2,3) + (-s + t)*(xi(3,3) - xi(4,3)))/(-1 + 2*t),
                                                                        (-((-1 + r + t)*xi(1,3)) + r*xi(2,3) - r*xi(3,3) - xi(4,3) + r*xi(4,3) + t*(-xi(2,3) + xi(3,3) + xi(4,3)))/(-1 + 2*t),
                                                                        -((r*xi(1,3) + xi(2,3) - r*xi(2,3) + r*xi(3,3) - (-1 + 2*r)*s*(xi(1,3) - xi(2,3) + xi(3,3) - xi(4,3)) + xi(4,3) - r*xi(4,3) -
                                                                           2*t*(xi(1,3) + xi(2,3) + xi(3,3) + xi(4,3) - 4*xi(5,3)) + 2*Power(t,2)*(xi(1,3) + xi(2,3) + xi(3,3) + xi(4,3) - 4*xi(5,3)) - 2*xi(5,3))/
                                                                          Power(1 - 2*t,2))));

        (void)dxidxij;

        /* Calculate M = A */
        for (int i=0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
              {
                if (ibasis == 4)
                  A(i,j) = dxidxij_top_vertex[i][j];
                else
                  A(i,j) = dxidxij_base_vertices[i][j];
              }
          }

        /* set inv(W) */
        DenseMatrix<3,3> WI;
        WI.zero();
        WI(0,0) = 1;
        WI(0,2) = -1./std::sqrt(2.);
        WI(1,1) = 1;
        WI(1,2) = -1./std::sqrt(2.);
        WI(2,2) = std::sqrt(2.);

        A = A * WI;

        /* Calculate det(M). */
        detJ = det(A);
        return detJ < 0.0;
      }

      inline bool jacobian_matrix_tri_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
      {

        /*****************************************************************************/
        /* The following set of functions reference triangular elements to an        */
        /* equilateral triangle in the plane defined by the normal.  They are        */
        /* used when assessing the quality of a triangular element.  A zero          */
        /* return value indicates success, while a nonzero value indicates failure.  */
        /*****************************************************************************/
        double *matr = &A(0,0);

        /* Calculate M = [A*inv(W) n] */
        matr[0] = x[1][0] - x[0][0];
        matr[1] = (2.0*x[2][0] - x[1][0] - x[0][0])*isqrt3;
        matr[2] = 0;

        matr[3] = x[1][1] - x[0][1];
        matr[4] = (2.0*x[2][1] - x[1][1] - x[0][1])*isqrt3;
        matr[5] = 0;

        matr[6] = 0; // x[1][2] - x[0][2];
        matr[7] = 0; // (2.0*x[2][2] - x[1][2] - x[0][2])*isqrt3;
        matr[8] = 1;

        /* Calculate det(M). */
        detJ = det(A);
        return detJ < 0.0;
      }

    };
    using JacobianUtil = JacobianUtilImpl<STKMesh>;

  }

#endif
#endif
