// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef VolumeUtil_hpp
#define VolumeUtil_hpp

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)


#include <percept/PerceptMesh.hpp>
#include <percept/math/DenseMatrix.hpp>

  namespace percept {

    class VolumeUtil
    {

    public:
      typedef double Vec3D[3];
      Vec3D mCoords[4];

      enum { NELEM_TYPES = 10, NNODES_MAX = 10 };

      double   m_detJ[NNODES_MAX];
      DenseMatrix<3,3> m_J[NNODES_MAX];
      DenseMatrix<3,3> m_dMetric_dA[NNODES_MAX];
      double m_grad[NNODES_MAX][NNODES_MAX][3];
      int m_num_nodes;
      bool m_use_approximate_quadratic_volume;

      VolumeUtil(bool use_approximate_quadratic_volume=true) :
        m_num_nodes(0),
        m_use_approximate_quadratic_volume(use_approximate_quadratic_volume)
      {
      }

      double average_metrics(const double detJ[NNODES_MAX], const unsigned n)
      {
        double sum=0.0;
        for (unsigned i=0; i < n; i++) sum += detJ[i];
        return sum/double(n);
      }

      // note: returns the Jacobian - must multiply by getJacobianToVolumeScale to get volume
      bool operator()(double& averageJ, PerceptMesh& eMesh, stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field,
                      const CellTopologyData * topology_data_in = 0 );
      bool operator()(double& averageJ, stk::mesh::Entity element, const stk::mesh::FieldBase *coord_field,
                      const CellTopologyData * topology_data_in = 0 )
      {
        PerceptMesh eMesh(&coord_field->mesh_meta_data(), &coord_field->get_mesh(), true);
        return this->operator()(averageJ, eMesh, element, coord_field, topology_data_in);
      }

      double getJacobianToVolumeScale(shards::CellTopology& cell_topo);


    private:
      void check_unhandled_topo(PerceptMesh& eMesh, const CellTopologyData * topology_data);

      inline bool volume_matrix_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
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
      bool volume_matrix_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);


      // following is from Mesquite.MeanRatioFunctions

      /*  1.0/sqrt(3.0)*/
#define isqrt3  5.77350269189625797959429519858e-01
      /*  2.0/sqrt(3.0)*/
#define tisqrt3  1.15470053837925159591885903972e+00
      /*  1.0/sqrt(6.0)*/
#define isqrt6   4.08248290463863052509822647505e-01
      /*  3.0/sqrt(6.0)*/
#define tisqrt6  1.22474487139158915752946794252e+00

      inline bool volume_matrix_tet_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        return volume_matrix_3D(detJ, A, x0, x1, x2, x3);
      }

      inline bool volume_matrix_wedge_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        return volume_matrix_3D(detJ, A, x0, x1, x2, x3);
      }

      inline bool volume_matrix_pyramid_3D(double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3)
      {
        return volume_matrix_3D(detJ, A, x0, x1, x2, x3);
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
      inline bool volume_matrix_pyramid_3D_new(const int ibasis, double &detJ, DenseMatrix<3,3>& A, const double *x0, const double *x1, const double *x2, const double *x3, const double *x4)
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

      inline bool volume_matrix_tri_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3])
      {
        return volume_matrix_2D(detJ, A, x);
      }

      bool volume_matrix_2D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);
      bool volume_matrix_1D_in_3D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);
      bool volume_matrix_1D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);

    };
  }

#endif
#endif
