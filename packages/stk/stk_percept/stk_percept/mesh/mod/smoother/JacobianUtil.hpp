/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef JacobianUtil_hpp
#define JacobianUtil_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__)


#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/math/DenseMatrix.hpp>

namespace stk {
  namespace percept {

    class JacobianUtil
    {

    public:
      typedef double Vec3D[3];
      Vec3D mCoords[4];

      enum { NELEM_TYPES = 10, NNODES_MAX = 8 };

      double   m_detJ[NNODES_MAX];
      DenseMatrix<3,3> m_J[NNODES_MAX];
      DenseMatrix<3,3> m_dMetric_dA[NNODES_MAX];
      double m_grad[NNODES_MAX][NNODES_MAX][3];
      int m_num_nodes;
      bool m_scale_to_unit;

      JacobianUtil() :
        m_num_nodes(0),
        m_scale_to_unit(false)
      {
      }

      double average_metrics(const double detJ[NNODES_MAX], const unsigned n)
      {
        double sum=0.0;
        for (unsigned i=0; i < n; i++) sum += detJ[i];
        return sum/double(n);
      }

      bool operator()(double& averageJ, PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::FieldBase *coord_field,
                      const CellTopologyData * topology_data_in = 0 );

      bool grad_metric_util( PerceptMesh& eMesh, stk::mesh::Entity element, stk::mesh::FieldBase *coord_field,
                             const CellTopologyData * topology_data );

    private:
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
        //if (m_scale_to_unit) scale_to_unit(A);

        detJ = det(A);
        return detJ < 0.0;
      }
      bool jacobian_matrix_2D(double &detJ, DenseMatrix<3,3>& A, const double *x[3]);

      void grad_util_2d(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);
      void grad_util(const DenseMatrix<3,3>& dMdA, double grad[NNODES_MAX][3], const int nnode, const int spd, const int *indices, const int nind);

    };
  }
}

#endif
#endif
