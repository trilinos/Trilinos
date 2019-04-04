// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef sgrid_metric_helpers_hpp
#define sgrid_metric_helpers_hpp

    KOKKOS_INLINE_FUNCTION
    void identity_dev( double I[3][3] )
    {
      I[0][0] = 1.0;
      I[0][1] = 0.0;
      I[0][2] = 0.0;

      I[1][0] = 0.0;
      I[1][1] = 1.0;
      I[1][2] = 0.0;

      I[2][0] = 0.0;
      I[2][1] = 0.0;
      I[2][2] = 1.0;
    }


    KOKKOS_INLINE_FUNCTION
    void matrix_transpose(double m[3][3], double result[3][3])
    {
        for(unsigned i=0;i<3;i++)
            for(unsigned j=0;j<3;j++)
                result[i][j]=m[j][i];
    }

    KOKKOS_INLINE_FUNCTION
    double det(double m[3][3])  {
        return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
    } //double det

    KOKKOS_INLINE_FUNCTION
    bool jacobian_matrix_3D(double &detJ, double A[3][3], const double * x0,
            const double * x1, const double * x2, const double * x3)  {
        A[0][0] = (x1[0] - x0[0]);
        A[0][1] = (x2[0] - x0[0]);
        A[0][2] = (x3[0] - x0[0]);

        A[1][0] = (x1[1] - x0[1]);
        A[1][1] = (x2[1] - x0[1]);
        A[1][2] = (x3[1] - x0[1]);

        A[2][0] = (x1[2] - x0[2]);
        A[2][1] = (x2[2] - x0[2]);
        A[2][2] = (x3[2] - x0[2]);

        detJ = det(A);

        return detJ < 0.0;
    } //double jacobian_matrix_3D

    KOKKOS_INLINE_FUNCTION
    bool sGridJacobianUtil(double detJ[8],
    double coords[8][3],Kokkos::Array<double[3][3], 8>& J)  {
        const int locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                4, 2 }, { 7, 5, 6, 3 } };

        bool metric_invalid = false;

        for (int i = 0; i < 8; ++i) {
            bool mi = jacobian_matrix_3D(detJ[i], J[i], coords[locs_hex_dev[i][0]],
                    coords[locs_hex_dev[i][1]], coords[locs_hex_dev[i][2]],
                    coords[locs_hex_dev[i][3]]);

            metric_invalid = metric_invalid || mi;
        }

        return metric_invalid;
    } //bool sGriJacobianUtil

    KOKKOS_INLINE_FUNCTION
    void matrix_inverse(const double * m[3], const double detm,
            double result[3][3])  {
        const double detInv = 1.0 / detm;
        result[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * detInv;
        result[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * detInv;
        result[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * detInv;

        result[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * detInv;
        result[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * detInv;
        result[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * detInv;

        result[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * detInv;
        result[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * detInv;
        result[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * detInv;
    }

    KOKKOS_INLINE_FUNCTION
    void matrix_product(const double * x[3], const double y[3][3],
            double z[3][3])  {
#define R(i,j)   z[i][j] = x[i][0]*y[0][j] + x[i][1]*y[1][j] + x[i][2]*y[2][j]
        R(0, 0);
        R(0, 1);
        R(0, 2);

        R(1, 0);
        R(1, 1);
        R(1, 2);

        R(2, 0);
        R(2, 1);
        R(2, 2);
#undef R
    }

    KOKKOS_INLINE_FUNCTION
    void matrix_product_mutator(double mutated[3][3], double mutator[3][3]
            )  {
        double temp[3][3];
#define R(i,j)   temp[i][j] = mutated[i][0]*mutator[0][j] + mutated[i][1]*mutator[1][j] + mutated[i][2]*mutator[2][j]
        R(0, 0);
        R(0, 1);
        R(0, 2);

        R(1, 0);
        R(1, 1);
        R(1, 2);

        R(2, 0);
        R(2, 1);
        R(2, 2);
#undef R
        for(unsigned i=0;i<3;i++)
            for(unsigned j=0;j<3;j++)
                mutated[i][j]=temp[i][j];
    }

    KOKKOS_INLINE_FUNCTION
    double matrix_sqr_Frobenius(const double m[3][3])  {
        double sum = 0.0;
#define R(i,j)   sum += m[i][j]*m[i][j]
        R(0, 0);
        R(0, 1);
        R(0, 2);

        R(1, 0);
        R(1, 1);
        R(1, 2);

        R(2, 0);
        R(2, 1);
        R(2, 2);
#undef R
        return sum;
    }

    KOKKOS_INLINE_FUNCTION
    void transpose_adj_dev(double m[][3], double result[3][3], double scalar=1.0 )
    {
      result[0][0] = scalar*(m[1][1]*m[2][2] - m[1][2]*m[2][1]);
      result[0][1] = scalar*(m[1][2]*m[2][0] - m[1][0]*m[2][2]);
      result[0][2] = scalar*(m[1][0]*m[2][1] - m[1][1]*m[2][0]);

      result[1][0] = scalar*(m[0][2]*m[2][1] - m[0][1]*m[2][2]);
      result[1][1] = scalar*(m[0][0]*m[2][2] - m[0][2]*m[2][0]);
      result[1][2] = scalar*(m[0][1]*m[2][0] - m[0][0]*m[2][1]);

      result[2][0] = scalar*(m[0][1]*m[1][2] - m[0][2]*m[1][1]);
      result[2][1] = scalar*(m[0][2]*m[1][0] - m[0][0]*m[1][2]);
      result[2][2] = scalar*(m[0][0]*m[1][1] - m[0][1]*m[1][0]);
    }


    KOKKOS_INLINE_FUNCTION
    void  grad_util(const double dMdA[3][3], double grad[3][3], const unsigned nnode, const unsigned *indices)
    {
      for (unsigned i=0; i < nnode; i++)
        for (unsigned j=0; j < 3; j++)
          grad[i][j]=0.0;

      grad[indices[1]][0] += dMdA[0][0]*(+1); grad[indices[0]][0] += dMdA[0][0]*(-1);
      grad[indices[2]][0] += dMdA[0][1]*(+1); grad[indices[0]][0] += dMdA[0][1]*(-1);
      grad[indices[3]][0] += dMdA[0][2]*(+1); grad[indices[0]][0] += dMdA[0][2]*(-1);

      grad[indices[1]][1] += dMdA[1][0]*(+1); grad[indices[0]][1] += dMdA[1][0]*(-1);
      grad[indices[2]][1] += dMdA[1][1]*(+1); grad[indices[0]][1] += dMdA[1][1]*(-1);
      grad[indices[3]][1] += dMdA[1][2]*(+1); grad[indices[0]][1] += dMdA[1][2]*(-1);

      grad[indices[1]][2] += dMdA[2][0]*(+1); grad[indices[0]][2] += dMdA[2][0]*(-1);
      grad[indices[2]][2] += dMdA[2][1]*(+1); grad[indices[0]][2] += dMdA[2][1]*(-1);
      grad[indices[3]][2] += dMdA[2][2]*(+1); grad[indices[0]][2] += dMdA[2][2]*(-1);

    }

    KOKKOS_INLINE_FUNCTION
    bool grad_metric_util(double dMetric_dA[8][3][3],double grads[8][8][3], const unsigned locs_hex_dev[8][4])
    {
        bool metric_valid  = false;
        for (unsigned i = 0; i < 8; ++i) {
          const unsigned *indices_hex = locs_hex_dev[i];
          grad_util(dMetric_dA[i], grads[i], 8, indices_hex);
        }

        return metric_valid;
    }
#endif
