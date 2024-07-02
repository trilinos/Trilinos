// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER

namespace Compadre {

    //! Metric factor (det(G)) at any point in the local chart
    KOKKOS_INLINE_FUNCTION
    double MetricFactor(const scratch_vector_type a_, const double h, const double u1, const double u2) {
        double det_g = 1.0;
        double q1 = 0;
        double q2 = 0;
        switch (a_.extent(0)) {
            case 45:
                q1 += (1.0/5040.0)*a_[36]*std::pow(u1, 7)/std::pow(h, 8) + (1.0/720.0)*a_[37]*std::pow(u1, 6)*u2/std::pow(h, 8) + (1.0/240.0)*a_[38]*std::pow(u1, 5)*std::pow(u2, 2)/std::pow(h, 8) + (1.0/144.0)*a_[39]*std::pow(u1, 4)*std::pow(u2, 3)/std::pow(h, 8) + (1.0/144.0)*a_[40]*std::pow(u1, 3)*std::pow(u2, 4)/std::pow(h, 8) + (1.0/240.0)*a_[41]*std::pow(u1, 2)*std::pow(u2, 5)/std::pow(h, 8) + (1.0/720.0)*a_[42]*u1*std::pow(u2, 6)/std::pow(h, 8) + (1.0/5040.0)*a_[43]*std::pow(u2, 7)/std::pow(h, 8);
                q2 += (1.0/5040.0)*a_[37]*std::pow(u1, 7)/std::pow(h, 8) + (1.0/720.0)*a_[38]*std::pow(u1, 6)*u2/std::pow(h, 8) + (1.0/240.0)*a_[39]*std::pow(u1, 5)*std::pow(u2, 2)/std::pow(h, 8) + (1.0/144.0)*a_[40]*std::pow(u1, 4)*std::pow(u2, 3)/std::pow(h, 8) + (1.0/144.0)*a_[41]*std::pow(u1, 3)*std::pow(u2, 4)/std::pow(h, 8) + (1.0/240.0)*a_[42]*std::pow(u1, 2)*std::pow(u2, 5)/std::pow(h, 8) + (1.0/720.0)*a_[43]*u1*std::pow(u2, 6)/std::pow(h, 8) + (1.0/5040.0)*a_[44]*std::pow(u2, 7)/std::pow(h, 8);
            case 36:
                q1 += (1.0/720.0)*a_[28]*std::pow(u1, 6)/std::pow(h, 7) + (1.0/120.0)*a_[29]*std::pow(u1, 5)*u2/std::pow(h, 7) + (1.0/48.0)*a_[30]*std::pow(u1, 4)*std::pow(u2, 2)/std::pow(h, 7) + (1.0/36.0)*a_[31]*std::pow(u1, 3)*std::pow(u2, 3)/std::pow(h, 7) + (1.0/48.0)*a_[32]*std::pow(u1, 2)*std::pow(u2, 4)/std::pow(h, 7) + (1.0/120.0)*a_[33]*u1*std::pow(u2, 5)/std::pow(h, 7) + (1.0/720.0)*a_[34]*std::pow(u2, 6)/std::pow(h, 7);
                q2 += (1.0/720.0)*a_[29]*std::pow(u1, 6)/std::pow(h, 7) + (1.0/120.0)*a_[30]*std::pow(u1, 5)*u2/std::pow(h, 7) + (1.0/48.0)*a_[31]*std::pow(u1, 4)*std::pow(u2, 2)/std::pow(h, 7) + (1.0/36.0)*a_[32]*std::pow(u1, 3)*std::pow(u2, 3)/std::pow(h, 7) + (1.0/48.0)*a_[33]*std::pow(u1, 2)*std::pow(u2, 4)/std::pow(h, 7) + (1.0/120.0)*a_[34]*u1*std::pow(u2, 5)/std::pow(h, 7) + (1.0/720.0)*a_[35]*std::pow(u2, 6)/std::pow(h, 7);
            case 28:
                q1 += (1.0/120.0)*a_[21]*std::pow(u1, 5)/std::pow(h, 6) + (1.0/24.0)*a_[22]*std::pow(u1, 4)*u2/std::pow(h, 6) + (1.0/12.0)*a_[23]*std::pow(u1, 3)*std::pow(u2, 2)/std::pow(h, 6) + (1.0/12.0)*a_[24]*std::pow(u1, 2)*std::pow(u2, 3)/std::pow(h, 6) + (1.0/24.0)*a_[25]*u1*std::pow(u2, 4)/std::pow(h, 6) + (1.0/120.0)*a_[26]*std::pow(u2, 5)/std::pow(h, 6);
                q2 += (1.0/120.0)*a_[22]*std::pow(u1, 5)/std::pow(h, 6) + (1.0/24.0)*a_[23]*std::pow(u1, 4)*u2/std::pow(h, 6) + (1.0/12.0)*a_[24]*std::pow(u1, 3)*std::pow(u2, 2)/std::pow(h, 6) + (1.0/12.0)*a_[25]*std::pow(u1, 2)*std::pow(u2, 3)/std::pow(h, 6) + (1.0/24.0)*a_[26]*u1*std::pow(u2, 4)/std::pow(h, 6) + (1.0/120.0)*a_[27]*std::pow(u2, 5)/std::pow(h, 6);
            case 21:
                q1 += (1.0/24.0)*a_[15]*std::pow(u1, 4)/std::pow(h, 5) + (1.0/6.0)*a_[16]*std::pow(u1, 3)*u2/std::pow(h, 5) + (1.0/4.0)*a_[17]*std::pow(u1, 2)*std::pow(u2, 2)/std::pow(h, 5) + (1.0/6.0)*a_[18]*u1*std::pow(u2, 3)/std::pow(h, 5) + (1.0/24.0)*a_[19]*std::pow(u2, 4)/std::pow(h, 5);
                q2 += (1.0/24.0)*a_[16]*std::pow(u1, 4)/std::pow(h, 5) + (1.0/6.0)*a_[17]*std::pow(u1, 3)*u2/std::pow(h, 5) + (1.0/4.0)*a_[18]*std::pow(u1, 2)*std::pow(u2, 2)/std::pow(h, 5) + (1.0/6.0)*a_[19]*u1*std::pow(u2, 3)/std::pow(h, 5) + (1.0/24.0)*a_[20]*std::pow(u2, 4)/std::pow(h, 5);
            case 15:
                q1 += (1.0/6.0)*a_[10]*std::pow(u1, 3)/std::pow(h, 4) + (1.0/2.0)*a_[11]*std::pow(u1, 2)*u2/std::pow(h, 4) + (1.0/2.0)*a_[12]*u1*std::pow(u2, 2)/std::pow(h, 4) + (1.0/6.0)*a_[13]*std::pow(u2, 3)/std::pow(h, 4);
                q2 += (1.0/6.0)*a_[11]*std::pow(u1, 3)/std::pow(h, 4) + (1.0/2.0)*a_[12]*std::pow(u1, 2)*u2/std::pow(h, 4) + (1.0/2.0)*a_[13]*u1*std::pow(u2, 2)/std::pow(h, 4) + (1.0/6.0)*a_[14]*std::pow(u2, 3)/std::pow(h, 4);
            case 10:
                q1 += (1.0/2.0)*a_[6]*std::pow(u1, 2)/std::pow(h, 3) + a_[7]*u1*u2/std::pow(h, 3) + (1.0/2.0)*a_[8]*std::pow(u2, 2)/std::pow(h, 3);
                q2 += (1.0/2.0)*a_[7]*std::pow(u1, 2)/std::pow(h, 3) + a_[8]*u1*u2/std::pow(h, 3) + (1.0/2.0)*a_[9]*std::pow(u2, 2)/std::pow(h, 3);
            case 6:
                q1 += a_[3]*u1/std::pow(h, 2) + a_[4]*u2/std::pow(h, 2);
                q2 += a_[4]*u1/std::pow(h, 2) + a_[5]*u2/std::pow(h, 2);
            case 3:
                q1 += a_[1]/h;
                q2 += a_[2]/h;
            case 1:
                break;
            default:
                compadre_kernel_assert_release(false && "curvature polynomial order is greater than 8.");

        }
        q1 = std::pow(q1, 2);
        q2 = std::pow(q2, 2);
        det_g += q1;
        det_g += q2;
        return det_g;
    }

    //! Gaussian curvature K at any point in the local chart
    KOKKOS_INLINE_FUNCTION
    double GaussianCurvature(const scratch_vector_type a_, const double h, const double u1, const double u2) {

        double q1=0, q2=0, q3=0;
        switch (a_.extent(0)) {
            case 45:
                q1 += (1.0/720.0)*a_[36]*std::pow(u1, 6)/std::pow(h, 8) + (1.0/120.0)*a_[37]*std::pow(u1, 5)*u2/std::pow(h, 8) + (1.0/48.0)*a_[38]*std::pow(u1, 4)*std::pow(u2, 2)/std::pow(h, 8) + (1.0/36.0)*a_[39]*std::pow(u1, 3)*std::pow(u2, 3)/std::pow(h, 8) + (1.0/48.0)*a_[40]*std::pow(u1, 2)*std::pow(u2, 4)/std::pow(h, 8) + (1.0/120.0)*a_[41]*u1*std::pow(u2, 5)/std::pow(h, 8) + (1.0/720.0)*a_[42]*std::pow(u2, 6)/std::pow(h, 8);
                q2 += (1.0/720.0)*a_[38]*std::pow(u1, 6)/std::pow(h, 8) + (1.0/120.0)*a_[39]*std::pow(u1, 5)*u2/std::pow(h, 8) + (1.0/48.0)*a_[40]*std::pow(u1, 4)*std::pow(u2, 2)/std::pow(h, 8) + (1.0/36.0)*a_[41]*std::pow(u1, 3)*std::pow(u2, 3)/std::pow(h, 8) + (1.0/48.0)*a_[42]*std::pow(u1, 2)*std::pow(u2, 4)/std::pow(h, 8) + (1.0/120.0)*a_[43]*u1*std::pow(u2, 5)/std::pow(h, 8) + (1.0/720.0)*a_[44]*std::pow(u2, 6)/std::pow(h, 8);
                q3 += (1.0/720.0)*a_[37]*std::pow(u1, 6)/std::pow(h, 8) + (1.0/120.0)*a_[38]*std::pow(u1, 5)*u2/std::pow(h, 8) + (1.0/48.0)*a_[39]*std::pow(u1, 4)*std::pow(u2, 2)/std::pow(h, 8) + (1.0/36.0)*a_[40]*std::pow(u1, 3)*std::pow(u2, 3)/std::pow(h, 8) + (1.0/48.0)*a_[41]*std::pow(u1, 2)*std::pow(u2, 4)/std::pow(h, 8) + (1.0/120.0)*a_[42]*u1*std::pow(u2, 5)/std::pow(h, 8) + (1.0/720.0)*a_[43]*std::pow(u2, 6)/std::pow(h, 8);
            case 36:
                q1 += (1.0/120.0)*a_[28]*std::pow(u1, 5)/std::pow(h, 7) + (1.0/24.0)*a_[29]*std::pow(u1, 4)*u2/std::pow(h, 7) + (1.0/12.0)*a_[30]*std::pow(u1, 3)*std::pow(u2, 2)/std::pow(h, 7) + (1.0/12.0)*a_[31]*std::pow(u1, 2)*std::pow(u2, 3)/std::pow(h, 7) + (1.0/24.0)*a_[32]*u1*std::pow(u2, 4)/std::pow(h, 7) + (1.0/120.0)*a_[33]*std::pow(u2, 5)/std::pow(h, 7);
                q2 += (1.0/120.0)*a_[30]*std::pow(u1, 5)/std::pow(h, 7) + (1.0/24.0)*a_[31]*std::pow(u1, 4)*u2/std::pow(h, 7) + (1.0/12.0)*a_[32]*std::pow(u1, 3)*std::pow(u2, 2)/std::pow(h, 7) + (1.0/12.0)*a_[33]*std::pow(u1, 2)*std::pow(u2, 3)/std::pow(h, 7) + (1.0/24.0)*a_[34]*u1*std::pow(u2, 4)/std::pow(h, 7) + (1.0/120.0)*a_[35]*std::pow(u2, 5)/std::pow(h, 7);
                q3 += (1.0/120.0)*a_[29]*std::pow(u1, 5)/std::pow(h, 7) + (1.0/24.0)*a_[30]*std::pow(u1, 4)*u2/std::pow(h, 7) + (1.0/12.0)*a_[31]*std::pow(u1, 3)*std::pow(u2, 2)/std::pow(h, 7) + (1.0/12.0)*a_[32]*std::pow(u1, 2)*std::pow(u2, 3)/std::pow(h, 7) + (1.0/24.0)*a_[33]*u1*std::pow(u2, 4)/std::pow(h, 7) + (1.0/120.0)*a_[34]*std::pow(u2, 5)/std::pow(h, 7);
            case 28:
                q1 += (1.0/24.0)*a_[21]*std::pow(u1, 4)/std::pow(h, 6) + (1.0/6.0)*a_[22]*std::pow(u1, 3)*u2/std::pow(h, 6) + (1.0/4.0)*a_[23]*std::pow(u1, 2)*std::pow(u2, 2)/std::pow(h, 6) + (1.0/6.0)*a_[24]*u1*std::pow(u2, 3)/std::pow(h, 6) + (1.0/24.0)*a_[25]*std::pow(u2, 4)/std::pow(h, 6);
                q2 += (1.0/24.0)*a_[23]*std::pow(u1, 4)/std::pow(h, 6) + (1.0/6.0)*a_[24]*std::pow(u1, 3)*u2/std::pow(h, 6) + (1.0/4.0)*a_[25]*std::pow(u1, 2)*std::pow(u2, 2)/std::pow(h, 6) + (1.0/6.0)*a_[26]*u1*std::pow(u2, 3)/std::pow(h, 6) + (1.0/24.0)*a_[27]*std::pow(u2, 4)/std::pow(h, 6);
                q3 += (1.0/24.0)*a_[22]*std::pow(u1, 4)/std::pow(h, 6) + (1.0/6.0)*a_[23]*std::pow(u1, 3)*u2/std::pow(h, 6) + (1.0/4.0)*a_[24]*std::pow(u1, 2)*std::pow(u2, 2)/std::pow(h, 6) + (1.0/6.0)*a_[25]*u1*std::pow(u2, 3)/std::pow(h, 6) + (1.0/24.0)*a_[26]*std::pow(u2, 4)/std::pow(h, 6);
            case 21:
                q1 += (1.0/6.0)*a_[15]*std::pow(u1, 3)/std::pow(h, 5) + (1.0/2.0)*a_[16]*std::pow(u1, 2)*u2/std::pow(h, 5) + (1.0/2.0)*a_[17]*u1*std::pow(u2, 2)/std::pow(h, 5) + (1.0/6.0)*a_[18]*std::pow(u2, 3)/std::pow(h, 5);
                q2 += (1.0/6.0)*a_[17]*std::pow(u1, 3)/std::pow(h, 5) + (1.0/2.0)*a_[18]*std::pow(u1, 2)*u2/std::pow(h, 5) + (1.0/2.0)*a_[19]*u1*std::pow(u2, 2)/std::pow(h, 5) + (1.0/6.0)*a_[20]*std::pow(u2, 3)/std::pow(h, 5);
                q3 += (1.0/6.0)*a_[16]*std::pow(u1, 3)/std::pow(h, 5) + (1.0/2.0)*a_[17]*std::pow(u1, 2)*u2/std::pow(h, 5) + (1.0/2.0)*a_[18]*u1*std::pow(u2, 2)/std::pow(h, 5) + (1.0/6.0)*a_[19]*std::pow(u2, 3)/std::pow(h, 5);
            case 15:
                q1 += (1.0/2.0)*a_[10]*std::pow(u1, 2)/std::pow(h, 4) + a_[11]*u1*u2/std::pow(h, 4) + (1.0/2.0)*a_[12]*std::pow(u2, 2)/std::pow(h, 4);
                q2 += (1.0/2.0)*a_[12]*std::pow(u1, 2)/std::pow(h, 4) + a_[13]*u1*u2/std::pow(h, 4) + (1.0/2.0)*a_[14]*std::pow(u2, 2)/std::pow(h, 4);
                q3 += (1.0/2.0)*a_[11]*std::pow(u1, 2)/std::pow(h, 4) + a_[12]*u1*u2/std::pow(h, 4) + (1.0/2.0)*a_[13]*std::pow(u2, 2)/std::pow(h, 4);
            case 10:
                q1 += a_[6]*u1/std::pow(h, 3) + a_[7]*u2/std::pow(h, 3);
                q2 += a_[8]*u1/std::pow(h, 3) + a_[9]*u2/std::pow(h, 3);
                q3 += a_[7]*u1/std::pow(h, 3) + a_[8]*u2/std::pow(h, 3);
            case 6:
                q1 += a_[3]/std::pow(h, 2);
                q2 += a_[5]/std::pow(h, 2);
                q3 += a_[4]/std::pow(h, 2);
            // case 3: nothing to do here
            case 1:
                break;
            default:
                compadre_kernel_assert_release(false && "curvature polynomial order is greater than 8.");

        }
        const double det_g = MetricFactor(a_, h, u1, u2);
        return (q1*q2-std::pow(q3,2))/(det_g*det_g);

    }

    //! Surface curl at any point in the local chart
    KOKKOS_INLINE_FUNCTION
    double SurfaceCurlOfScalar(const scratch_vector_type a_, const double h, const double u1, const double u2, int x_pow, int y_pow, const int component) {

        const double factorial[15] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};

        const double det_g = MetricFactor(a_, h, u1, u2);

        int partial_direction = -1;
        if (component == 0) { // comp 0 of surface curl, i.e. dy(p)
            partial_direction = 1;
        } else if (component == 1) { // comp 1 of surface curl, i.e. -dx(p)
            partial_direction = 0;
        } else {
            compadre_kernel_assert_release(false);
        }
        const double sign_change = (component == 1) ? -1 : 1;
        int n_x_pow = (partial_direction == 0) ? x_pow-1 : x_pow;
        int n_y_pow = (partial_direction == 1) ? y_pow-1 : y_pow;

        double return_val = 0;
        if (n_x_pow<0 || n_y_pow<0) {
            return 0;
        } else {
            double alphaf = factorial[n_x_pow]*factorial[n_y_pow];
            return_val = 1./h 
                        *std::pow(u1/h,n_x_pow)
                        *std::pow(u2/h,n_y_pow)/alphaf;
        }
        return_val *= sign_change;
        return_val /= sqrt(det_g);
        return return_val;

    }

} // Compadre namespace
