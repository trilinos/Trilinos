// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_SCALAR_TAYLORPOLYNOMIAL_BASIS_HPP_
#define _COMPADRE_SCALAR_TAYLORPOLYNOMIAL_BASIS_HPP_

#include "Compadre_GMLS.hpp"

namespace Compadre {
/*!  \brief Definition of scalar Taylor polynomial basis
   
   For order P, the sum of the basis is defined by:
   \li in 1D) \f[\sum_{n=0}^{n=P} \frac{(x/h)^n}{n!}\f]
   \li in 2D) \f[\sum_{n=0}^{n=P} \sum_{k=0}^{k=n} \frac{(x/h)^{n-k}*(y/h)^k}{(n-k)!k!}\f]
   \li in 3D) \f[\sum_{p=0}^{p=P} \sum_{k_1+k_2+k_3=n} \frac{(x/h)^{k_1}*(y/h)^{k_2}*(z/h)^{k_3}}{k_1!k_2!k_3!}\f]
*/
namespace ScalarTaylorPolynomialBasis {
   
    /*! \brief Returns size of basis
        \param degree               [in] - highest degree of polynomial
        \param dimension            [in] - spatial dimension to evaluate
    */
    KOKKOS_INLINE_FUNCTION
    int getSize(const int degree, const int dimension) {
        if (dimension == 3) return (degree+1)*(degree+2)*(degree+3)/6;
        else if (dimension == 2) return (degree+1)*(degree+2)/2;
        else return degree+1;
    }

    /*! \brief Evaluates the scalar Taylor polynomial basis
     *  delta[j] = weight_of_original_value * delta[j] + weight_of_new_value * (calculation of this function)
        \param delta                [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _basis_multipler*the dimension of the polynomial basis.
        \param workspace            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param dimension                [in] - spatial dimension to evaluate
        \param max_degree               [in] - highest degree of polynomial
        \param h                        [in] - epsilon/window size
        \param x                        [in] - x coordinate (already shifted by target)
        \param y                        [in] - y coordinate (already shifted by target)
        \param z                        [in] - z coordinate (already shifted by target)
        \param starting_order           [in] - polynomial order to start with (default=0)
        \param weight_of_original_value [in] - weighting to assign to original value in delta (default=0, replace, 1-increment current value)
        \param weight_of_new_value      [in] - weighting to assign to new contribution        (default=1, add to/replace)
    */
    KOKKOS_INLINE_FUNCTION
    void evaluate(const member_type& teamMember, double* delta, double* workspace, const int dimension, const int max_degree, const double h, const double x, const double y, const double z, const int starting_order = 0, const double weight_of_original_value = 0.0, const double weight_of_new_value = 1.0) {
        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
            const double factorial[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
            if (dimension==3) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                scratch_vector_type z_over_h_to_i(workspace+2*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                z_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                    z_over_h_to_i[i] = z_over_h_to_i[i-1]*(z/h);
                }
                // (in 3D) \sum_{p=0}^{p=P} \sum_{k1+k2+k3=n} (x/h)^k1*(y/h)^k2*(z/h)^k3 / (k1!k2!k3!)
                int alphax, alphay, alphaz;
                double alphaf;
                int i=0, s=0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphaz = 0; alphaz <= n; alphaz++){
                        s = n - alphaz;
                        for (alphay = 0; alphay <= s; alphay++){
                            alphax = s - alphay;
                            alphaf = factorial[alphax]*factorial[alphay]*factorial[alphaz];
                            *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * x_over_h_to_i[alphax]*y_over_h_to_i[alphay]*z_over_h_to_i[alphaz]/alphaf;
                            i++;
                        }
                    }
                }
            } else if (dimension==2) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                }
                // (in 2D) \sum_{n=0}^{n=P} \sum_{k=0}^{k=n} (x/h)^(n-k)*(y/h)^k / ((n-k)!k!)
                int alphax, alphay;
                double alphaf;
                int i = 0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        alphax = n - alphay;
                        alphaf = factorial[alphax]*factorial[alphay];
                        *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * x_over_h_to_i[alphax]*y_over_h_to_i[alphay]/alphaf;
                        i++;
                    }
                }
            } else {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                x_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                }
                // (in 1D) \sum_{n=0}^{n=P} (x/h)^n / n!
                for (int i=starting_order; i<=max_degree; ++i) {
                    *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * x_over_h_to_i[i]/factorial[i];
                }
            }
        });
    }

    /*! \brief Evaluates the first partial derivatives of scalar Taylor polynomial basis
     *  delta[j] = weight_of_original_value * delta[j] + weight_of_new_value * (calculation of this function)
        \param delta                [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _basis_multipler*the dimension of the polynomial basis.
        \param workspace            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param dimension                [in] - spatial dimension to evaluate
        \param max_degree               [in] - highest degree of polynomial
        \param partial_direction        [in] - direction in which to take partial derivative
        \param h                        [in] - epsilon/window size
        \param x                        [in] - x coordinate (already shifted by target)
        \param y                        [in] - y coordinate (already shifted by target)
        \param z                        [in] - z coordinate (already shifted by target)
        \param starting_order           [in] - polynomial order to start with (default=0)
        \param weight_of_original_value [in] - weighting to assign to original value in delta (default=0, replace, 1-increment current value)
        \param weight_of_new_value      [in] - weighting to assign to new contribution        (default=1, add to/replace)
    */
    KOKKOS_INLINE_FUNCTION
    void evaluatePartialDerivative(const member_type& teamMember, double* delta, double* workspace, const int dimension, const int max_degree, const int partial_direction, const double h, const double x, const double y, const double z, const int starting_order = 0, const double weight_of_original_value = 0.0, const double weight_of_new_value = 1.0) {
        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
            const double factorial[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
            if (dimension==3) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                scratch_vector_type z_over_h_to_i(workspace+2*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                z_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                    z_over_h_to_i[i] = z_over_h_to_i[i-1]*(z/h);
                }
                // (in 3D) \sum_{p=0}^{p=P} \sum_{k1+k2+k3=n} (x/h)^k1*(y/h)^k2*(z/h)^k3 / (k1!k2!k3!)
                int alphax, alphay, alphaz;
                double alphaf;
                int i=0, s=0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphaz = 0; alphaz <= n; alphaz++){
                        s = n - alphaz;
                        for (alphay = 0; alphay <= s; alphay++){
                            alphax = s - alphay;

                            int var_pow[3] = {alphax, alphay, alphaz};
                            var_pow[partial_direction]--;

                            if (var_pow[0]<0 || var_pow[1]<0 || var_pow[2]<0) {
                                *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                            } else {
                                alphaf = factorial[var_pow[0]]*factorial[var_pow[1]]*factorial[var_pow[2]];
                                *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]*y_over_h_to_i[var_pow[1]]*z_over_h_to_i[var_pow[2]]/alphaf;
                            }
                            i++;
                        }
                    }
                }
            } else if (dimension==2) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                }
                // (in 3D) \sum_{p=0}^{p=P} \sum_{k1+k2+k3=n} (x/h)^k1*(y/h)^k2*(z/h)^k3 / (k1!k2!k3!)
                int alphax, alphay;
                double alphaf;
                int i=0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        alphax = n - alphay;

                        int var_pow[2] = {alphax, alphay};
                        var_pow[partial_direction]--;

                        if (var_pow[0]<0 || var_pow[1]<0) {
                            *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                        } else {
                            alphaf = factorial[var_pow[0]]*factorial[var_pow[1]];
                            *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]*y_over_h_to_i[var_pow[1]]/alphaf;
                        }
                        i++;
                    }
                }
            } else {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                x_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                }
                int alphax;
                double alphaf;
                int i=0;
                for (int n = starting_order; n <= max_degree; n++){
                    alphax = n;
                
                    int var_pow[1] = {alphax};
                    var_pow[partial_direction]--;

                    if (var_pow[0]<0) {
                        *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                    } else {
                        alphaf = factorial[var_pow[0]];
                        *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]/alphaf;
                    }
                    i++;
                }
            }
        });
    }

    /*! \brief Evaluates the second partial derivatives of scalar Taylor polynomial basis
     *  delta[j] = weight_of_original_value * delta[j] + weight_of_new_value * (calculation of this function)
        \param delta                [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _basis_multipler*the dimension of the polynomial basis.
        \param workspace            [in/out] - scratch space that is allocated so that each thread has its own copy. Must be at least as large as the _poly_order*the spatial dimension of the polynomial basis.
        \param dimension                [in] - spatial dimension to evaluate
        \param max_degree               [in] - highest degree of polynomial
        \param partial_direction_1      [in] - direction in which to take first partial derivative
        \param partial_direction_2      [in] - direction in which to take second partial derivative
        \param h                        [in] - epsilon/window size
        \param x                        [in] - x coordinate (already shifted by target)
        \param y                        [in] - y coordinate (already shifted by target)
        \param z                        [in] - z coordinate (already shifted by target)
        \param starting_order           [in] - polynomial order to start with (default=0)
        \param weight_of_original_value [in] - weighting to assign to original value in delta (default=0, replace, 1-increment current value)
        \param weight_of_new_value      [in] - weighting to assign to new contribution        (default=1, add to/replace)
    */
    KOKKOS_INLINE_FUNCTION
    void evaluateSecondPartialDerivative(const member_type& teamMember, double* delta, double* workspace, const int dimension, const int max_degree, const int partial_direction_1, const int partial_direction_2, const double h, const double x, const double y, const double z, const int starting_order = 0, const double weight_of_original_value = 0.0, const double weight_of_new_value = 1.0) {
        Kokkos::single(Kokkos::PerThread(teamMember), [&] () {
            const double factorial[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200};
            if (dimension==3) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                scratch_vector_type z_over_h_to_i(workspace+2*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                z_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                    z_over_h_to_i[i] = z_over_h_to_i[i-1]*(z/h);
                }
                // (in 3D) \sum_{p=0}^{p=P} \sum_{k1+k2+k3=n} (x/h)^k1*(y/h)^k2*(z/h)^k3 / (k1!k2!k3!)
                int alphax, alphay, alphaz;
                double alphaf;
                int i=0, s=0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphaz = 0; alphaz <= n; alphaz++){
                        s = n - alphaz;
                        for (alphay = 0; alphay <= s; alphay++){
                            alphax = s - alphay;

                            int var_pow[3] = {alphax, alphay, alphaz};
                            var_pow[partial_direction_1]--;
                            var_pow[partial_direction_2]--;

                            if (var_pow[0]<0 || var_pow[1]<0 || var_pow[2]<0) {
                                *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                            } else {
                                alphaf = factorial[var_pow[0]]*factorial[var_pow[1]]*factorial[var_pow[2]];
                                *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]*y_over_h_to_i[var_pow[1]]*z_over_h_to_i[var_pow[2]]/alphaf;
                            }
                            i++;
                        }
                    }
                }
            } else if (dimension==2) {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                scratch_vector_type y_over_h_to_i(workspace+1*(max_degree+1), max_degree+1);
                x_over_h_to_i[0] = 1;
                y_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                    y_over_h_to_i[i] = y_over_h_to_i[i-1]*(y/h);
                }
                // (in 3D) \sum_{p=0}^{p=P} \sum_{k1+k2+k3=n} (x/h)^k1*(y/h)^k2*(z/h)^k3 / (k1!k2!k3!)
                int alphax, alphay;
                double alphaf;
                int i=0;
                for (int n = starting_order; n <= max_degree; n++){
                    for (alphay = 0; alphay <= n; alphay++){
                        alphax = n - alphay;

                        int var_pow[2] = {alphax, alphay};
                        var_pow[partial_direction_1]--;
                        var_pow[partial_direction_2]--;

                        if (var_pow[0]<0 || var_pow[1]<0) {
                            *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                        } else {
                            alphaf = factorial[var_pow[0]]*factorial[var_pow[1]];
                            *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]*y_over_h_to_i[var_pow[1]]/alphaf;
                        }
                        i++;
                    }
                }
            } else {
                scratch_vector_type x_over_h_to_i(workspace, max_degree+1);
                x_over_h_to_i[0] = 1;
                for (int i=1; i<=max_degree; ++i) {
                    x_over_h_to_i[i] = x_over_h_to_i[i-1]*(x/h);
                }
                int alphax;
                double alphaf;
                int i=0;
                for (int n = starting_order; n <= max_degree; n++){
                    alphax = n;
                
                    int var_pow[1] = {alphax};
                    var_pow[partial_direction_1]--;
                    var_pow[partial_direction_2]--;

                    if (var_pow[0]<0) {
                        *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 0.0;
                    } else {
                        alphaf = factorial[var_pow[0]];
                        *(delta+i) = weight_of_original_value * *(delta+i) + weight_of_new_value * 1./h * x_over_h_to_i[var_pow[0]]/alphaf;
                    }
                    i++;
                }
            }
        });
    }

} // ScalarTaylorPolynomialBasis

} // Compadre

#endif
