// @HEADER
// *****************************************************************************
//                   KLU2: A Direct Linear Solver package
//
// Copyright 2011 NTESS and the KLU2 contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <complex>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "klu2_defaults.hpp"
#include "klu2_analyze.hpp"
#include "klu2_factor.hpp"
#include "klu2_solve.hpp"
#include "klu2_tsolve.hpp"
#include "klu2_free_symbolic.hpp"
#include "klu2_free_numeric.hpp"

/**
 * @file
 *
 * This file contains a standalone test for the @c KLU2 solver.
 */

/**
 * @name Data.
 *
 * This simple example will solve the following linear system:
 *
 * @verbatim
 *
 *          Matrix               Solution   Right-hand side
 *
 * | 2. |  3. |     |    |    |    [ 1 ]       [  8 ]
 * | 3. |     |  4. |    | 6. |    [ 2 ]       [ 45 ]
 * |    | -1. | -3. | 2. |    |    [ 3 ]     = [ -3 ]
 * |    |     |  1. |    |    |    [ 4 ]       [  3 ]
 * |    |  4. |  2. |    | 1. |    [ 5 ]       [ 19 ]
 *
 * @endverbatim
 *
 * @warning The algorithm works with the compressed column storage format (CCS).
 */
///@{
//! Size of the system.
const int size = 5;
//! Column offsets.
int Ap [ ] = {0, 2, 5, 9, 10, 12};
//! Row values.
int Ai [ ] = { 0 , 1 , 0 ,  2 , 4 , 1 ,  2 , 3 , 4 , 2 , 1 , 4 };

//! Get entries values.
template <typename T>
auto get_Ax() {
    return std::vector<T>{ 2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.};
}

//! Get right-hand side values.
template <typename T>
auto get_b() {
    return std::vector<T>{8., 45., -3., 3., 19.};
}

//! Get solution.
template <typename T>
auto get_x() {
    return std::vector<T>{1., 2., 3., 4., 5.};
}

//! Get solution for transposed matrix.
template <typename T>
auto get_x_transpose() {
    return std::vector<T>{69. / 38., 83./57., 3./2., -2833/114., 195./19.};
}
///@}

//! Check that @p a and @p b are equal.
#define ASSERT_EQ(a, b) if(a != b) throw std::runtime_error(#a " is different from " #b);

//! Check that @p a and @p b are equal to the given asolute tolerance.
#define ASSERT_TOL_EQ(a, b, tol) if( std::abs(a-b) > tol ) throw std::runtime_error(#a " is different from " #b);

//! Print vectors and compare them.
template <typename T, typename tol_t>
void print_and_compare(const std::vector<T>& values_a, const std::vector<T>& values_b, const tol_t tol)
{
    std::cout << "Comparing vector:" << std::endl << "\t" << std::scientific;
    for(size_t i = 0; i < values_a.size(); ++i) std::cout << values_a.at(i) << ",";
    std::cout << std::endl << "with vector" << std::endl << "\t";
    for(size_t i = 0; i < values_b.size(); ++i) std::cout << values_b.at(i) << ",";
    std::cout << std::endl;

    ASSERT_EQ(values_a.size(), values_b.size());

    for(size_t i = 0; i < values_a.size(); ++i) ASSERT_TOL_EQ(values_a.at(i), values_b.at(i), tol);
}

/**
 * Solve and check the solution.
 */
template <typename T, typename D = int, typename tol_t>
void solve_and_check(const tol_t tol)
{
    klu_symbolic<T, D> *Symbolic;
    klu_numeric <T, D> *Numeric;
    klu_common  <T, D>  Common;

    klu_defaults<T, D> (&Common);

    //! Query data for @tparam T.
          auto Ax   = get_Ax         <T>();
    const auto x    = get_x          <T>();
    const auto x_tr = get_x_transpose<T>();

    //! Symbolic step.
    Symbolic = klu_analyze<T, D>(size, Ap, Ai, &Common);

    //! Numeric step.
    Numeric = klu_factor<T, D>(Ap, Ai, Ax.data(), Symbolic, &Common);

    //! Solve and check against expected solution. Note that @ref klu_solve writes the solution into @c b.
    {
        auto b = get_b<T>();
        klu_solve<T, D>(Symbolic, Numeric, size, 1, b.data(), &Common);
        print_and_compare<T>(x, b, tol);
    }

    //! Transpose solve. Note that @ref klu_solve writes the solution into @c b.
    {
        auto b = get_b<T>();
        klu_tsolve<T, D>(Symbolic, Numeric, size, 1, b.data(), &Common);
        print_and_compare<T>(x_tr, b, tol);
    }

    //! Free memory.
    klu_free_symbolic<T, D>(&Symbolic, &Common);
    klu_free_numeric <T, D>(&Numeric, &Common);
}

int main (void)
{
    //! Solve for @c float.
    solve_and_check<float>(float (1.e-6));

    //! Solve for @c double.
    solve_and_check<double>(double(1.e-15));

    //! Solve for @c std::complex<float>. @todo This does not compile for now!
    // solve_and_check<std::complex<float>>(float(1.e-6));

#ifdef HAVE_TEUCHOS_COMPLEX
    //! Solve for @c std::complex<double>.
    solve_and_check<std::complex<double>>(double(1.e-15));
#endif

    return EXIT_SUCCESS;
}
