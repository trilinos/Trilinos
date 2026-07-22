// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/*! \file monitor.h
 *  \brief Monitor for iterative solver convergence with multiple right hand sides
 */

#pragma once

#include <cusp/detail/config.h>
#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/array2d.h>
#include <cusp/print.h>

#include <limits>
#include <iostream>
#include <iomanip>

// Classes to monitor iterative solver progress, check for convergence, etc.
// Follows the implementation of Iteration in the ITL:
//   http://www.osl.iu.edu/research/itl/doc/Iteration.html

namespace cusp
{
/*! \addtogroup iterative_solvers Iterative Solvers
 *  \addtogroup monitors Monitors
 *  \ingroup iterative_solvers
 *  \{
 */

/*! \p default_monitor : Implements standard convergence criteria
 * and reporting for iterative solvers.
 *
 * \tparam ValueType scalar type used in the solver (e.g. \c float or \c cusp::complex<double>).
 *
 *  \see \p verbose_monitor
 *
 */
template <typename ValueType>
class default_block_monitor
{
public:
  typedef typename norm_type<ValueType>::type Real;

  /*! Construct a \p default_monitor for a given right-hand-side \p b
   *
   *  The \p default_monitor terminates iteration when the residual norm
   *  satisfies the condition
   *       ||b - A x|| <= absolute_tolerance + relative_tolerance * ||b||
   *  or when the iteration limit is reached.
   *
   *  \param b right-hand-side of the linear system A x = b
   *  \param iteration_limit maximum number of solver iterations to allow
   *  \param relative_tolerance determines convergence criteria
   *  \param absolute_tolerance determines convergence criteria
   *
   *  \tparam VectorType vector
   */
  template <typename MV>
  default_block_monitor(const MV& b,
                        size_t iteration_limit = 500,
                        Real absolute_tolerance = 1e-6,
                        Real relative_tolerance = 1e-6,
                        bool verbose = true) :
    numRHS(b.num_cols),
    iteration_limit_(iteration_limit),
    iteration_count_(0),
    relative_tolerance_(relative_tolerance),
    absolute_tolerance_(absolute_tolerance),
    verbose_(verbose),
    b_norm(b.num_cols)
    {
      for (int i = 0; i < numRHS; i++)
        b_norm[i] = cusp::blas::nrm2(b.column(i));
    }

  /*! increment the iteration count
   */
  void operator++(void) {  ++iteration_count_; } // prefix increment

  /*! applies convergence criteria to determine whether iteration is finished
   *
   *  \param r residual vector of the linear system (r = b - A x)
   *  \tparam Vector vector
   */
  template <typename MV>
  bool finished(const MV& r)
  {

    if (converged(r))
    {
      if (verbose_) {
        cusp::array1d<ValueType, cusp::host_memory> resid(numRHS);
        std::cout << "Successfully converged after " << iteration_count() << " iterations to tolerance " << tolerance(0) << std::endl;
        std::cout << "with max residual norm ";
        Real norm_max = 0;
        for (int i = 0; i < numRHS; i++) {
          resid[i] = cusp::blas::nrm2(r.column(i));
          if (resid[i] >  norm_max) norm_max = resid[i];
        }
        std::cout << norm_max << std::endl;

        //cusp::print(resid);
      }

      return true;
    }
    else if (iteration_count() >= iteration_limit())
    {
      if (verbose_) {
        cusp::array1d<ValueType, cusp::host_memory> resid(numRHS);
        std::cout << "Failed to converge after " << iteration_count() << " iterations." << std::endl;
        std::cout << "with max residual norm ";
        Real norm_max = 0;
        for (int i = 0; i < numRHS; i++) {
          resid[i] = cusp::blas::nrm2(r.column(i));
          if (resid[i] >  norm_max) norm_max = resid[i];
        }
        std::cout << norm_max << std::endl;

        //cusp::print(resid);
      }

      return true;
    }
    else
    {
      return false;
    }


  }

  /*! whether the last tested residual satifies the convergence tolerance
   */
  template <typename MV>
  bool converged(MV& r) const
  {
    for (int i = 0; i < numRHS; i++){

      if (cusp::blas::nrm2(r.column(i)) > tolerance(i)){
        return false;
      }
    }

    return true;
  }

  /*! number of iterations
   */
  size_t iteration_count() const { return iteration_count_; }

  /*! maximum number of iterations
   */
  size_t iteration_limit() const { return iteration_limit_; }

  /*! relative tolerance
   */
  Real relative_tolerance() const { return relative_tolerance_; }

  /*! absolute tolerance
   */
  Real absolute_tolerance() const { return absolute_tolerance_; }

  /*! tolerance
   *
   *  Equal to absolute_tolerance() + relative_tolerance() * ||b||
   *
   */
  Real tolerance(int i) const { return absolute_tolerance() + relative_tolerance() * b_norm[i]; }

protected:

  Real relative_tolerance_;
  Real absolute_tolerance_;
  bool verbose_;
  size_t numRHS;
  size_t iteration_limit_;
  size_t iteration_count_;
  cusp::array1d<ValueType, cusp::host_memory> b_norm;

};

/*! \p verbose_monitor is similar to \p default monitor except that
 * it displays the solver status during iteration and reports a
 * summary after iteration has stopped.
 *
 * \tparam ValueType scalar type used in the solver (e.g. \c float or \c cusp::complex<double>).
 *
 * \see \p default_monitor
 */

/*
template <typename ValueType>
class verbose_block_monitor : public default_block_monitor<ValueType>
{
    typedef typename norm_type<ValueType>::type Real;
    typedef cusp::default_block_monitor<ValueType> super;

    public:

    *! Construct a \p verbose_monitor for a given right-hand-side \p b
     *
     *  The \p verbose_monitor terminates iteration when the residual norm
     *  satisfies the condition
     *       ||b - A x|| <= absolute_tolerance + relative_tolerance * ||b||
     *  or when the iteration limit is reached.
     *
     *  \param b right-hand-side of the linear system A x = b
     *  \param iteration_limit maximum number of solver iterations to allow
     *  \param relative_tolerance determines convergence criteria
     *  \param absolute_tolerance determines convergence criteria
     *
     *  \tparam VectorType vector
     *


    template <typename MV>
    verbose_block_monitor(const MV& b, size_t iteration_limit = 500, Real relative_tolerance = 1e-5, Real absolute_tolerance = 0)
        : super(b, iteration_limit, relative_tolerance, absolute_tolerance)
    {
        std::cout << "Solver will continue until ";
        std::cout << "residual norm" << super::tolerance() << " or reaching ";
        std::cout << super::iteration_limit() << " iterations " << std::endl;
        std::cout << "  Iteration Number  | Residual Norm" << std::endl;
    }

    template <typename MV>
    bool finished(const MV& r)
    {
        for (int i = 0; i < r.num_cols; i++)
                super::r_norm[i] = cusp::blas::nrm2(r.column(i));

        std::cout << "       "  << std::setw(10) << super::iteration_count();
        std::cout << "       "  << std::setw(10) << std::scientific << super::residual_norm_average() << std::endl;

        if (super::converged(r))
        {
            std::cout << "Successfully converged after " << super::iteration_count() << " iterations." << std::endl;
            return true;
        }
        else if (super::iteration_count() >= super::iteration_limit())
        {
            std::cout << "Failed to converge after " << super::iteration_count() << " iterations." << std::endl;
            return true;
        }
        else
        {
            return false;
        }
    }
};

*/

/*! \}
 */

} // end namespace cusp
