// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_TRANSFER_UTIL_LEAST_SQUARES_HPP
#define STK_TRANSFER_UTIL_LEAST_SQUARES_HPP

#include <vector>
#include "stk_search/DistanceComparison.hpp" // for stk::search::less_than

namespace stk {
namespace mesh {
class BulkData;
class FieldBase;
struct Entity;
} // namespace mesh
} // namespace stk

namespace stk {
namespace transfer {

constexpr unsigned col_major_index(unsigned i, unsigned j, unsigned nx, unsigned ny) { return (i + j * nx); }

static constexpr int INT_ONE = 1;
static constexpr double REAL_ONE = 1.0;
static constexpr double ZERO = 0.0;
static constexpr char NOTRANSPOSE = 'N';
static constexpr char TRANSPOSE = 'T';

enum PatchRecoveryEvaluationType {
  LINEAR_LEAST_SQUARES = 0,
  LINEAR_MOVING_LEAST_SQUARES,
  QUADRATIC_LEAST_SQUARES,
  QUADRATIC_MOVING_LEAST_SQUARES,
  UNDEFINED_PATCH_RECOVERY_EVALUATION_TYPE = 0xff
};

class LeastSquares {
 public:
  LeastSquares(const unsigned ncomp, const unsigned nsamp, const unsigned nbasis)
    : m_numComponents(ncomp)
    , m_numSamples(nsamp)
    , m_numBasis(nbasis)
    , m_numData(ncomp * nsamp)
  {
    resize_data(ncomp, nsamp, nbasis);
  }

  virtual ~LeastSquares(){};

  virtual int least_squares(const unsigned nrhs, const std::vector<double>& fieldVal,
                            const std::vector<double>& basisVal, double* coeff);

  virtual int least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal,
                                 const std::vector<double>& basisVal, const double rcond, double* coeff);

  std::vector<double>& get_double_scratch_space() { return m_doubleScratchSpace; }

  void resize_data(const unsigned ncomp, const unsigned nsamp, const unsigned nbasis);

  unsigned get_num_components() const { return m_numComponents; }
  unsigned get_num_samples() const { return m_numSamples; }
  unsigned get_num_basis() const { return m_numBasis; }

 protected:
  unsigned m_numComponents;
  unsigned m_numSamples;
  unsigned m_numBasis;
  unsigned m_numData;

  std::vector<int> m_integerScratchSpace;
  std::vector<double> m_doubleScratchSpace;

  std::vector<int> m_integerCondScratchSpace;
  std::vector<double> m_doubleCondScratchSpace;

  std::vector<double> m_tempVector;

  void fill_covariance_matrix(const std::vector<double>& A);

  int compute_least_squares_solution(const unsigned nrhs, const std::vector<double>& basisVal,
                                     const std::vector<double>& fieldVal, double* coeff);

  int compute_least_squares_solution_with_factorization(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                        const std::vector<double>& basisVal, double* coeff);

  double compute_one_norm();

  int compute_factorization_and_condition_number(double& cond);
};

class GeometricMovingLeastSquares : public LeastSquares {
 public:
  GeometricMovingLeastSquares(const unsigned ncomp, const unsigned nsamp, const unsigned nbasis,
                              const stk::mesh::BulkData& bulk, const std::vector<stk::mesh::Entity>& patch,
                              const stk::mesh::FieldBase* coordField, const std::vector<double>& evalPoint);

  int least_squares(const unsigned nrhs, const std::vector<double>& fieldVal, const std::vector<double>& basisVal,
                    double* coeff) override;

  int least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal, const std::vector<double>& basisVal,
                         const double rcond, double* coeff) override;

  std::vector<double> compute_weights(double EPS = stk::search::STK_COMPARISON_EPSILON);

 private:
  const std::vector<stk::mesh::Entity>& m_patch;
  const stk::mesh::FieldBase* m_coordField;
  const std::vector<double>& m_evalPoint;

  std::vector<double> m_weights;

  void fill_diagonal_scaled_vector(const std::vector<double>& D, const double* f, std::vector<double>& x);

  void fill_covariance_matrix(const std::vector<double>& A, const std::vector<double>& W);

  int moving_least_squares(const unsigned nrhs, const std::vector<double>& fieldVal,
                           const std::vector<double>& basisVal, const std::vector<double>& weight, double* coeff);

  int moving_least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal,
                                const std::vector<double>& basisVal, const std::vector<double>& weight,
                                const double rcond, double* coeff);

  int compute_moving_least_squares_solution(const unsigned nrhs, const std::vector<double>& fieldVal,
                                            const std::vector<double>& basisVal, const std::vector<double>& weight,
                                            double* coeff);

  int compute_moving_least_squares_solution_with_factorization(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                               const std::vector<double>& basisVal,
                                                               const std::vector<double>& weight, double* coeff);
};

} // namespace transfer
} // namespace stk

#endif /* STK_TRANSFER_UTIL_LEAST_SQUARES_HPP */
