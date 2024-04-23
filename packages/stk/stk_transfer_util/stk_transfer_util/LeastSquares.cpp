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

#include "stk_util/environment/Env.hpp"            // for parallel_rank
#include "stk_util/environment/RuntimeMessage.hpp" // for MessageCode
#include "stk_util/environment/RuntimeWarning.hpp"
#include "stk_util/util/ReportHandler.hpp" // for STK_ThrowRequireMsg, etc
#include "stk_transfer_util/LeastSquares.hpp"
#include "stk_search_util/ObjectCoordinates.hpp"

#include "stk_util/util/BlasLapack.hpp"

namespace stk {
namespace transfer {

int LeastSquares::least_squares(const unsigned nrhs, const std::vector<double>& fieldVal,
                                const std::vector<double>& basisVal, double* coeff)
{
  STK_ThrowRequireMsg(fieldVal.size() >= m_numData * nrhs, "Insufficient size for fieldVal");
  STK_ThrowRequireMsg(basisVal.size() >= m_numData * m_numBasis, "Insufficient size for basisVal");

  fill_covariance_matrix(basisVal);
  return compute_least_squares_solution_with_factorization(nrhs, fieldVal, basisVal, coeff);
}

int LeastSquares::least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal,
                                     const std::vector<double>& basisVal, const double rcond, double* coeff)
{
  STK_ThrowRequireMsg(fieldVal.size() >= m_numData * nrhs, "Insufficient size for fieldVal");
  STK_ThrowRequireMsg(basisVal.size() >= m_numData * m_numBasis, "Insufficient size for basisVal");

  fill_covariance_matrix(basisVal);
  double cond;
  int info = compute_factorization_and_condition_number(cond);
  if(info != 0) return info;

  if(cond > rcond) {
    // Condition number is okay.
    // Can procede with LU
    info = compute_least_squares_solution(nrhs, fieldVal, basisVal, coeff);
  }
  else {
    info = -1;
  }

  return info;
}


int LeastSquares::compute_factorization_and_condition_number(double& cond)
{
  int info = 0;

  // Get the 1-norm of doubleScratchSpace
  double wNorm = compute_one_norm();

  // Replace covariance matrix with its LU factorization
  // using partial pivoting with row interchanges [LAPACK].
  int dim1 = m_numBasis;
  SIERRA_FORTRAN(dgetrf)(&dim1, &dim1, m_doubleScratchSpace.data(), &dim1, m_integerScratchSpace.data(), &info);

  // Check for error code.
  // If the matrix is exactly singular or there are an insufficient
  // number of sampling points, return.
  if(m_numData < m_numBasis) info = -2;
  if(info != 0) return info;

  // Get the condition number of doubleScratchSpace
  const char norm = '1';
  SIERRA_FORTRAN(dgecon)(&norm, &dim1, m_doubleScratchSpace.data(), &dim1, &wNorm, &cond, m_doubleCondScratchSpace.data(),
                         m_integerCondScratchSpace.data(), &info);
  return info;
}

double LeastSquares::compute_one_norm()
{
  // Get the 1-norm of doubleScratchSpace
  double wNorm = 0.0;
  double testNorm = 0.0;
  for(unsigned j = 0; j < m_numBasis; j++) {
    testNorm = 0.0;
    for(unsigned i = 0; i < m_numBasis; i++) {
      unsigned workIndex = col_major_index(i, j, m_numBasis, m_numBasis);
      testNorm += std::fabs(m_doubleScratchSpace[workIndex]);
    }
    wNorm = std::max(wNorm, testNorm);
  }
  return wNorm;
}

void LeastSquares::fill_covariance_matrix(const std::vector<double>& A)
{
  // Get the covariance matrix = transpose(A) * A
  // where matrix "A" has dimension(m_numData,m_numBasis).
  // Result is stored in m_doubleScratchSpace
  int dim1 = m_numBasis;
  int dim2 = m_numData;
  SIERRA_FORTRAN(dgemm)(&TRANSPOSE, &NOTRANSPOSE, &dim1, &dim1, &dim2, &REAL_ONE, A.data(), &dim2,
                        A.data(), &dim2, &ZERO, m_doubleScratchSpace.data(), &dim1);
}

int LeastSquares::compute_least_squares_solution_with_factorization(const unsigned nrhs,
                                                                    const std::vector<double>& fieldVal,
                                                                    const std::vector<double>& basisVal, double* coeff)
{
  // Replace covariance matrix with its LU factorization
  // using partial pivoting with row interchanges [LAPACK].
  int dim1 = m_numBasis;
  int info = 0;
  SIERRA_FORTRAN(dgetrf)(&dim1, &dim1, m_doubleScratchSpace.data(), &dim1, m_integerScratchSpace.data(), &info);

  if(info != 0) return info;

  return compute_least_squares_solution(nrhs, fieldVal, basisVal, coeff);
}

int LeastSquares::compute_least_squares_solution(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                 const std::vector<double>& basisVal, double* coeff)
{
  int dim1 = m_numBasis;
  int dim2 = m_numData;

  for(unsigned irhs = 0; irhs < nrhs; ++irhs) {
    // Compute the beta vector, beta = transpose(A) * f
    // Store beta in the coeff vector.
    unsigned fieldIndex = col_major_index(0, irhs, m_numData, m_numComponents);
    unsigned coeffIndex = col_major_index(0, irhs, m_numBasis, nrhs);

    SIERRA_FORTRAN(dgemv)(&TRANSPOSE, &dim2, &dim1, &REAL_ONE, basisVal.data(), &dim2,
                          &fieldVal[fieldIndex], &INT_ONE, &ZERO, &coeff[coeffIndex], &INT_ONE);
  }
  // Solve system of equations using the LU factorization [LAPACK].
  // Replace beta vector with solution(s).
  int dim3 = nrhs;
  int info = 0;
  SIERRA_FORTRAN(dgetrs)(&NOTRANSPOSE, &dim1, &dim3, m_doubleScratchSpace.data(), &dim1,
                         m_integerScratchSpace.data(), coeff, &dim1, &info);
  return info;
}

void LeastSquares::resize_data(const unsigned ncomp, const unsigned nsamp, const unsigned nbasis)
{
  m_numComponents = ncomp;
  m_numSamples = nsamp;
  m_numBasis = nbasis;
  m_numData = ncomp * nsamp;

  m_integerScratchSpace.resize(nbasis);
  m_doubleScratchSpace.resize(nbasis * nbasis);
  m_integerCondScratchSpace.resize(nbasis);
  m_doubleCondScratchSpace.resize(4 * nbasis);
  m_tempVector.resize(ncomp * nsamp);
}

GeometricMovingLeastSquares::GeometricMovingLeastSquares(const unsigned ncomp, const unsigned nsamp,
                                                         const unsigned nbasis, const stk::mesh::BulkData& bulk,
                                                         const std::vector<stk::mesh::Entity>& patch,
                                                         const stk::mesh::FieldBase* coordField,
                                                         const std::vector<double>& evalPoint)
  : LeastSquares(ncomp, nsamp, nbasis)
  , m_patch(patch)
  , m_coordField(coordField)
  , m_evalPoint(evalPoint)
{
  m_weights = compute_weights();
}

int GeometricMovingLeastSquares::least_squares(const unsigned nrhs, const std::vector<double>& fieldVal,
                                               const std::vector<double>& basisVal, double* coeff)
{
  return moving_least_squares(nrhs, fieldVal, basisVal, m_weights, coeff);
}

int GeometricMovingLeastSquares::least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                    const std::vector<double>& basisVal, const double rcond,
                                                    double* coeff)
{
  return moving_least_squares_cond(nrhs, fieldVal, basisVal, m_weights, rcond, coeff);
}

int GeometricMovingLeastSquares::moving_least_squares(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                      const std::vector<double>& basisVal,
                                                      const std::vector<double>& weight, double* coeff)
{
  STK_ThrowRequireMsg(fieldVal.size() >= m_numData * nrhs, "Insufficient size for fieldVal");
  STK_ThrowRequireMsg(basisVal.size() >= m_numData * m_numBasis, "Insufficient size for basisVal");
  STK_ThrowRequireMsg(weight.size() >= m_numData, "Insufficient size for weight");

  fill_covariance_matrix(basisVal, weight);
  return compute_moving_least_squares_solution_with_factorization(nrhs, fieldVal, basisVal, weight, coeff);
}

int GeometricMovingLeastSquares::moving_least_squares_cond(const unsigned nrhs, const std::vector<double>& fieldVal,
                                                           const std::vector<double>& basisVal,
                                                           const std::vector<double>& weight, const double rcond,
                                                           double* coeff)
{
  m_numData = m_numComponents * m_numSamples;

  STK_ThrowRequireMsg(fieldVal.size() >= m_numData * nrhs, "Insufficient size for fieldVal");
  STK_ThrowRequireMsg(basisVal.size() >= m_numData * m_numBasis, "Insufficient size for basisVal");
  STK_ThrowRequireMsg(weight.size() >= m_numData, "Insufficient size for weight");

  fill_covariance_matrix(basisVal, weight);

  double cond;
  int info = compute_factorization_and_condition_number(cond);
  if(info != 0) return info;

  if(cond > rcond) {
    // Condition number is okay.
    // Can proceed with LU
    info = compute_moving_least_squares_solution(nrhs, fieldVal, basisVal, weight, coeff);
  }
  else {
    info = -1;
  }

  return info;
}

void GeometricMovingLeastSquares::fill_covariance_matrix(const std::vector<double>& A, const std::vector<double>& W)
{
  // Get the covariance matrix = transpose(A) * W * A with dimension(m_numBasis,m_numBasis)
  // where matrix "A" has dimension(m_numData,m_numBasis)
  // and diagonal matrix W stored as a vector has dimension(m_numData)
  // There is no real optimized BLAS routine to do this so we make use of the
  // symmetry and diagonality of available matrices
  // Result is stored in m_doubleScratchSpace
  std::fill(m_doubleScratchSpace.begin(), m_doubleScratchSpace.end(), 0.0);

  for(unsigned row = 0; row < m_numData; row++) {
    // Compute M(row,:)'*M(row,:)
    for(unsigned i = 0; i < m_numBasis; i++) {
      for(unsigned j = 0; j < m_numBasis; j++) {
        unsigned resultIndex = col_major_index(i, j, m_numBasis, m_numBasis);
        unsigned sourceRowIndex = col_major_index(row, i, m_numData, m_numBasis);
        unsigned sourceColIndex = col_major_index(row, j, m_numData, m_numBasis);

        m_doubleScratchSpace[resultIndex] += A[sourceRowIndex] * W[row] * A[sourceColIndex];
      }
    }
  }
}

void GeometricMovingLeastSquares::fill_diagonal_scaled_vector(const std::vector<double>& D, const double* f,
                                                              std::vector<double>& x)
{
  // Get the scaled vector x = D * f.
  // where the diagonal matrix D is stored as a vector with dimension(m_numData)
  // and f is a vector with dimension(m_numData)

  for(unsigned i = 0; i < m_numData; ++i) {
    x[i] = D[i] * f[i];
  }
}

int GeometricMovingLeastSquares::compute_moving_least_squares_solution_with_factorization(
    const unsigned nrhs, const std::vector<double>& fieldVal, const std::vector<double>& basisVal,
    const std::vector<double>& weight, double* coeff)
{
  // Replace covariance matrix with its LU factorization
  // using partial pivoting with row interchanges [LAPACK].
  int dim1 = m_numBasis;
  int info = 0;
  SIERRA_FORTRAN(dgetrf)(&dim1, &dim1, m_doubleScratchSpace.data(), &dim1, m_integerScratchSpace.data(), &info);

  if(info != 0) return info;

  return compute_moving_least_squares_solution(nrhs, fieldVal, basisVal, weight, coeff);
}

int GeometricMovingLeastSquares::compute_moving_least_squares_solution(const unsigned nrhs,
                                                                       const std::vector<double>& fieldVal,
                                                                       const std::vector<double>& basisVal,
                                                                       const std::vector<double>& weight, double* coeff)
{
  int dim1 = m_numBasis;
  int dim2 = m_numData;

  for(unsigned irhs = 0; irhs < nrhs; ++irhs) {
    //  Compute the beta vector, beta = transpose(A) * W * f
    //  Store beta in the coeff vector.
    unsigned fieldIndex = col_major_index(0, irhs, m_numData, m_numComponents);
    unsigned coeffIndex = col_major_index(0, irhs, m_numBasis, nrhs);
    fill_diagonal_scaled_vector(weight, &fieldVal[fieldIndex], m_tempVector);

    SIERRA_FORTRAN(dgemv)
    (&TRANSPOSE, &dim2, &dim1, &REAL_ONE, basisVal.data(), &dim2, m_tempVector.data(), &INT_ONE, &ZERO,
        coeff + coeffIndex, &INT_ONE);
  }
  // Solve system of equations using the LU factorization [LAPACK].
  // Replace beta vector with solution(s).
  int dim3 = nrhs;
  int info = 0;
  SIERRA_FORTRAN(dgetrs)(&NOTRANSPOSE, &dim1, &dim3, m_doubleScratchSpace.data(), &dim1,
                         m_integerScratchSpace.data(), coeff, &dim1, &info);
  return info;
}

std::vector<double> GeometricMovingLeastSquares::compute_weights(double EPS)
{
  std::vector<double> weights;
  std::vector<std::vector<double> > samplePoints;
  std::vector<double> centroid;
  for(stk::mesh::Entity entity : m_patch) {
    centroid.clear();
    stk::search::compute_entity_centroid(entity, *m_coordField, centroid);

    STK_ThrowRequireMsg(m_evalPoint.size() == centroid.size(), "Dimension of evalPoint does not match sample points");
    samplePoints.push_back(centroid);
  }

  double maxDistance = 0.0;
  std::vector<double> sampleDistances;
  for(const std::vector<double>& samplePoint : samplePoints) {
    double sampleDistance = stk::search::distance(m_evalPoint.size(), m_evalPoint.data(), samplePoint.data());
    sampleDistances.push_back(sampleDistance);
    maxDistance = std::max(maxDistance, sampleDistance);
  }

  if(maxDistance == 0.0) {
    maxDistance = EPS;
  }

  double radius = 1.05 * maxDistance;
  double maxLocalWeight = 1.0 - EPS;
  for(double sampleDistance : sampleDistances) {
    double localWeight = std::pow(1.0 - sampleDistance / radius, 2) * (1.0 + 2.0 * sampleDistance / radius);
    double clippedWeight = std::min(localWeight, maxLocalWeight);
    double sampleWeight = clippedWeight / (1.0 - clippedWeight);

    for(unsigned i = 0; i < m_numComponents; ++i) {
      weights.push_back(sampleWeight);
    }
  }

  return weights;
}

} // namespace transfer
} // namespace stk
