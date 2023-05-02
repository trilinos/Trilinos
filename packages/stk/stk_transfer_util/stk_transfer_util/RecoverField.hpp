/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_TRANSFER_RECOVER_FIELD_HPP
#define STK_TRANSFER_RECOVER_FIELD_HPP

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_transfer_util/LeastSquares.hpp>
#include <stk_transfer_util/Patch.hpp>              // for entity_patch_if
#include <cmath>                                     // for fabs
#include <limits>                                    // for numeric_limits
#include <ostream>                                   // for endl, etc
#include <stk_util/util/Fortran.hpp>                 // for SIERRA_FORTRAN
#include <vector>                                    // for vector
#include "stk_mesh/base/BulkData.hpp"                // for BulkData
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_util/environment/Env.hpp"              // for parallel_rank
#include "stk_util/util/ReportHandler.hpp"    // for STK_ThrowRequireMsg, etc
#include "stk_util/environment/RuntimeMessage.hpp"   // for MessageCode
#include "stk_util/environment/RuntimeWarning.hpp"
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################
/**
 * Support for least-squares recovery of a low-order field
 *
 * @par Description
 * This class provides an algorithm for computing a recovered (i.e.,
 * higher-order) field from a lower-order field. The lower-order field
 * is associated with an element, face, or edge. The recovered field
 * is returned in form of coefficients of a given polynomial basis.
 *
 * @par Requirements
 * Before you call the algorithm, you must construct a
 * RecoverField::RecoverInfo object, and you must implement a
 * derived class using RecoverField as the bases class. The
 * derived class must implement the 'sample_patch' virtual method.
 */

namespace stk {
namespace transfer {

inline void evaluate_trilinear_basis(const double x, const double y, const double z, double* basisSample, int stride)
{
  basisSample[0         ] = 1.0;

  basisSample[    stride] = x;
  basisSample[2 * stride] = y;
  basisSample[3 * stride] = z;

  basisSample[4 * stride] = x * y;
  basisSample[5 * stride] = y * z;
  basisSample[6 * stride] = x * z;

  basisSample[7 * stride] = x * y * z;
}

inline void evaluate_triquadratic_basis(const double x, const double y, const double z, double* basisSample, int stride)
{
  basisSample[0         ] = 1.0;

  basisSample[    stride] = x;
  basisSample[2 * stride] = y;
  basisSample[3 * stride] = z;

  basisSample[4 * stride] = x * y;
  basisSample[5 * stride] = y * z;
  basisSample[6 * stride] = x * z;
  basisSample[7 * stride] = x * x;
  basisSample[8 * stride] = y * y;
  basisSample[9 * stride] = z * z;

  basisSample[10 * stride] = x * y * y;
  basisSample[11 * stride] = x * z * z;
  basisSample[12 * stride] = x * x * y;
  basisSample[13 * stride] = y * z * z;
  basisSample[14 * stride] = x * x * z;
  basisSample[15 * stride] = y * y * z;
  basisSample[16 * stride] = x * y * z;

  basisSample[17 * stride] = x * y * z * z;
  basisSample[18 * stride] = x * y * y * z;
  basisSample[19 * stride] = x * x * y * z;
  basisSample[20 * stride] = x * x * y * y;
  basisSample[21 * stride] = y * y * z * z;
  basisSample[22 * stride] = x * x * z * z;

  basisSample[23 * stride] = x * y * y * z * z;
  basisSample[24 * stride] = x * x * y * z * z;
  basisSample[25 * stride] = x * x * y * y * z;

  basisSample[26 * stride] = x * x * y * y * z * z;
}

inline void evaluate_tricubic_basis(const double x, const double y, const double z, double* basisSample, int stride)
{
  basisSample[0         ] = 1.0;

  basisSample[    stride] = x;
  basisSample[2 * stride] = y;
  basisSample[3 * stride] = z;

  basisSample[ 4 * stride] = x * x;
  basisSample[ 5 * stride] = x * y;
  basisSample[ 6 * stride] = x * z;
  basisSample[ 7 * stride] = y * y;
  basisSample[ 8 * stride] = y * z;
  basisSample[ 9 * stride] = z * z;

  basisSample[10 * stride] = x * x * x;
  basisSample[11 * stride] = x * x * y;
  basisSample[12 * stride] = x * x * z;
  basisSample[13 * stride] = x * y * y;
  basisSample[14 * stride] = x * y * z;
  basisSample[15 * stride] = x * z * z;
  basisSample[16 * stride] = y * y * y;
  basisSample[17 * stride] = y * y * z;
  basisSample[18 * stride] = y * z * z;
  basisSample[19 * stride] = z * z * z;

  basisSample[20 * stride] = x * x * x * y;
  basisSample[21 * stride] = x * x * x * z;
  basisSample[22 * stride] = x * x * y * y;
  basisSample[23 * stride] = x * x * y * z;
  basisSample[24 * stride] = x * x * z * z;
  basisSample[25 * stride] = x * y * y * y;
  basisSample[26 * stride] = x * y * y * z;
  basisSample[27 * stride] = x * y * z * z;
  basisSample[28 * stride] = x * z * z * z;
  basisSample[29 * stride] = y * y * y * z;
  basisSample[30 * stride] = y * y * z * z;
  basisSample[31 * stride] = y * z * z * z;

  basisSample[32 * stride] = x * x * x * y * y;
  basisSample[33 * stride] = x * x * x * y * z;
  basisSample[34 * stride] = x * x * x * z * z;
  basisSample[35 * stride] = x * x * y * y * y;
  basisSample[36 * stride] = x * x * y * y * z;
  basisSample[37 * stride] = x * x * y * z * z;
  basisSample[38 * stride] = x * x * z * z * z;
  basisSample[39 * stride] = x * y * y * y * z;
  basisSample[40 * stride] = x * y * y * z * z;
  basisSample[41 * stride] = x * y * z * z * z;
  basisSample[42 * stride] = y * y * y * z * z;
  basisSample[43 * stride] = y * y * z * z * z;

  basisSample[44 * stride] = x * x * x * y * y * y;
  basisSample[45 * stride] = x * x * x * y * y * z;
  basisSample[46 * stride] = x * x * x * y * z * z;
  basisSample[47 * stride] = x * x * x * z * z * z;
  basisSample[48 * stride] = x * x * y * y * y * z;
  basisSample[49 * stride] = x * x * y * y * z * z;
  basisSample[50 * stride] = x * x * y * z * z * z;
  basisSample[51 * stride] = x * y * y * y * z * z;
  basisSample[52 * stride] = x * y * y * z * z * z;
  basisSample[53 * stride] = y * y * y * z * z * z;

  basisSample[54 * stride] = x * x * x * y * y * y * z;
  basisSample[55 * stride] = x * x * x * y * y * z * z;
  basisSample[56 * stride] = x * x * x * y * z * z * z;
  basisSample[57 * stride] = x * x * y * y * y * z * z;
  basisSample[58 * stride] = x * x * y * y * z * z * z;
  basisSample[59 * stride] = x * y * y * y * z * z * z;

  basisSample[60 * stride] = x * x * x * y * y * y * z * z;
  basisSample[61 * stride] = x * x * x * y * y * z * z * z;
  basisSample[62 * stride] = x * x * y * y * y * z * z * z;

  basisSample[63 * stride] = x * x * x * y * y * y * z * z * z;
}

class RecoverField {
 public:
  /**
   * The value of RecoveryType is equal to the number of unknown coefficients in the
   * recovered polynomial. The assumed ordering for these coefficients is as indicated.
   */
  enum RecoveryType {
    LINEAR_2D = 3,     ///< a_0 + a_1*x + a_2*y
    LINEAR_3D = 4,     ///< a_0 + a_1*x + a_2*y + a_3*z
    BILINEAR = 4,      ///< a_0 + a_1*x + a_2*y + a_3*x*y
    TRILINEAR = 8,     ///< see evaluate_trilinear_basis function above
    TRIQUADRATIC = 27, ///< see evaluate_triquadratic_basis function above
    TRICUBIC = 64,     ///< see evaluate_tricubic_basis function above
    MAX_COEFF = 64
  };

  /** Typedef for pairing field object (the variable to be recovered) with stateissue */

  /**
   * Constructor for the recovery of a field, assuming the recovered
   * function will be stored as a (non-registered) list of polynomial basis coefficients
   *
   * @param rec_type the recovery type (of type RecoverField::RecoveryType)
   * @param rec_vars the (lower-order) list of fields (with states) we wish to compute a recovery for
   * @param rec_node_var the field variable (and state) containing geometric node coordinates
   * @param nsamp_elem the number of sampling points per patch object
   */
  RecoverField(RecoverField::RecoveryType recType, int nSampElem)
    : m_nSampleElements(nSampElem)
    , m_totalVarComponents(0)
    , m_recoveryType(recType)
  {

  }

  virtual ~RecoverField() {}

  /**
   * Perform a least-squares recovery of a field onto a function defined by
   * info.recovery_type(). The recovery_type is NOT guaranteed. That is, if the
   * least-squares projection matrix is rank-deficient (as estimated using the
   * value of diagonal_tol), then a projection onto a
   * lower-order polynomial will systematically be attempted. If all efforts fail
   * at a projection, then the function returns 0; otherwise the rank of the
   * matrix used to determine the coefficients is returned.
   * If projection onto a lower-order polynomial was necessary, then those
   * coefficients corresponding to the removed-monomials will be zero. If the
   * projection failed (returning 0), then all coefficients will be NAN.
   *
   * This method calls 'sample_patch', which is an abstract virtual method in
   * this class. It must be implemented in an application-derived class.
   *
   * To use this method, support for dynamic mesh modifications must be enabled
   * in the region, via the call region->support_dynamic_mesh_modifications(true);
   * (this constructs the back-relations, i.e., node-to-element connectivity).
   *
   * If diagonalTol is zero, then the check for a rank-deficient system will
   * be skipped, potentially making the algorithm more efficient (i.e., the
   * condition number will not be computed). However, such a choice may also
   * cause the algorithm to fail if the condition number actually is bad.
   */
  template <typename FILTER>
  int compute_recovered_coeff(const stk::mesh::BulkData& stkMesh, stk::mesh::Entity entity, const FILTER& patchFilter,
                              const stk::mesh::Selector& selector,
                              const int numRecCoeff, double* recoveredCoeff, const double diagonalTol = 0.0) const
  {
    STK_ThrowRequireMsg(diagonalTol >= 0.0, "The value supplied for 'diagonalTol' is negative\n");

    STK_ThrowRequireMsg(diagonalTol <= 1.0,
                    "The value supplied for 'diagonalTol' is greater than one; therefore it is not usable\n"
                        << "Please supply a tolerance greater than or equal to 0.0 and less than or equal to 1.0.\n\n");

    // Construct the object-centered patch
    std::vector<stk::mesh::Entity> patch;
    std::vector<stk::mesh::Entity> nodesInPatch;

    entity_patch_if(stkMesh, entity, patchFilter, selector, patch, nodesInPatch, stkMesh.entity_rank(entity));

    STK_ThrowRequireMsg(!patch.empty(), "Mesh entity " << stkMesh.identifier(entity) << " on processor '"
                                                   << sierra::Env::parallel_rank() << "' has an empty patch\n");

    const int numSampPtsPatch = patch.size() * m_nSampleElements;
    const int basisSize = (int)m_recoveryType;

    STK_ThrowRequireMsg(numRecCoeff == basisSize * m_totalVarComponents,
                    "numRecCoeff is " << numRecCoeff << " which is incompatible with expected value '"
                                      << basisSize * m_totalVarComponents << "\n");

    int ncomp = 1;
    LeastSquares leastSquaresCalculator(ncomp, numSampPtsPatch, basisSize);

    return recover_coeff(stkMesh, entity, patch, leastSquaresCalculator, recoveredCoeff, diagonalTol);
  }

  int recover_coeff(const stk::mesh::BulkData& stkMesh, stk::mesh::Entity entity,
                    const std::vector<stk::mesh::Entity>& patch, LeastSquares& leastSquaresCalculator,
                    double* recoveredCoeff, const double inputDiagonalTol = 0.0) const
  {
    const int numSampPtsPatch = leastSquaresCalculator.get_num_samples();
    const int basisSize = leastSquaresCalculator.get_num_basis();
    int rank = basisSize;

    // Allocate arrays and variables
    int ncomp = leastSquaresCalculator.get_num_components();
    std::vector<double> fieldSample(numSampPtsPatch * m_totalVarComponents);
    std::vector<double> basisSample(numSampPtsPatch * basisSize);

    int ierr = -1;

    // 'sample_patch' implemented by application-derived class
    sample_patch(patch, numSampPtsPatch, fieldSample.data(), basisSample.data());

    if(inputDiagonalTol == 0.0 /* skip check of condition number */) {
      ierr = leastSquaresCalculator.least_squares(m_totalVarComponents, fieldSample, basisSample, recoveredCoeff);

      if(0 != ierr) {
        static stk::MessageCode report_id(1);
        stk::RuntimeWarning(report_id) << "Least-squares projection failed for "
                                       << "Mesh entity '" << stkMesh.identifier(entity) << "' on processor '"
                                       << sierra::Env::parallel_rank() << "'\n" << "Suggest recomputing with nonzero 'inputDiagonalTol'\n";

        double realMax = std::numeric_limits<double>::max();
        for(int kk = 0; kk < basisSize * m_totalVarComponents; ++kk)
          recoveredCoeff[kk] = realMax;

        return 0;
      }
    }

    else // inputDiagonalTol > 0.0; perform condition-number check
    {
      double rcond = 1.0e-10;
      ierr =
          leastSquaresCalculator.least_squares_cond(m_totalVarComponents, fieldSample, basisSample, rcond, recoveredCoeff);

      std::vector<double>& work = leastSquaresCalculator.get_double_scratch_space();

      if(0 != ierr) {
        // One or more columns is linearly dependent

        // Scale diagonalTol by largest diagonal

        double diagonalTol(-1);
        for(int i(0); i < basisSize; ++i) {
          diagonalTol =
              diagonalTol < std::fabs(work[i * (1 + basisSize)]) ? std::fabs(work[i * (1 + basisSize)]) : diagonalTol;
        }

        STK_ThrowAssert(0 < diagonalTol);

        diagonalTol *= inputDiagonalTol;

        std::vector<int> badCols(basisSize);

        int numBad(0);

        // Loop over diagonals in work.
        // Zero entries are bad
        for(int i(0); i < basisSize; ++i) {
          if(diagonalTol > std::fabs(work[i * (1 + basisSize)])) {
            badCols[numBad++] = i;
          }
        }

        if(numBad == basisSize /* bail out here */) {
          static stk::MessageCode report_id(1);
          stk::RuntimeWarning(report_id) << "Least-squares projection failed for "
                                         << "Mesh entity '" << stkMesh.identifier(entity) << " on processor '"
                                         << stkMesh.parallel_rank() << "'" << std::endl
                                         << "The number of bad columns:" << numBad
                                         << " is equal to the basis size:" << basisSize << std::endl
                                         << "The largest scaled diagonal entry:" << diagonalTol
                                         << " ascale value:" << inputDiagonalTol << std::endl;
          double realMax = std::numeric_limits<double>::max();
          for(int kk = 0; kk < basisSize * m_totalVarComponents; ++kk)
            recoveredCoeff[kk] = realMax;

          return 0;
        }
        else if(numBad) {
          // Remove the bad column by shifting the other columns to the left.
          // iCurrent -> the column into which another column will be copied
          // iCheck   -> the column we are checking
          // iBad     -> the entry into badCols, badCols[iBad] = next bad column
          for(int iCurrent(badCols[0]), iCheck(iCurrent + 1), iBad(1); iCheck < basisSize; ++iCheck) {
            if(iBad < numBad && iCheck == badCols[iBad]) {
              // The iCheck column is a bad one.
              // Increment iBad so we check the next entry in badCols next time.
              ++iBad;
            }
            else {
              // The iCheck column is a good one.
              double* good = &basisSample[numSampPtsPatch * iCheck];

              // Get the column to fill.
              double* current = &basisSample[numSampPtsPatch * iCurrent];

              // Increment iCurrent so we fill the next column next time.
              ++iCurrent;

              // Copy the good column to the current column.
              for(int k(0); k < numSampPtsPatch; ++k) {
                current[k] = good[k];
              }
            }
          }
        }

        rank = basisSize - numBad;

        leastSquaresCalculator.resize_data(ncomp, numSampPtsPatch, rank);
        ierr = leastSquaresCalculator.least_squares(m_totalVarComponents, fieldSample, basisSample, recoveredCoeff);

        // If the following if-test ever fails, perhaps we should call
        //  least_squares_cond() above instead of least_squares()
        if(0 != ierr) {
          static stk::MessageCode report_id(1);
          stk::RuntimeWarning x(report_id);
          x << "Least-squares projection failed for "
            << "Mesh entity '" << stkMesh.identifier(entity) << sierra::Env::parallel_rank() << "'" << std::endl
            << "The Error code returned from least_squares:" << ierr << " is not zero." << std::endl
            << "The rank of the offending matrix is:" << rank << std::endl
            << "The number of sampling points in patch:" << numSampPtsPatch << std::endl
            << std::endl << "Basis Matrix:\n";
          for(int j(0); j < numSampPtsPatch; ++j) {
            for(int i(0); i < rank; ++i) {
              x << basisSample[i * numSampPtsPatch + j] << " ";
            }
            x << std::endl;
          }
          x << "Result of (Transpose Basis)*(Basis) Matrix.\n"
            << "If a positive error value was returned from least_squares\n"
            << "then this matrix somehow became singular which is an error.\n";
          for(int i(0); i < rank; ++i) {
            for(int j(0); j < rank; ++j) {
              x << work[i * rank + j] << " ";
            }
            x << std::endl;
          }
          x << std::endl;
          double realMax = std::numeric_limits<double>::max();
          for(int kk = 0; kk < basisSize * m_totalVarComponents; ++kk)
            recoveredCoeff[kk] = realMax;

          return 0;
        }

        // We now have rank coefficients.
        // Put them into the correct positions with zeroes elsewhere.

        // For each component
        for(int i(m_totalVarComponents - 1); i > -1; --i) {
          double* pFinalCoeff = &recoveredCoeff[basisSize * i];
          double* pCurrentCoeff = &recoveredCoeff[rank * i];

          // Copy values.  Use work as scratch space.
          for(int j(0); j < rank; ++j) {
            work[j] = pCurrentCoeff[j];
          }

          // iCheck   -> the column we are checking
          // iBad     -> the entry into badCols, badCols[iBad] = next bad column
          int iGood(0);
          for(int iCheck(0), iBad(0); iCheck < basisSize; ++iCheck) {
            if(iBad < numBad && iCheck == badCols[iBad]) {
              // The iCheck row is a bad one.
              pFinalCoeff[iCheck] = 0.0;
              // Increment iBad so we check the next entry in badCols next time.
              ++iBad;
            }
            else {
              // The iCheck row is a good one.
              pFinalCoeff[iCheck] = work[iGood++];
            }
          }

          STK_ThrowAssert(iGood == rank);
        }

      } // end if linearly dependent columns
    }   // end if diagonal_tol

    return rank;
  }

  /**
   * Given a patch of mesh objects, compute the values of the to-be-recovered field and
   * the basis at each of the sampling points. This method must be implemented in a
   * derived class. The derived class has knowledge of the sampling points and how to
   * compute the field and basis values. The base class, along with 'nsamp_patch' which
   * is passed in, contains the information in its protected data on the dimensions of
   * field_sample and basis_sample arrays.
   *
   * @param patch the list of mesh objects in the patch
   * @param nsamp_patch total number of sampling points in the patch (= patch.size() * m_nSampleElements)
   * @param field_sample array of field values to be computed; dimension of array is (fortran-ordered)
   *        (nsamp_patch,m_totalVarComponents)
   * @param basis_sample array of basis values to be computed; dimension of array is (fortran-ordered)
   *        (nsamp_patch*recoveryType), where 'm_recoveryType' is the number of basis vectors including the
   *        constant
   */
  virtual void sample_patch(const std::vector<stk::mesh::Entity>& patch, int nsampPatch, double* fieldSample,
                            double* basisSample) const = 0;

  /// Query total number of field components to be recovered (e.g., == 3 for a 3D vector, == 6 for two 3D vectors, etc.)
  int total_recov_var_comp() const { return m_totalVarComponents; }

 protected:
  int m_nSampleElements;
  int m_totalVarComponents;

  RecoveryType m_recoveryType;
};


} // namespace transfer
} // namespace stk

#endif // STK_TRANSFER_RECOVER_FIELD_HPP
