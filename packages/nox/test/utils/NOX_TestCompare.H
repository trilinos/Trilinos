// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_TESTCOMPARE_H
#define NOX_TESTCOMPARE_H

#include "NOX_Common.H"
#include "NOX_Utils.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_MultiVector.H"

namespace NOX {

  //! Utility class to compare values for testing purposes
  class TestCompare {

  public:

    //! Type of comparison
    enum CompareType {

      //! Absolute comparison
      Absolute = 0,

      //! Relative comparison
      Relative = 1
    };

  public:

    //! Constructor
    /*!
     * @param outputStream [in] Stream to print test details
     * @param utilities [in] Utils object to determine verbosity level
     */
    TestCompare(std::ostream& outputStream, const NOX::Utils& utilities);

    //! Destructor
    virtual ~TestCompare();

    //! Compare two values
    /*!
     * @param \c value [in] Computed value
     * @param \c value_expected [in] Value to be compared to
     * @param \c tolerence [in] Tolerance by which values are allowed to differ
     * @param \c name [in] Name for this test when printing details
     * @param \c compareType [in] Type of comparison (Absolute or Relative)
     * @return Returns 0 for is values are within tolerance, 1 otherwise.
     */
    virtual int testValue(double value, double value_expected,
                          double tolerance, const std::string& name,
                          NOX::TestCompare::CompareType compareType);

    //! Compare two vectors
    /*!
     * @param \c vec [in] Computed vector
     * @param \c vec_expected [in] Vector to be compared to
     * @param \c rtol [in] Relative tolerance \f$\epsilon_r\f$
     * @param \c atol [in] Absolute tolerance \f$\epsilon_a\f$
     * @param \c name [in] Name for this test when printing details
     * @return Returns 0 if
     * \f[
     *     \frac{|vec_i - vec_expected_i|}{\epsilon_a + \epsilon_r|vec_expected_i|} < 1
     * \f]
     * for each $i$, and 1 otherwise.
     */
    virtual int testVector(const NOX::Abstract::Vector& vec,
                           const NOX::Abstract::Vector& vec_expected,
                           double rtol, double atol, const std::string& name);

    //! Compare two matrices
    /*!
     * @param \c mat [in] Computed matrix
     * @param \c mat_expected [in] Matrix to be compared to
     * @param \c rtol [in] Relative tolerance \f$\epsilon_r\f$
     * @param \c atol [in] Absolute tolerance \f$\epsilon_a\f$
     * @param \c name [in] Name for this test when printing details
     * @return Returns 0 if
     * \f[
     *     \frac{|mat_{ij} - mat_expected_{ij}|}{\epsilon_a + \epsilon_r|mat_expected_{ij}|} < 1
     * \f]
     * for each $i$, $j$, and 1 otherwise.
     */
    virtual int testMatrix(const NOX::Abstract::MultiVector::DenseMatrix& mat,
                           const NOX::Abstract::MultiVector::DenseMatrix& mat_expected,
                           double rtol, double atol, const std::string& name);

    //! Compare two multi-vectors
    /*!
     * @param \c mvec [in] Computed multi-vector
     * @param \c mvec_expected [in] Multi-vector to be compared to
     * @param \c rtol [in] Relative tolerance \f$\epsilon_r\f$
     * @param \c atol [in] Absolute tolerance \f$\epsilon_a\f$
     * @param \c name [in] Name for this test when printing details
     * @return Returns 0 if
     * \f[
     *     \frac{|mvec_{ij} - mvec_expected_{ij}|}{\epsilon_a + \epsilon_r|mvec_expected_{ij}|} < 1
     * \f]
     * for each $i$, $j$, and 1 otherwise.
     */
    virtual int testMultiVector(const NOX::Abstract::MultiVector& mvec,
                                const NOX::Abstract::MultiVector& mvec_expected,
                                double rtol, double atol, const std::string& name);

    /*!
     * Computes the vector norm
     * \f[
     *     \max_i \frac{|vec_i - vec_expected_i|}{\epsilon_a + \epsilon_r|vec_expected_i|}
     * \f]
     */
    virtual double computeVectorNorm(const NOX::Abstract::Vector& vec,
                                     const NOX::Abstract::Vector& vec_expected,
                                     double rtol, double atol);

  protected:

    //! Output stream to print any messages
    std::ostream& os;

    //! Utils object determing what to print out
    const NOX::Utils utils;

  };

}

#endif
