// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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

  protected:
    
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
