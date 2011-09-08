// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef THYRA_LINEAR_OP_TESTER_DECL_HPP
#define THYRA_LINEAR_OP_TESTER_DECL_HPP


#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_MultiVectorRandomizerBase.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_FancyOStream.hpp"


namespace Thyra {


/** \brief Testing class for <tt>LinearOpBase</tt>.
 *
 * This testing class performs many different tests just given a
 * <tt>LinearOpBase</tt> object using the function <tt>check()</tt>.
 *
 * This testing class also can check if two linear operators are the same
 * using random vectors by using the function <tt>compare()</tt>.
 *
 * ToDo: Finish documentation!
 *
 * The default compiler-generated copy constructor and assignment operators
 * are allowed since they have the correct semantics which are to simply copy
 * control parameters.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class LinearOpTester {
public:

  /** \brief Local typedef for promoted scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Default constructor which sets default parameter values.
   *
   * See the implementation of this function for the defaults that get set.
   */
  LinearOpTester();
  
  /** \brief Set if to check for linear properties <tt>alpha*op*(x + y) ==
   * op(alpha*x) + op(alpha*y)</tt>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_linear_properties );

  /** \brief Set the tolerance above which a relative error will generate a
   * warning message for the check of the linear properties.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, linear_properties_warning_tol );

  /** \brief Set the tolerance above which a relative error will generate a
   * error message and result in test failure for the check of the linear
   * properties.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, linear_properties_error_tol );

  /** \brief Set if to check for adjoint property <tt>x'*(op*y) ==
   * y'*(op'*x)</tt> if adjoint is supported.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_adjoint  );

  /** \brief Set the tolerance above which a relative error will generate a
   * warning message for the check of the adjoint.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, adjoint_warning_tol );

  /** \brief Set the tolerance above which a relative error will generate a
   * error message and result in test failure for the check of the adjoint.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, adjoint_error_tol );

  /** \brief Set if to check for symmetry property <tt>x'*(op*y) ==
   * y'*(op*x)</tt> for symmetric operators.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, check_for_symmetry  );

  /** \brief Set the tolerance above which a relative error will generate a
   * warning message for the check of symmetry.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, symmetry_warning_tol );

  /** \brief Set the tolerance above which a relative error will generate a
   * error message and result in test failure for the check of symmetry.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( ScalarMag, symmetry_error_tol );

  /** \brief Set the number random vectors that is generated during each test.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_random_vectors );

  /** \brief Set if all tests are shown or just summaries.
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, show_all_tests );

  /** \brief Set if all of the vectors are dumped or not (only relevant if
   * <tt>show_all_tests()==true</tt>).
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, dump_all );

  /** \brief Set the number of right-hand-sides in the multivectors
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( int, num_rhs );

  /** \brief Enable or disable all tests.
   *
   * Postconditions:<ul>
   * <li><tt>this->check_linear_properties()==enable_all_tests</tt>
   * <li><tt>this->check_adjoint()==enable_all_tests</tt>
   * <li><tt>this->check_for_symmetry()==enable_all_tests</tt>
   * </ul>
   */
  void enable_all_tests( const bool enable_all_tests );

  /** \brief Set all the warning tolerances to the same value.
   *
   * Postconditions:<ul>
   * <li><tt>this->linear_properties_warning_tol()==warning_tol</tt>
   * <li><tt>this->adjoint_warning_tol()==warning_tol</tt>
   * <li><tt>this->symmetry_warning_tol()==warning_tol</tt>
   * </ul>
   */
  void set_all_warning_tol( const ScalarMag warning_tol );

  /** \brief Set all the error tolerances to the same value.
   *
   * Postconditions:<ul>
   * <li><tt>this->linear_properties_error_tol()==error_tol</tt>
   * <li><tt>this->adjoint_error_tol()==error_tol</tt>
   * <li><tt>this->symmetry_error_tol()==error_tol</tt>
   * </ul>
   */
  void set_all_error_tol( const ScalarMag error_tol );

  /** \brief Check a linear operator.
   *
   * \param op [in] The linear operator to check.
   *
   * \param out [in/out] If <tt>out!=NULL</tt> then trace output about the
   * tests performed will be sent to <tt>*out</tt>.
   *
   * \param rangeRandomizer [in] Randomizer strategy object for creating
   * random vectors in the range of the operator <tt>op</tt>.  If
   * <tt>NULL</tt> then <tt>UniveralMultiVectorRandomizer</tt> is used intead.
   *
   * \param domainRandomizer [in] Randomizer strategy object for creating
   * random vectors in the domain of the operator <tt>op</tt>.  If
   * <tt>NULL</tt> then <tt>UniveralMultiVectorRandomizer</tt> is used intead.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>rangeRandomizer!=NULL</tt>] <tt>rangeRandomizer->isCompatible(*op.range())==true</tt>
   * <li>[<tt>domainRandomizer!=NULL</tt>] <tt>domainRandomizer->isCompatible(*op.domain())==true</tt>
   * </ul>
   *
   * This function performs a number of tests on <tt>op</tt>:<ul>
   *
   * <li>Checks that the domain and range spaces are valid.
   *
   * <li>Creates temporary vectors using the domain and range spaces.
   *
   * <li>If <tt>this->check_linear_properties()==true</tt> then checks that
   *     \f$\alpha A ( u  + v ) = \alpha A u + \alpha A v\f$ to a
   *     relative tolerance defined by <tt>error_tol()</tt>.  Note, if the client
   *     wants to check the linear properties of the adjoint then the client
   *     should first create an implicit adjoint using <tt>adjoint()</tt>
   *     (or transpose using <tt>transpose()</tt>) which wraps the operation in
   *     a <tt>DefaultScaledAdjointLinearOp</tt>.  Using this method a client can check
   *     all the various values of the enum <tt>EOpTransp</tt>.
   *
   * <li>If <tt>this->check_adjoint()==true</tt> then, checks that the non-transposed
   *     operator and the adjoint operator agree.
   *     The operator and adjoint operator must obey the defined
   *     scalar product.  Specifically, for any two vectors \f$w\f$ (in the
   *     domain space \f$\mathcal{D}\f$) and \f$u\f$ (in the range space
   *     \f$\mathcal{R}\f$), the adjoint operation must obey:
   *     \f[<u,A v>_{\mathcal{R}} = <A^H u, v>_{\mathcal{D}}\f]
   *     where \f$<.,.>_{\mathcal{R}}\f$ is the scalar product defined by
   *     <tt>op.range()->scalarProd()</tt> and \f$<.,.>_{\mathcal{D}}\f$
   *     is the scalar product defined by
   *     <tt>op.domain()->scalarProd()</tt>.
   *
   * <li>If <tt>this->check_for_symmetry()==true</tt> the the operator will
   *     be checked to see if it is symmetric.  Specifically, for any two random
   *     vectors \f$w\f$ and \f$u\f$ in the operator's space \f$\mathcal{S}\f$,
   *     the following property is checked:
   *     \f[<u,A v>_{\mathcal{S}} = <A u, v>_{\mathcal{S}}\f]
   *     where \f$<.,.>_{\mathcal{S}}\f$ is the scalar product defined by
   *     <tt>op.domain()->scalarProd()</tt> and <tt>op.domain()->scalarProd()</tt>.
   *
   * </ul>
   *
   * All relative errors that exceed <tt>xxx_warning_tol()</tt> but do not
   * exceed <tt>xxx_error_tol</tt> will result in special warning messages
   * printed to <tt>*out</tt> (if <tt>out!=NULL</tt>).
   *
   * @return The function returns <tt>true</tt> if all of the tests
   * where within the <tt>xxx_error_tol()</tt> and returns <tt>false</tt>
   * if not.
   *
   * The best way to see what this testing function is doing is to run
   * the test with <tt>out!=NULL</tt> and to look at the
   * implementation by clicking on the following link to the source code:
   */
  bool check(
    const LinearOpBase<Scalar> &op,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &rangeRandomizer,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &domainRandomizer,
    const Ptr<Teuchos::FancyOStream> &out
    ) const;

  /** \brief Calls <tt>this->check(op,null,null,out,leadingIndent,indentSpacer)</tt> */
  bool check(
    const LinearOpBase<Scalar> &op,
    const Ptr<Teuchos::FancyOStream> &out
    ) const;

  /** \brief Check if two linear operators are the same or not.
   *
   * \param op1 [in] The first linear operator
   *
   * \param op2 [in] The second linear operator
   *
   * \param domainRandomizer [in] Randomizer strategy object for creating
   * random vectors in the domain of the operator <tt>op</tt>.  If
   * <tt>NULL</tt> then <tt>UniveralMultiVectorRandomizer</tt> is used intead.
   *
   * \param out [in/out] If <tt>out!=NULL</tt> then trace output about the
   * tests performed will be sent to <tt>*out</tt>.
   *
   * \param leadingIndent [in] All output to <tt>*out</tt> will insert this
   * spacer before each new line is printed.  Default value <tt>""</tt>.
   *
   * \param indentSpacer [in] All output to <tt>*out</tt> that is further
   * indented will use this indentation.  Default value <tt>" "</tt>.
   *
   * This function checks if <tt>op1</tt> and <tt>op2</tt> are the same by
   * checking that the range and domain spaces are compatible and then
   * checking that <tt>sum(op1*v) == sum(op2*v)</tt> for various random
   * vectors.  The allowed warning and error tolerances are taken from
   * <tt>linear_properties_warning_tol()</tt> and
   * <tt>linear_properties_error_tol()</tt>.
   *
   * All relative errors that exceed <tt>xxx_warning_tol()</tt> but do not
   * exceed <tt>xxx_error_tol</tt> will result in special warning messages
   * printed to <tt>*out</tt> (if <tt>out!=NULL</tt>).
   *
   * @return The function returns <tt>true</tt> if all of the tests
   * where within the <tt>xxx_error_tol()</tt> and returns <tt>false</tt>
   * if not.
   *
   * The best way to see what this testing function is doing is to run
   * the test with <tt>out!=NULL</tt> and to look at the
   * implementation by clicking on the following link to the source code:
   */
  bool compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    const Ptr<MultiVectorRandomizerBase<Scalar> > &domainRandomizer,
    const Ptr<Teuchos::FancyOStream> &out_arg
    ) const;
 
  /** \brief Calls
   * <tt>this->compare(op1,op2,NULL,out,leadingIndent,indentSpacer)</tt>.
   */
  bool compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    const Ptr<Teuchos::FancyOStream> &out_arg
    ) const;

  /** \brief Deprecated. */
  //@{

  /** \brief Deprecated. */
  bool check(
    const LinearOpBase<Scalar> &op,
    MultiVectorRandomizerBase<Scalar> *rangeRandomizer,
    MultiVectorRandomizerBase<Scalar> *domainRandomizer,
    Teuchos::FancyOStream *out
    ) const
    {
      using Teuchos::ptr;
      return check(op, ptr(rangeRandomizer), ptr(domainRandomizer), ptr(out));
    }

  /** \brief Deprecated. */
  bool check(
    const LinearOpBase<Scalar> &op,
    Teuchos::FancyOStream *out
    ) const
    {
      return check(op, Teuchos::ptr(out));
    }

  /** \brief Deprecated. */
  bool compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    MultiVectorRandomizerBase<Scalar> *domainRandomizer,
    Teuchos::FancyOStream *out_arg
    ) const
    {
      using Teuchos::ptr;
      return compare(op1, op2, ptr(domainRandomizer), ptr(out_arg));
    }
 
  /** \brief Deprecated. */
  bool compare(
    const LinearOpBase<Scalar> &op1,
    const LinearOpBase<Scalar> &op2,
    Teuchos::FancyOStream *out_arg
    ) const
    {
      return compare(op1, op2, Teuchos::ptr(out_arg));
    }

  //@}

private:

  void setDefaultTols();

}; // class LinearOpTester


} // namespace Thyra


#endif // THYRA_LINEAR_OP_TESTER_DECL_HPP
