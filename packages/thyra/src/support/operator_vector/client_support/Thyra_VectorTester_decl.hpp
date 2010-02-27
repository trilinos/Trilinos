// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_VECTOR_TESTER_DECL_HPP
#define THYRA_VECTOR_TESTER_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_MultiVectorTester.hpp"

namespace Thyra {

/** \brief Unit testing class for a <tt>VectorBase</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class VectorTester {
public:

  /** \brief Local typedef for scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Expose a non-<tt>const</tt> reference to the <tt>MultiVectorTester</tt>
   * object used to test the <tt>MultiVectorBase</tt> interface supported by
   * <tt>VectorBase</tt>.
   *
   * Clients can use this interface to override options directly set on
   * <tt>*this</tt>.
   */
  MultiVectorTester<Scalar>& multiVectorTester();

  /** \brief Expose a <tt>const</tt> reference to the <tt>MultiVectorTester</tt>
   * object used to test the <tt>MultiVectorBase</tt> interface supported by
   * <tt>VectorBase</tt>.
   *
   * Clients can use this interface to query (but not change) options.
   */
  const MultiVectorTester<Scalar>& multiVectorTester() const;

  /** \brief Set the tolerance above which a relative error will generate a
   * warning message.
   *
   * Also calls <tt>this->multiVectorTester().set_all_warning_tol(warning_tol).
   */
  void warning_tol( const ScalarMag &warning_tol );

  /** \brief Return the warning tolerance for <tt>*this</tt>. */
  ScalarMag warning_tol() const;

  /** \brief Set the error above which a relative error will generate a
   * an message and cause the test to fail.
   *
   * Also calls <tt>this->multiVectorTester().set_all_error_tol(error_tol).
   */
  void error_tol( const ScalarMag &error_tol );

  /** \brief Return the error tolerance for <tt>*this</tt>. */
  ScalarMag error_tol() const;

  /** \brief Set the number random vectors that is generated during each test.
   *
   * Also calls <tt>this->multiVectorTester().num_random_vectors(num_random_vectors).
   */
  void num_random_vectors( const int num_random_vectors );

  /** \brief Return the number of random vectors used for <tt>*this</tt> objects
   * tests.
   */
  int num_random_vectors() const;

  /** \brief Set whether all of the tests will be printed independent if they
   * pass or fail.
   *
   * Also calls <tt>this->multiVectorTester().show_all_tests(show_all_tests).
   */
  void show_all_tests( const bool show_all_tests );

  /** \brief Return the number of random vectors used for <tt>*this</tt> objects
   * tests.
   */
  bool show_all_tests() const;

  /** \brief Set whether all of the vectors and multi-vectors will be dumped
   * or not.
   *
   * Also calls <tt>this->multiVectorTester().dump_all(dump_all).
   */
  void dump_all( const bool dump_all );

  /** \brief Return the number of random vectors used for <tt>*this</tt> objects
   * tests.
   */
  bool dump_all() const;

  /** \brief Default constructor which sets default parameter values.
   *
   * Note: It is not recommended that the client pass in values in this
   * constructor since the argument list may change in the near future but
   * instead use the above set functions to change an option after
   * construction.
   *
   * Postconditions:<ul>
   * <li><tt>this->warning_tol() == warning_tol</tt>
   * <li><tt>this->error_tol() == error_tol</tt>
   * <li><tt>this->num_random_vectors() == num_random_vectors</tt>
   * <li><tt>this->show_all_tests() == show_all_tests</tt>
   * <li><tt>this->dump_all() == dump_all</tt>
   * <li>Calls <tt>this->multiVectorTester().set_all_warning_tol(warning_tol)</tt>
   * <li>Calls <tt>this->multiVectorTester().set_all_error_tol(error_tol)</tt>
   * <li>Calls <tt>this->multiVectorTester().num_random_vectors(num_random_vectors)</tt>
   * <li>Calls <tt>this->multiVectorTester().show_all_tests(show_all_tests)</tt>
   * <li>Calls <tt>this->multiVectorTester().dump_all(dump_all)</tt>
   * </ul>
   */
  VectorTester(
    const ScalarMag     warning_tol            = 1e-13
    ,const ScalarMag    error_tol              = 1e-10
    ,const int          num_random_vectors     = 1
    ,const bool         show_all_tests         = false
    ,const bool         dump_all               = false
    );

  /** \brief Check a vector object in a set of comprehensive tests.
   *
   * @param  v      [in] The vector object to test.
   * @param  out    [in/out] If <tt>out != NULL</tt> then output will be sent to <tt>*out</tt>.
   *
   * The behavior of this function greatly depends on a number of options (see
   * <tt>VectorTester()</tt> for the default values for these options):
   *
   * <ul>
   * <li> <b><tt>print_all_tests(bool)</tt></b>:  If <tt>print_all_tests() == true</tt>, then some output will be sent to
   *      <tt>*out</tt> for every test performed.  This is useful to see all of tests that are performed and
   *      in debugging.
   * <li> <b><tt>dump_all(bool)</tt></b>:  If <tt>dump_all() == true</tt>, then all of the vectors will be printed
   *      that are created during the tests.  This option is really only needed during initial debugging
   *      and should only be used with small vector spaces since it will produce a lot of <tt>O(space.dim())</tt>
   *      output.
   * <li> <b><tt>num_random_tests(int)</tt></b>:  This is the number of random tests to perform per category of test.
   *      A higher number will result is better validation but will consume more CPU time.
   * <li> <b><tt>warning_tol(ScalarMag)</tt></b>:  Any test with a relative error greater than <tt>warning_tol()</tt> will
   *      result in a warning message printed to <tt>*out</tt> but will not result in a filed test.
   * <li> <b><tt>error_tol(Scalar)</tt></b>:  Any test with a relative error greater than <tt>error_tol()</tt> will
   *      result in an error message printed to <tt>*out</tt> and will result in a failed test.
   * </ul>
   *
   * @return The function returns <tt>true</tt> if all of the tests where
   * within the <tt>error_tol()</tt> and returns <tt>false</tt> if not.
   *
   * The best way to see what this testing function is doing is to run the
   * test with <tt>out!=NULL</tt> and to look at the implementation by
   * clicking on the following link to the source code:
   */
  bool check(
    const VectorBase<Scalar>       &v
    ,Teuchos::FancyOStream         *out
    ) const;

private:

  MultiVectorTester<Scalar> multiVectorTester_;

  ScalarMag    warning_tol_;
  ScalarMag    error_tol_;
  int          num_random_vectors_;
  bool         show_all_tests_;
  bool         dump_all_;

}; // class VectorTester

// //////////////////////////////////
// Inline members

template<class Scalar>
inline
MultiVectorTester<Scalar>& VectorTester<Scalar>::multiVectorTester()
{
  return multiVectorTester_;
}

template<class Scalar>
inline
const MultiVectorTester<Scalar>& VectorTester<Scalar>::multiVectorTester() const
{
  return multiVectorTester_;
}

template<class Scalar>
inline
void VectorTester<Scalar>::warning_tol( const ScalarMag &warning_tol_in )
{
  warning_tol_ = warning_tol_in;
  multiVectorTester_.warning_tol(warning_tol_in);
}

template<class Scalar>
inline
typename VectorTester<Scalar>::ScalarMag
VectorTester<Scalar>::warning_tol() const
{
  return warning_tol_;
}

template<class Scalar>
inline
void VectorTester<Scalar>::error_tol( const ScalarMag &error_tol_in )
{
  error_tol_ = error_tol_in;
  multiVectorTester_.error_tol(error_tol_in);
}

template<class Scalar>
inline
typename VectorTester<Scalar>::ScalarMag
VectorTester<Scalar>::error_tol() const
{
  return error_tol_;
}

template<class Scalar>
inline
void VectorTester<Scalar>::num_random_vectors( const int num_random_vectors_in )
{
  num_random_vectors_ = num_random_vectors_in;
  multiVectorTester_.num_random_vectors(num_random_vectors_in);
}

template<class Scalar>
inline
int VectorTester<Scalar>::num_random_vectors() const
{
  return num_random_vectors_;
}

template<class Scalar>
inline
void VectorTester<Scalar>::show_all_tests( const bool show_all_tests_in )
{
  show_all_tests_ = show_all_tests_in;
  multiVectorTester_.show_all_tests(show_all_tests_in);
}

template<class Scalar>
inline
bool VectorTester<Scalar>::show_all_tests() const
{
  return show_all_tests_;
}

template<class Scalar>
inline
void VectorTester<Scalar>::dump_all( const bool dump_all_in )
{
  dump_all_ = dump_all_in;
  multiVectorTester_.dump_all(dump_all_in);
}

template<class Scalar>
inline
bool VectorTester<Scalar>::dump_all() const
{
  return dump_all_;
}

} // namespace Thyra

#endif // THYRA_VECTOR_TESTER_DECL_HPP
