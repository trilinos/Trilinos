// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_TESTER_DECL_HPP
#define THYRA_MULTI_VECTOR_TESTER_DECL_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Teuchos_Describable.hpp"

namespace Thyra {

/** \brief Unit testing class for a <tt>MultiVectorBase</tt> object.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template <class Scalar>
class MultiVectorTester : public Teuchos::Describable {
 public:
  /** \brief Local typedef for scalar magnitude */
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  /** \brief Expose a non-<tt>const</tt> reference to the <tt>LinearOpTester</tt>
   * object used to test the <tt>LinearOpBase</tt> interface supported by
   * <tt>MultiVectorBase</tt>.
   *
   * Clients can use this interface to override options directly set on
   * <tt>*this</tt>.
   */
  LinearOpTester<Scalar> &linearOpTester();

  /** \brief Expose a <tt>const</tt> reference to the <tt>LinearOpTester</tt>
   * object used to test the <tt>LinearOpBase</tt> interface supported by
   * <tt>MultiVectorBase</tt>.
   *
   * Clients can use this interface to query (but not change) options.
   */
  const LinearOpTester<Scalar> &linearOpTester() const;

  /** \brief Set the tolerance above which a relative error will generate a
   * warning message.
   *
   * Also calls <tt>this->linearOpTester().set_all_warning_tol(warning_tol).
   */
  void warning_tol(const ScalarMag &warning_tol);

  /** \brief Return the warning tolerance for <tt>*this</tt>. */
  ScalarMag warning_tol() const;

  /** \brief Set the error above which a relative error will generate a
   * an message and cause the test to fail.
   *
   * Also calls <tt>this->linearOpTester().set_all_error_tol(error_tol).
   */
  void error_tol(const ScalarMag &error_tol);

  /** \brief Return the error tolerance for <tt>*this</tt>. */
  ScalarMag error_tol() const;

  /** \brief Set the number random vectors that is generated during each test.
   *
   * Also calls <tt>this->linearOpTester().num_random_vectors(num_random_vectors).
   */
  void num_random_vectors(const int num_random_vectors);

  /** \brief Return the number of random vectors used for <tt>*this</tt> objects
   * tests.
   */
  int num_random_vectors() const;

  /** \brief Set whether all of the tests will be printed independent if they
   * pass or fail.
   *
   * Also calls <tt>this->linearOpTester().show_all_tests(show_all_tests).
   */
  void show_all_tests(const bool show_all_tests);

  /** \brief Return the number of random vectors used for <tt>*this</tt> objects
   * tests.
   */
  bool show_all_tests() const;

  /** \brief Set whether all of the vectors and multi-vectors will be dumped
   * or not.
   *
   * Also calls <tt>this->linearOpTester().dump_all(dump_all).
   */
  void dump_all(const bool dump_all);

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
   * <li>Calls <tt>this->linearOpTester().set_all_warning_tol(warning_tol)</tt>
   * <li>Calls <tt>this->linearOpTester().set_all_error_tol(error_tol)</tt>
   * <li>Calls <tt>this->linearOpTester().num_random_vectors(num_random_vectors)</tt>
   * <li>Calls <tt>this->linearOpTester().show_all_tests(show_all_tests)</tt>
   * <li>Calls <tt>this->linearOpTester().dump_all(dump_all)</tt>
   * </ul>
   */
  MultiVectorTester(
      const ScalarMag warning_tol  = 1e-13,
      const ScalarMag error_tol    = 1e-10,
      const int num_random_vectors = 1,
      const bool show_all_tests    = false,
      const bool dump_all          = false);

  /** \brief Check a multi-vector as created by a VectorSpaceBase object. */
  bool checkMultiVector(const VectorSpaceBase<Scalar> &vs,
                        const Ptr<Teuchos::FancyOStream> &out) const;

  /** \brief Check a multi-vector object in a set of comprehensive teats.
   *
   * @param  mv     [in] The multi-vector object to test.
   * @param  out    [in/out] If <tt>out != NULL</tt> then output will be sent to <tt>*out</tt>.
   *
   * The behavior of this function greatly depends on a number of options (see
   * <tt>MultiVectorTester()</tt> for the default values for these options):
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
      const MultiVectorBase<Scalar> &mv,
      const Ptr<Teuchos::FancyOStream> &out) const;

 private:
  LinearOpTester<Scalar> linearOpTester_;

  ScalarMag warning_tol_;
  ScalarMag error_tol_;
  int num_random_vectors_;
  bool show_all_tests_;
  bool dump_all_;

};  // class MultiVectorTester

// //////////////////////////////////
// Inline members

template <class Scalar>
inline LinearOpTester<Scalar> &MultiVectorTester<Scalar>::linearOpTester() {
  return linearOpTester_;
}

template <class Scalar>
inline const LinearOpTester<Scalar> &MultiVectorTester<Scalar>::linearOpTester() const {
  return linearOpTester_;
}

template <class Scalar>
inline void MultiVectorTester<Scalar>::warning_tol(const ScalarMag &warning_tol_in) {
  warning_tol_ = warning_tol_in;
  linearOpTester_.set_all_warning_tol(warning_tol_in);
}

template <class Scalar>
inline
    typename MultiVectorTester<Scalar>::ScalarMag
    MultiVectorTester<Scalar>::warning_tol() const {
  return warning_tol_;
}

template <class Scalar>
inline void MultiVectorTester<Scalar>::error_tol(const ScalarMag &error_tol_in) {
  error_tol_ = error_tol_in;
  linearOpTester_.set_all_error_tol(error_tol_in);
}

template <class Scalar>
inline
    typename MultiVectorTester<Scalar>::ScalarMag
    MultiVectorTester<Scalar>::error_tol() const {
  return error_tol_;
}

template <class Scalar>
inline void MultiVectorTester<Scalar>::num_random_vectors(const int num_random_vectors_in) {
  num_random_vectors_ = num_random_vectors_in;
  linearOpTester_.num_random_vectors(num_random_vectors_in);
}

template <class Scalar>
inline int MultiVectorTester<Scalar>::num_random_vectors() const {
  return num_random_vectors_;
}

template <class Scalar>
inline void MultiVectorTester<Scalar>::show_all_tests(const bool show_all_tests_in) {
  show_all_tests_ = show_all_tests_in;
  linearOpTester_.show_all_tests(show_all_tests_in);
}

template <class Scalar>
inline bool MultiVectorTester<Scalar>::show_all_tests() const {
  return show_all_tests_;
}

template <class Scalar>
inline void MultiVectorTester<Scalar>::dump_all(const bool dump_all_in) {
  dump_all_ = dump_all_in;
  linearOpTester_.dump_all(dump_all_in);
}

template <class Scalar>
inline bool MultiVectorTester<Scalar>::dump_all() const {
  return dump_all_;
}

}  // namespace Thyra

#endif  // THYRA_MULTI_VECTOR_TESTER_DECL_HPP
