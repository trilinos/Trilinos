// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Range1D class used for representing a range of positive integers.
// Its primary usage is in accessing vectors and matrices by subregions
// of rows and columns
//

#ifndef TEUCHOS_RANGE1D_HPP
#define TEUCHOS_RANGE1D_HPP

/*! \file Teuchos_Range1D.hpp
    \brief .
*/

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Assert.hpp"


namespace Teuchos {


/** \brief Subregion Index Range Class.
 *
 * The class <tt>%Range1D</tt> encapsulates a 1-D, zero-based, range of
 * non-negative indexes.  It is used to index into vectors and matrices and
 * return subregions of them respectively.
 *
 * Constructing using <tt>Range1D()</tt> yields a range that represents the
 * entire dimension of an object <tt>[0, max_ubound]</tt> (an entire
 * std::vector, all the rows in a matrix, or all the columns in a matrix
 * etc.).
 *
 * Constructing using <tt>\ref Range1D::Range1D "Range1D(INVALID)"</tt> yields
 * an invalid range <tt>[0,-2]</tt> with <tt>size() == -1</tt>.  Once
 * constructed with <tt>Range1D(INVALID)</tt>, a <tt>%Range1D</tt> object can
 * pass through many other operations that may change <tt>%lbound()</tt> and
 * <tt>%ubound()</tt> but will never change <tt>size()==-1</tt>.
 *
 * Constructing using <tt>\ref Range1D::Range1D "Range1D(lbound,ubound)"</tt>
 * yields a finite-dimensional zero-based range.  The validity of constructed
 * range will only be checked if <tt>TEUCHOS_DEBUG</tt> is defined.
 *
 * There are many \ref Range1D_funcs_grp "non-member functions" that can be
 * used with <tt>%Range1D</tt> objects.
 *
 * The default copy constructor and assignment operator functions are allowed
 * since they have the correct semantics.
 */
class Range1D {
public:

  /** \brief Deprecated. */
  typedef Teuchos_Ordinal  Index;

  /** \brief Deprecated. */
  typedef Teuchos_Ordinal  Ordinal;

  /** \brief . */
  enum EInvalidRange { INVALID };

  /** \brief Used for Range1D(INVALID) */
  static const Range1D Invalid;

  /** \brief Construct a full range.
   *
   * Postconditions: <ul>
   * <li> <tt>this->full_range()==true</tt>
   * <li> <tt>this->size()</tt> is a very large number
   * <li> <tt>this->lbound()==0</tt>
   * <li> <tt>this->ubound()</tt> is a very large number
   * </ul>
   */
  inline Range1D();

  /** \brief Constructs an invalid range.
   *
   * Postconditions: <ul>
   * <li> <tt>this->full_range() == false</tt>
   * <li> <tt>this->size() == -1</tt>
   * <li> <tt>this->lbound()==0</tt>
   * <li> <tt>this->ubound()==-2</tt>
   * </ul>
   */
  inline Range1D(EInvalidRange);

  /** \brief Construct a finite range <tt>[lbound, ubound]</tt>.
   *
   * Preconditions: <ul>
   * <li> <tt>lbound >= 0</tt> (throw \c out_of_range)
   * <li> <tt>ubound >= lbound-1</tt> (throw \c out_of_range)
   * </ul>
   *
   * Postconditions: <ul>
   * <li> <tt>this->full_range() == false</tt>
   * <li> <tt>this->size() == ubound - lbound + 1</tt>
   * <li> <tt>this->lbound() == lbound</tt>
   * <li> <tt>this->ubound() == ubound</tt>
   * </ul>
   *
   * \note It is allowed for <tt>ubound == lbound-1</tt> which yields a
   * zero-sized range.  There are use cases where this is useful so it is
   * allowed.
   */
  inline Range1D(Ordinal lbound, Ordinal ubound);

  /** \brief Returns \c true if the range represents the entire region. */
  inline bool full_range() const;

  /** \brief Return lower bound of the range */
  inline Ordinal lbound() const;

  /** \brief Return upper bound of the range */
  inline Ordinal ubound() const;

  /** \brief Return the size of the range (<tt>ubound() - lbound() + 1</tt>) */
  inline Ordinal size() const;

  /** \brief Return true if the index is in range */
  inline bool in_range(Ordinal i) const;

  /** \brief Increment the range by a constant
   *
   * \precondition <tt>this->lbound() + incr >= 0</tt> (throws \c out_of_range)
   */
  inline Range1D& operator+=( Ordinal incr );

  /** \brief Deincrement the range by a constant.
   *
   * \precondition <tt>this->lbound() - incr >= 0</tt> (throws \c out_of_range)
   */
  inline Range1D& operator-=( Ordinal incr );

private:

  Ordinal lbound_;
  Ordinal ubound_;

  inline void assert_valid_range(Ordinal lbound, Ordinal ubound) const;

}; // end class Range1D


/** \brief rng1 == rng2.
 *
 * @return Returns <tt>rng1.lbound() == rng2.ubound() && rng1.ubound() == rng2.ubound()</tt>.
 *
 * \relates Range1D
 */
inline bool operator==(const Range1D& rng1, const Range1D& rng2 )
{
  return rng1.lbound() == rng2.lbound() && rng1.ubound() == rng2.ubound();
}


/** \brief rng1 == rng2.
 *
 * @return Returns <tt>rng1.lbound() == rng2.ubound() && rng1.ubound() == rng2.ubound()</tt>.
 *
 * \relates Range1D
 */
inline bool operator!=(const Range1D& rng1, const Range1D& rng2 )
{
  return !(rng1 == rng2);
}


/** \brief rng_lhs = rng_rhs + i.
  *
  * Increments the upper and lower bounds by a constant.
  *
  * Postcondition: <ul>
  *	<li> <tt>rng_lhs.lbound() == rng_rhs.lbound() + i</tt>
  *	<li> <tt>rng_lhs.ubound() == rng_rhs.ubound() + i</tt>
  *	</ul>
 *
 * \relates Range1D
  */
inline Range1D operator+(const Range1D &rng_rhs, Range1D::Ordinal i)
{
    return Range1D(i+rng_rhs.lbound(), i+rng_rhs.ubound());
}


/** \brief rng_lhs = i + rng_rhs.
  *
  * Increments the upper and lower bounds by a constant.
  *
  * Postcondition: <ul>
  *	<li> <tt>rng_lhs.lbound() == i + rng_rhs.lbound()</tt>
  *	<li> <tt>rng_lhs.ubound() == i + rng_rhs.ubound()</tt>
  *	</ul>
 *
 * \relates Range1D
  */
inline Range1D operator+(Range1D::Ordinal i, const Range1D &rng_rhs)
{
    return Range1D(i+rng_rhs.lbound(), i+rng_rhs.ubound());
}


/** \brief rng_lhs = rng_rhs - i.
  *
  * Deincrements the upper and lower bounds by a constant.
  *
  * Postcondition: <ul>
  *	<li> <tt>rng_lhs.lbound() == rng_rhs.lbound() - i</tt>
  *	<li> <tt>rng_lhs.ubound() == rng_rhs.ubound() - i</tt>
  *	</ul>
 *
 * \relates Range1D
  */
inline Range1D operator-(const Range1D &rng_rhs, Range1D::Ordinal i)
{
    return Range1D(rng_rhs.lbound()-i, rng_rhs.ubound()-i);
}


/** \brief Return a bounded index range from a potentially unbounded index
  * range.
  *
  * Return a index range of lbound to ubound if rng.full_range() == true
  * , otherwise just return a copy of rng.
  *
  * Postconditions: <ul>
  *	<li> [<tt>rng.full_range() == true</tt>] <tt>return.lbound() == lbound</tt>
  *	<li> [<tt>rng.full_range() == true</tt>] <tt>return.ubound() == ubound</tt>
  *	<li> [<tt>rng.full_range() == false</tt>] <tt>return.lbound() == rng.lbound()</tt>
  *	<li> [<tt>rng.full_range() == false</tt>] <tt>return.ubound() == rng.ubound()</tt>
  *	</ul>
 *
 * \relates Range1D
  */
inline Range1D full_range(const Range1D &rng, Range1D::Ordinal lbound, Range1D::Ordinal ubound)
{	return rng.full_range() ? Range1D(lbound,ubound) : rng; }


/** \brief Print out to ostream.
 *
 * \relates Range1D
 */
TEUCHOSCORE_LIB_DLL_EXPORT
std::ostream& operator<<(std::ostream &out, const Range1D& rng);


// //////////////////////////////////////////////////////////
// Inline members

inline
Range1D::Range1D()
  : lbound_(0), ubound_(std::numeric_limits<Ordinal>::max()-1)
{}

inline
Range1D::Range1D( EInvalidRange )
  : lbound_(0), ubound_(-2)
{}


inline
Range1D::Range1D(Ordinal lbound_in, Ordinal ubound_in)
  : lbound_(lbound_in), ubound_(ubound_in)
{
  assert_valid_range(lbound_in,ubound_in);
}

inline
bool Range1D::full_range() const {
  return (lbound_ == 0 && ubound_ == std::numeric_limits<Ordinal>::max()-1);
}

inline
Range1D::Ordinal Range1D::lbound() const {
  return lbound_;
}

inline
Range1D::Ordinal Range1D::ubound() const {
  return ubound_;
}

inline
Range1D::Ordinal Range1D::size() const {
  return ubound_ - lbound_ + 1;
}

inline
bool Range1D::in_range(Ordinal i) const {
  return lbound_ <= i && i <= ubound_;
}

inline
Range1D& Range1D::operator+=( Ordinal incr ) {
  assert_valid_range( lbound_ + incr, ubound_ + incr );
  lbound_ += incr;
  ubound_ += incr;
  return *this;
}

inline
Range1D& Range1D::operator-=( Ordinal incr )
{
  assert_valid_range( lbound_ - incr, ubound_ - incr );
  lbound_ -= incr;
  ubound_ -= incr;
  return *this;
}


// See Range1D.cpp
inline
void Range1D::assert_valid_range(Ordinal lbound_in, Ordinal ubound_in) const
{
  (void)lbound_in; (void)ubound_in;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_INEQUALITY(lbound_in, >=, 0);
  TEUCHOS_ASSERT_INEQUALITY(ubound_in, >=, lbound_in - 1);
#endif
}

} // end namespace Teuchos

#endif // end TEUCHOS_RANGE1D_HPP
