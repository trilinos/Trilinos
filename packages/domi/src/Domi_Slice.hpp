// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_SLICE_HPP
#define DOMI_SLICE_HPP

/** \file Domi_Slice.hpp
 *
 *  \brief A SLice defines a subset of a container.
 */

#include "Domi_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"

namespace Domi
{

/**
 * \brief A Slice contains a start, stop, and step index, describing a
 *        subset of a container.
 *
 * A Slice is an alternative to an index and specifies a range of
 * values with a start value, a stop value, and a step interval.
 * Containers that support indexing with a Slice (such as the
 * <tt>MDArray</tt> classes) typically overload the
 * <tt>operator[]</tt> method and return a view into the same type of
 * container being referenced.
 *
 * The formal constructors for Slices are of the form
 * <tt>Domi::Slice(...)</tt>. To ease their use as index-like objects,
 * it is recommended that the user define a typedef to provide a
 * shorthand:
 
 \code
 typedef Domi::Slice slc;
 \endcode

 * Thus if you wanted a view into an existing <tt>MDArray</tt> (named
 * <tt>array</tt>) consisting of the 5th through 10th elements of that
 * array, you could specify

 \code
 view = array[slc(5,10)];
 \endcode

 * and you would obtain an <tt>MDArrayView</tt>.  Note that the
 * meaning of the stop index is up to the container class, not the
 * Slice struct.  It is recommended that the stop index be
 * non-inclusive.  Thus in the example above, the length of view would
 * be 10 - 5 = 5 and would support integer indexes 0, 1, ... 4.

 * There are several Slice constructors, taking 0-3 arguments and
 * supporting both implicit and explicit default values.  For
 * explicitly requesting a default value, you can use Slice::Default.
 * It is also recommended that the user define a shorthand for this as
 * well.  For example,

 \code
 const Domi_Ordinal & dflt = Domi::Slice::Default;
 \endcode

 * For the common case of a positive step value, the default start
 * index is zero and the default stop index is the size of the
 * dimension being indexed.  If the step value is negative, these
 * defaults are reversed.  The default step value is one.  Given these
 * rules, plus the 3-argument constructor meaning <tt>Slice(start,
 * stop, step)</tt>, the following constructors have the following
 * meanings:

 * <tt>Slice()</tt> is equivalent to <tt>Slice(0,Default,1)</tt>.

 * <tt>Slice(3)</tt> is equivalent to <tt>Slice(0,3,1)</tt>.

 * <tt>Slice(1,4)</tt> is equivalent to <tt>Slice(1,4,1)</tt>.

 * <tt>Slice(Default,5,2)</tt> is equivalent to <tt>Slice(0,5,2)</tt>.

 * <tt>Slice(Default,Default,-1)</tt> is equivalent to
 * <tt>Slice(Default,0,-1)</tt>.

 * Note again that it is up to the container class to recognize that a
 * value of <tt>Default</tt> refers to the size of the container, and
 * not to the literal value of <tt>Default</tt> (which would be set,
 * for practical purposes, to the maximum value supported by the
 * Ordinal, which might be something like 2**63 as an example).

 * A container class can easily convert a Slice object that
 * (potentially) has default values to a Slice object that has
 * concrete values by calling the <tt>bounds()</tt> method, which
 * takes as its single argument the size of the container.

 */
struct Slice
{
public:

  /** \brief The type for start indexes, stop indexes, and step
      intervals
   */
  typedef Teuchos_Ordinal Ordinal;

  /** \brief Default value for Slice constructors
   *
   *  The <tt>Default</tt> value can be provided to Slice constructors
   *  to set start, stop, and/or step to default values.  The default
   *  <tt>start</tt> index is zero.  The default <tt>stop</tt> index
   *  is the maximum index of the container the slice ultimately
   *  references.  These default values are reversed if the
   *  <tt>step</tt> is negative.  The default <tt>step</tt> is one.
   */
  static const Ordinal Default;

  /** \name Public data members
   *
   * Note that these public data members are const, can only be set
   * via a constructor, and that the Slice struct is <i>immutable</i>.
   */
  //@{

  /** \brief Start index
   *
   * If <tt>start</tt> is a non-negative ordinal, then it is a
   * concrete value.  If <tt>start</tt> is negative, it is interpreted
   * to represent <tt>size-start</tt>, where <tt>size</tt> is the size
   * of the container.  If <tt>start</tt> equals <tt>Default</tt>,
   * then it is interpreted to represent the size of the container.
   */
  const Ordinal start;

  /** \brief Stop index
   *
   * If <tt>stop</tt> is a non-negative ordinal, then it is a concrete
   * value.  If <tt>stop</tt> is negative, it is interpreted to
   * represent <tt>size-stop</tt>, where <tt>size</tt> is the size of
   * the container.  If <tt>stop</tt> equals <tt>Default</tt>, then it
   * is interpreted to represent the size of the container.
   */
  const Ordinal stop;

  /** \brief Step interval
   *
   * If <tt>step</tt> is a non-zero ordinal, then is is a concrete
   * value.  If <tt>step</tt> equals <tt>0</tt>, the constructor will
   * throw an exception.  If <tt>step</tt> equals <tt>Default</tt>,
   * the constructor will convert it to a value of one.
   */
  const Ordinal step;

  //@}

  /** \name Constructors and destructor */

  //@{

  /** \brief Default constructor
   *
   * Returns <tt>start == 0</tt>, <tt>step == Default</tt>, <tt>step == 1</tt>.
   */
  inline Slice();

  /** \brief One argument constructor
   *
   * \param stopVal [in] The stop index of the slice.
   *
   * Returns Slice with <tt>start == 0</tt>, <tt>step == stopVal</tt>,
   * <tt>step == 1</tt>.
   */
  inline Slice(Ordinal stopVal);

  /** \brief Two or three argument constructor
   *
   * \param startVal [in] The start index of the slice.
   *
   * \param stopVal [in] The stop index of the slice.
   *
   * \param stepVal [in] The step interval of the slice.
   *
   * Returns Slice with <tt>start == startVal</tt>, <tt>step ==
   * stopVal</tt>, <tt>step == stepVal</tt> (default 1).  If
   * <tt>stepVal</tt> is positive, then <tt>startVal == Default</tt>
   * is converted to <tt>start = 0</tt>.  If <tt>stepVal</tt> is
   * negative, then <tt>stopVal == Default</tt> is converted to
   * <tt>stop = 0</tt>.  If <tt>stepVal</tt> is zero, an exception is
   * thrown.
   */
  inline Slice(Ordinal startVal, Ordinal stopVal, Ordinal stepVal=1);

  /** \brief Destructor
   */
  ~Slice() { }

  //@}

  /** \name Slice operators */
  //@{

  /** \brief Equals operator
   */
  inline bool operator==(const Slice & slice) const;

  /** \brief Inequality operator
   */
  inline bool operator!=(const Slice & slice) const;

  //@}

  /** \brief Return a <tt>Slice</tt> with concrete <tt>start</tt> and
   * <tt>stop</tt> values
   *
   * \param size [in] The size of the container
   *
   * The <tt>Slice</tt> returned by the <tt>bounds()</tt> method will
   * be safe and accurate for a subsequent <tt>for</tt> loop.  For
   * example, if <tt>s</tt> is a <tt>Slice</tt> that may contain
   * <tt>Default</tt> values or negative indexes (representing
   * distance from the end of the container), and <tt>size</tt> is the
   * size of a container, the following code is valid:
   *
   * \code
   * Slice bounds = s.bounds(size);
   * for (Teuchos_Ordinal i = bounds.start; i != bounds.stop; i += bounds.step)
   * { ... }
   * \endcode
   *
   * Note that in order to accommodate both positive and negative
   * <tt>step</tt>s, the <tt>for</tt> loop continue condition is <tt>(i
   * != bounds.stop)</tt>.  This requires that <tt>bounds()</tt> return
   * precisely the first ordinal outside the bounds that will be
   * returned by <tt>(i += bounds.step)</tt>.
   */
  Slice bounds(Ordinal size) const;

  /** \brief Return a string representation of the Slice
   *
   * Return a string, encapsulated by "[" and "]", with the start,
   * stop, and step values separated by colons (":").  If the start or
   * stop values are default values, represent them with the null
   * string.  If the step value is 1, do not include the step or the
   * second colon.
   */
  std::string toString() const;

  /** \brief Streaming output operator
   *
   * Friend operator for streaming output
   */
  friend std::ostream & operator<<(std::ostream & os,
                                   const Slice & slice);

private:

  // Boolean flag indicating whether the step is positive and the stop
  // index is concrete (i.e. it is not Default and not negative)
  const bool _bounded_pos;

  // Boolean flag indicating whether the step is negative and the
  // start index is concrete (i.e. it is not Default and not negative)
  const bool _bounded_neg;

  // This is private and never implemented, thus eliminating it from
  // the interface.  Since Slices are immutable, the assignment
  // operator is of little value.
  Slice operator=(const Slice & source);

};

/////////////////////////
// Inline implementations
/////////////////////////

Slice::Slice() :
  start(0),
  stop(Slice::Default),
  step(1),
  _bounded_pos(false),
  _bounded_neg(false)
{
}

Slice::Slice(Ordinal stopVal) :
  start(0),
  stop(stopVal),
  step(1),
  _bounded_pos((stop >= 0 && stop != Slice::Default)),
  _bounded_neg(false)
{
}

Slice::Slice(Ordinal startVal, Ordinal stopVal, Ordinal stepVal) :
  start(((startVal==Slice::Default) && (stepVal > 0)) ? 0 : startVal),
  stop( ((stopVal ==Slice::Default) && (stepVal < 0)) ? 0 : stopVal ),
  step(  (stepVal ==Slice::Default) ? 1 : stepVal),
  _bounded_pos(((step > 0) && (stop  >= 0 && stop  != Slice::Default))),
  _bounded_neg(((step < 0) && (start >= 0 && start != Slice::Default)))
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (step==0), std::invalid_argument, "Slice step interval cannot be zero"
    );
}

bool Slice::operator==(const Slice & slice) const
{
  return ((start == slice.start) &&
          (stop  == slice.stop ) &&
          (step  == slice.step )    );
}

bool Slice::operator!=(const Slice & slice) const
{
  return (not operator==(slice));
}

}  // namespace Domi

#endif
