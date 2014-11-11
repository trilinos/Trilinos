// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
 *  \brief A Slice defines a subset of a container.
 */

// Teuchos includes
#include "Teuchos_Assert.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"

namespace Domi
{

/**
 * \brief A Slice contains a start, stop, and step index, describing a
 *        subset of an ordered container.
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
 *
 *   \code
 *   typedef Domi::Slice slc;
 *   \endcode
 *
 * Thus if you wanted a view into an existing <tt>MDArray</tt> (named
 * <tt>array</tt>) consisting of the 5th through 10th elements of that
 * array, you could specify
 *
 *   \code
 *   view = array[slc(5,10)];
 *   \endcode
 *
 * and you would obtain an <tt>MDArrayView</tt>.  Note that the
 * meaning of the stop index is up to the container class, not the
 * Slice struct.  It is recommended that the stop index be
 * non-inclusive.  Thus in the example above, the length of view would
 * be 10 - 5 = 5 and would support integer indexes 0, 1, ... 4.
 *
 * There are several Slice constructors, taking 0-3 arguments and
 * supporting both implicit and explicit default values.  For
 * explicitly requesting a default value, you can use Slice::Default.
 * It is also recommended that the user define a shorthand for this as
 * well.  For example,
 *
 *   \code
 *   const Domi::dim_type & dflt = Domi::Slice::Default;
 *   \endcode
 *
 * For the common case of a positive step value, the default start
 * index is zero and the default stop index is the size of the
 * dimension being indexed.  If the step value is negative, these
 * defaults are reversed.  The default step value is one.  Given these
 * rules, plus the 3-argument constructor meaning <tt>Slice(start,
 * stop, step)</tt>, the following constructors have the following
 * meanings:
 *
 *   <tt>Slice()</tt> is equivalent to <tt>Slice(0,Default,1)</tt>.
 *
 *   <tt>Slice(3)</tt> is equivalent to <tt>Slice(0,3,1)</tt>.
 *
 *   <tt>Slice(1,4)</tt> is equivalent to <tt>Slice(1,4,1)</tt>.
 *
 *   <tt>Slice(Default,5,2)</tt> is equivalent to <tt>Slice(0,5,2)</tt>.
 *
 *   <tt>Slice(Default,Default,-1)</tt> is equivalent to
 *   <tt>Slice(Default,0,-1)</tt>.
 *
 * Note again that it is up to the container class to recognize that a
 * value of <tt>Default</tt> refers to the size of the container, and
 * not to the literal value of <tt>Default</tt> (which would be set,
 * for practical purposes, to the maximum value supported by the
 * dim_type, which might be something like 2**31 as an example).
 *
 * A container class can easily convert a Slice object that
 * (potentially) has default values to a Slice object that has
 * concrete values by calling the <tt>bounds()</tt> method, which
 * takes as its single argument the size of the container.
 *
 */
struct Slice
{
public:

  /** \brief The type for start indexes, stop indexes, and step
   *         intervals
   */
  typedef Domi::dim_type dim_type;

  /** \brief Default value for Slice constructors
   *
   *  The <tt>Default</tt> value can be provided to Slice constructors
   *  to set start, stop, and/or step to default values.  The default
   *  <tt>start</tt> index is zero.  The default <tt>stop</tt> index
   *  is the maximum index of the container the slice ultimately
   *  references.  These default values are reversed if the
   *  <tt>step</tt> is negative.  The default <tt>step</tt> is one.
   */
  static const dim_type Default;

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
  inline Slice(dim_type stopVal);

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
  inline Slice(dim_type startVal, dim_type stopVal, dim_type stepVal=1);

  /** \brief Copy constructor
   */
  inline Slice(const Slice & source);

  /** \brief Destructor
   */
  virtual ~Slice() { }

  //@}

  /** \name Accessor methods
   */
  //@{

  /** \brief Start index
   *
   * If <tt>start</tt> is a non-negative ordinal, then it is a
   * concrete value.  If <tt>start</tt> is negative, it is interpreted
   * to represent <tt>size-start</tt>, where <tt>size</tt> is the size
   * of the container.  If <tt>start</tt> equals <tt>Default</tt>, and
   * <tt>step</tt> is positive, then it is set to zero.  If
   * <tt>start</tt> equals <tt>Default</tt>, and <tt>step</tt> is
   * negative, then it is interpreted to represent the size of the
   * container.
   */
  inline const dim_type start() const;

  /** \brief Stop index
   *
   * If <tt>stop</tt> is a non-negative ordinal, then it is a concrete
   * value.  If <tt>stop</tt> is negative, it is interpreted to
   * represent <tt>size-stop</tt>, where <tt>size</tt> is the size of
   * the container.  If <tt>stop</tt> equals <tt>Default</tt>, and
   * <tt>step</tt> is positive, then it is interpreted to represent
   * the size of the container.  If <tt>start</tt> equals
   * <tt>Default</tt>, and <tt>step</tt> is negative, then it is set
   * to zero.
   */
  inline const dim_type stop() const;

  /** \brief Step interval
   *
   * If <tt>step</tt> is a non-zero ordinal, then it is a concrete
   * value.  If <tt>step</tt> equals <tt>0</tt>, the constructor will
   * throw an exception.  If <tt>step</tt> equals <tt>Default</tt>,
   * the constructor will convert it to a value of one.
   */
  inline const dim_type step() const;

  //@}

  /** \name Slice operators */
  //@{

  /** \brief Assignment operator
   */
  inline Slice & operator=(const Slice & source);

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
   *   \code
   *   Slice bounds = s.bounds(size);
   *   for (Domi::dim_type i = bounds.start(); i != bounds.stop();
   *        i += bounds.step())
   *   { ... }
   *   \endcode
   *
   * Note that in order to accommodate both positive and negative
   * <tt>step</tt>s, the <tt>for</tt> loop continue condition is
   * <tt>(i != bounds.stop())</tt>.  This requires that
   * <tt>bounds()</tt> return precisely the first ordinal outside the
   * bounds that will be returned by <tt>(i += bounds.step())</tt>.
   */
  virtual Slice bounds(dim_type size) const;

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

  // Start index
  dim_type _start;

  // Stop index
  dim_type _stop;

  // Step interval
  dim_type _step;

  // Boolean flag indicating whether the step is positive and the stop
  // index is concrete (i.e. it is not Default and not negative)
  bool _bounded_pos;

  // Boolean flag indicating whether the step is negative and the
  // start index is concrete (i.e. it is not Default and not negative)
  bool _bounded_neg;

};

/**
 * \brief A ConcreteSlice is a Slice that does not accept Default or
 * negative start or stop values.
 *
 * By ensuring that the start, stop and step values are all concrete,
 * the bounds() method can just return *this, which is much more
 * efficient than the base class Slice implementation of bounds().
 *
 * It is not expected that the end user should ever have to deal
 * directly with ConcreteSlice.  The Slice bounds() method returns a
 * Slice object that is in fact a ConcreteSlice.  Any further calls to
 * bounds() will therefore be efficient.
 */
struct ConcreteSlice : public Slice
{
public:
  /** \brief One argument constructor
   *
   *  \param stopVal [in] The stop index of the Slice
   *
   * Returns Slice with <tt>start == 0</tt>, <tt>step == stopVal</tt>,
   * <tt>step == 1</tt>.
   */
  ConcreteSlice(dim_type stopVal);

  /** \brief Two or three argument constructor
   *
   * \param startVal [in] The start index of the slice.
   *
   * \param stopVal [in] The stop index of the slice.
   *
   * \param stepVal [in] The step interval of the slice.
   *
   * Returns Slice with <tt>start == startVal</tt>, <tt>step ==
   * stopVal</tt>, <tt>step == stepVal</tt> (default 1).
   */
  ConcreteSlice(dim_type startVal, dim_type stopVal, dim_type stepVal=1);

  /** \brief Destructor
   */
  virtual ~ConcreteSlice() { }

  /** \brief Simply return this ConcreteSlice
   */
  inline Slice bounds(dim_type size) const { return *this; }

private:

  // Private and not implemented
  ConcreteSlice();
};

/////////////////////////
// Inline implementations
/////////////////////////

Slice::Slice() :
  _start(0),
  _stop(Slice::Default),
  _step(1),
  _bounded_pos(false),
  _bounded_neg(false)
{
}

////////////////////////////////////////////////////////////////////////

Slice::Slice(dim_type stopVal) :
  _start(0),
  _stop(stopVal),
  _step(1),
  _bounded_pos((_stop >= 0 && _stop != Slice::Default)),
  _bounded_neg(false)
{
}

////////////////////////////////////////////////////////////////////////

Slice::Slice(dim_type startVal, dim_type stopVal, dim_type stepVal) :
  _start(((startVal==Slice::Default) && (stepVal > 0)) ? 0 : startVal),
  _stop( ((stopVal ==Slice::Default) && (stepVal < 0)) ? 0 : stopVal ),
  _step(  (stepVal ==Slice::Default) ? 1 : stepVal),
  _bounded_pos(((_step > 0) && (_stop  >= 0 && _stop  != Slice::Default))),
  _bounded_neg(((_step < 0) && (_start >= 0 && _start != Slice::Default)))
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_step == 0),
    InvalidArgument, "Slice step interval cannot be zero"
    );
}

////////////////////////////////////////////////////////////////////////

Slice::Slice(const Slice & source) :
  _start(source._start),
  _stop(source._stop),
  _step(source._step),
  _bounded_pos(source._bounded_pos),
  _bounded_neg(source._bounded_neg)
{
}

////////////////////////////////////////////////////////////////////////

const dim_type Slice::start() const
{
  return _start;
}

////////////////////////////////////////////////////////////////////////

const dim_type Slice::stop() const
{
  return _stop;
}

////////////////////////////////////////////////////////////////////////

const dim_type Slice::step() const
{
  return _step;
}

////////////////////////////////////////////////////////////////////////

Slice & Slice::operator=(const Slice & slice)
{
  _start       = slice._start      ;
  _stop        = slice._stop       ;
  _step        = slice._step       ;
  _bounded_pos = slice._bounded_pos;
  _bounded_neg = slice._bounded_neg;
  return *this;
}

////////////////////////////////////////////////////////////////////////

bool Slice::operator==(const Slice & slice) const
{
  return ((_start == slice._start) &&
          (_stop  == slice._stop ) &&
          (_step  == slice._step )    );
}

////////////////////////////////////////////////////////////////////////

bool Slice::operator!=(const Slice & slice) const
{
  return (not operator==(slice));
}

}  // namespace Domi

#endif
