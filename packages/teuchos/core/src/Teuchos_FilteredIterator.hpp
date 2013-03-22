// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_FILTERED_ITERATOR_HPP
#define TEUCHOS_FILTERED_ITERATOR_HPP


#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Exceptions.hpp"


namespace Teuchos {


/** \brief C++ Standard Library compatable filtered iterator.
 *
 * Provides a bidirectional iterator providing a filtered view of an
 * underlying iterator (defined by a predicate policy object).
 */
template<class IteratorType, class Predicate>
class FilteredIterator {
public:

  /** \name public typedefs */
  //@{

  /** \brief . */
  typedef	std::bidirectional_iterator_tag iterator_category;
  /** \brief . */
  typedef	typename std::iterator_traits<IteratorType>::value_type	value_type;
  /** \brief . */
  typedef typename std::iterator_traits<IteratorType>::reference reference;
  /** \brief . */
  typedef typename std::iterator_traits<IteratorType>::pointer pointer;
  /** \brief . */
  typedef	typename std::iterator_traits<IteratorType>::difference_type difference_type;
  
  //@}
    
  /** \name Constructors. */
  //@{

  /** \brief construct to a null iterator. */
  FilteredIterator()
    {}
  /** \brief Construct with iterator and range.
   * \todo Assert valid range for random-access iterators.
   */
  FilteredIterator(IteratorType current_in, IteratorType begin_in, IteratorType end_in,
    Predicate pred_in = Predicate()
    )
    :current_(current_in), begin_(begin_in), end_(end_in), pred_(pred_in)
    { advanceForwardToValid(); }	
  /** \brief Convert type of iterators (mainly for non-const to const). */
  template<class IteratorType2, class Predicate2>
  FilteredIterator(const FilteredIterator<IteratorType2,Predicate2>& rhs)
    :current_(rhs.current()), begin_(rhs.begin()), end_(rhs.end()), pred_(rhs.pred())
    {}
  /** \brief Assign different types of iterators (mainly for non-const to const). */
  template<class IteratorType2, class Predicate2>
  FilteredIterator & operator=(const FilteredIterator<IteratorType2,Predicate2>& rhs)
    {
      current_ = rhs.current();
      begin_ = rhs.begin();
      end_ = rhs.end();
      pred_ = rhs.pred();
      return *this;
    }
  
  //@}

  /** @name Access */
  //@{

  /** \brief itr* */
  reference operator*()	const 
    { return *current_; }
  /** \brief itr->member */
  pointer operator->() const
    { return current_.operator->(); }

  //@}

  /** @name Incrementation */
  //@{

  /** \brief ++itr */
  FilteredIterator& operator++()
    {
      assertNotIterateForwardPastEnd();
      ++current_;
      advanceForwardToValid();
      return *this;
    }
  /** \brief itr++ */
  const FilteredIterator operator++(int)
    {
      FilteredIterator tmp = *this;
      ++*this;
      return tmp;
    }
  /** \brief --itr */
  FilteredIterator&	operator--()
    {
      assertNotIterateBackwardPastBegin();
      --current_;
      advanceBackwardToValid();
      return *this;
    }
  /** \brief itr-- */
  const FilteredIterator operator--(int) 
    {
      FilteredIterator tmp = *this;
      --*this;
      return tmp;
    }

  //@}

  /** @name Implementation access */
  //@{

  /** \brief . */
  IteratorType	current() const { return current_; }
  /** \brief . */
  IteratorType	begin() const { return begin_; }
  /** \brief . */
  IteratorType	end() const { return end_; }
  /** \brief . */
  Predicate pred() const{ return pred_; }

  //@}

private: // Data members

  /** \brief . */
  IteratorType current_;
  /** \brief . */
  IteratorType begin_;
  /** \brief . */
  IteratorType end_;
  /** \brief . */
  Predicate pred_;

private: // Functions

  /** \brief . */
  void advanceForwardToValid();
  /** \brief . */
  void advanceBackwardToValid();
  /** \brief . */
  void assertNotIterateForwardPastEnd()
#ifndef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    {}
#else
  ;
#endif
  /** \brief . */
  void assertNotIterateBackwardPastBegin()
#ifndef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    {}
#else
  ;
#endif

};


/** \brief itr1 == itr2.
 * \relates FilteredIterator
 */
template<class IteratorType, class Predicate>
inline bool operator==(const FilteredIterator<IteratorType,Predicate>& itr1, 
  const FilteredIterator<IteratorType,Predicate>& itr2)
{
  return itr1.current() == itr2.current();
}


/** \brief itr1 != itr2.
 * \relates FilteredIterator
 */
template<class IteratorType, class Predicate>
inline bool operator!=(const FilteredIterator<IteratorType,Predicate>& itr1, 
  const FilteredIterator<IteratorType,Predicate>& itr2)
{
  return itr1.current() != itr2.current();
}


/** \brief ostream operator.
 *
 * WARNING: This requires that IteratorType also support operator<<().
 *
 * \relates FilteredIterator
 */
template<class IteratorType, class Predicate>
std::ostream& operator<<(std::ostream &out, const FilteredIterator<IteratorType,Predicate>& itr)
{
  out << "FilteredIterator{current=???, end=???, pred="<<TypeNameTraits<Predicate>::name()<<"}";
  return out;
}

//
// Template definitions
//


template<class IteratorType, class Predicate>
void FilteredIterator<IteratorType,Predicate>::advanceForwardToValid()
{
  while (current_ != end_ && !pred_(*current_)) {
    ++current_;
  }
}


template<class IteratorType, class Predicate>
void FilteredIterator<IteratorType,Predicate>::advanceBackwardToValid()
{
  while (current_ != begin_ && !pred_(*current_)) {
    --current_;
  }
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


template<class IteratorType, class Predicate>
void FilteredIterator<IteratorType,Predicate>::assertNotIterateForwardPastEnd()
{
  const bool current_is_at_end = (current_ == end_);
  TEUCHOS_TEST_FOR_EXCEPTION( current_is_at_end, RangeError,
    "Error, trying to iterate " << *this << " forward ++ past end!");
}


template<class IteratorType, class Predicate>
void FilteredIterator<IteratorType,Predicate>::assertNotIterateBackwardPastBegin()
{
  const bool current_is_at_begin = (current_ == begin_);
  TEUCHOS_TEST_FOR_EXCEPTION( current_is_at_begin, RangeError,
    "Error, trying to iterate " << *this << " backward -- past begin!");
}


#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK


}	// end namespace Teuchos


#endif // TEUCHOS_FILTERED_ITERATOR_HPP
