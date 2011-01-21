// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_SPARSEARRAY_HPP
#define STOKHOS_SPARSEARRAY_HPP

#include <algorithm>
#include "Teuchos_Array.hpp"

namespace Stokhos {

  template <typename index_iterator, typename value_iterator> 
  struct SparseArrayIterator;

  template <typename index_iterator, typename value_iterator> 
  struct SparseArrayReverseIterator;
    
  //! Container for a "sparse" array
  template <typename ordinal_type, typename val_type>
  struct SparseArray {

    //! Indices of nonzeros
    Teuchos::Array<ordinal_type> indices;

    //! Nonzero values
    Teuchos::Array<val_type> values;
    
    //! Return size of array
    ordinal_type size() const { return indices.size(); }
    
    //! Resize array
    void resize(ordinal_type sz) {
      indices.resize(sz);
      values.resize(sz);
    }
    
    typedef typename Teuchos::Array<ordinal_type>::const_iterator index_const_iterator;
    typedef typename Teuchos::Array<val_type>::const_iterator value_const_iterator;
    typedef SparseArrayIterator<index_const_iterator, value_const_iterator> const_iterator;
    typedef SparseArrayReverseIterator<index_const_iterator, value_const_iterator> const_reverse_iterator;
    
    //! Iterator pointing to beginning of array
    const_iterator begin() const { 
      return const_iterator(indices.begin(), values.begin()); 
    }

    //! Iterator pointing to end of array
    const_iterator end() const { 
      return const_iterator(indices.end(), values.end()); 
    }

    //! Reverse iterator pointing to end of array
    const_reverse_iterator rbegin() const { 
      return const_reverse_iterator(end()); 
    }

    //! Reverse iterator pointing to begining of array
    const_reverse_iterator rend() const { 
      return const_reverse_iterator(begin()); 
    }

    //! Return iterator pointing to given index \c i
    /*!
     * Returns \c end() if array doesn't contain index \c i
     */
    const_iterator find(ordinal_type i) const {
      const_iterator it = std::lower_bound(begin(), end(), i);
      bool found = it != end() && *it == i;
      if (found)
	return it;
      else
	return end();
    }

  };

  //! Bi-directional iterator for traversing a sparse array
  /*!
   * The "value" of the iterator is the index, which is what you get from
   * dereferencing the iterator (* or ->).  There is also a method called
   * value() that gives you the value of the sparse array pointed at by the
   * iterator.
   *
   * This could easily be a random access iterator.  I just haven't 
   * implemented those methods.
   */
  template <typename index_iterator_type, typename value_iterator_type>
  class SparseArrayIterator : 
    public std::iterator<
      std::bidirectional_iterator_tag, 
      typename std::iterator_traits<index_iterator_type>::value_type,
      typename std::iterator_traits<index_iterator_type>::difference_type,
      typename std::iterator_traits<index_iterator_type>::pointer,
      typename std::iterator_traits<index_iterator_type>::reference > {
  public:

    typedef std::iterator<
      std::bidirectional_iterator_tag, 
      typename std::iterator_traits<index_iterator_type>::value_type,
      typename std::iterator_traits<index_iterator_type>::difference_type,
      typename std::iterator_traits<index_iterator_type>::pointer,
      typename std::iterator_traits<index_iterator_type>::reference > base_type;
    typedef typename base_type::iterator_category iterator_category;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;
    typedef typename std::iterator_traits<value_iterator_type>::reference value_reference;

    //! Default constructor
    SparseArrayIterator() : index_iterator(), value_iterator() {}
    
    //! Constructor
    SparseArrayIterator(index_iterator_type index_it, 
			value_iterator_type value_it) :
      index_iterator(index_it), value_iterator(value_it) {}
    
    //! == operator
    bool operator==(const SparseArrayIterator& it) const {
      return index_iterator == it.index_iterator && 
	value_iterator == it.value_iterator;
    }
    
    //! != operator
    bool operator!=(const SparseArrayIterator& it) const {
      return index_iterator != it.index_iterator || 
	value_iterator != it.value_iterator;
    }
    
    //! * operator 
    reference
    operator*() const {
      return *index_iterator;
    }
    
    //! -> operator
    pointer
    operator->() const {
      return index_iterator.operator->();
    }
    
    //! Prefix ++
    SparseArrayIterator& operator++() {
	index_iterator++;
	value_iterator++;
      return *this;
    }
    
    //! Postfix ++
    SparseArrayIterator operator++(int) {
      SparseArrayIterator tmp(*this);
      ++(*this);
      return tmp;
    }

    //! Prefix --
    SparseArrayIterator& operator--() {
	index_iterator--;
	value_iterator--;
      return *this;
    }
    
    //! Postfix --
    SparseArrayIterator operator--(int) {
      SparseArrayIterator tmp(*this);
      --(*this);
      return tmp;
    }

    //! Return true of *this < b
    bool operator<(const SparseArrayIterator& b) const {
      return index_iterator < b.index_iterator && 
	value_iterator < b.value_iterator;
    }
    
    //! Return value associated with iterator
    value_reference value() const {
      return *value_iterator;
    }
    
  protected:
    
    //! Index iterator
    index_iterator_type index_iterator;

    //! Value iterator
    value_iterator_type value_iterator;
    
  };

  //! Bi-directional reverse iterator for traversing a sparse array
  template <typename index_iterator_type, typename value_iterator_type>
  class SparseArrayReverseIterator : 
    public std::reverse_iterator< SparseArrayIterator<index_iterator_type,
						      value_iterator_type> > {
  public:

    typedef SparseArrayIterator<index_iterator_type, value_iterator_type> iterator_type;
    typedef std::reverse_iterator<iterator_type> base_type;
    typedef typename base_type::iterator_category iterator_category;
    typedef typename base_type::value_type value_type;
    typedef typename base_type::difference_type difference_type;
    typedef typename base_type::reference reference;
    typedef typename base_type::pointer pointer;
    typedef typename std::iterator_traits<value_iterator_type>::reference value_reference;

    //! Default constructor
    SparseArrayReverseIterator() : base_type() {}
    
    //! Constructor
    SparseArrayReverseIterator(iterator_type it) : base_type(it) {}
    
    //! Return value associated with iterator
    value_reference value() const {
      iterator_type tmp = this->base();
      --tmp;
      return tmp.value();
    }
    
  };

} // namespace Stokhos

#endif // STOKHOS_SPARSEARRAY_HPP
