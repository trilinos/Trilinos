// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SPARSE3TENSOR_HPP
#define STOKHOS_SPARSE3TENSOR_HPP

#include <ostream>
#include <map>
#include "Stokhos_SparseArray.hpp"

namespace Stokhos {

  /*! 
   * \brief Data structure storing a sparse 3-tensor C(i,j,k) in a 
   * a compressed format.
   */
  template <typename ordinal_type, typename value_type>
  class Sparse3Tensor {

  private:

    typedef std::map<const ordinal_type, value_type> i_map;
    typedef std::map<const ordinal_type, i_map> ji_map;
    typedef std::map<const ordinal_type, ji_map> kji_map;

    typedef std::map<const ordinal_type, value_type> j_map;
    typedef std::map<const ordinal_type, j_map> kj_map;
    typedef std::map<const ordinal_type, kj_map> ikj_map;

    typedef SparseArray<ordinal_type, value_type> i_sparse_array;
    typedef SparseArray<ordinal_type, i_sparse_array> ji_sparse_array;
    typedef SparseArray<ordinal_type, ji_sparse_array> kji_sparse_array;

    typedef SparseArray<ordinal_type, value_type> j_sparse_array;
    typedef SparseArray<ordinal_type, j_sparse_array> kj_sparse_array;
    typedef SparseArray<ordinal_type, kj_sparse_array> ikj_sparse_array;

  public:

    //! Iterator for looping over k entries
    typedef typename kji_sparse_array::const_iterator k_iterator;

    //! Iterator for looping over k entries in reverse
    typedef typename kji_sparse_array::const_reverse_iterator k_reverse_iterator;

    //! Iterator for looping over j entries given k
    typedef typename ji_sparse_array::const_iterator kj_iterator;

    //! Iterator for looping over i entries given k and j
    typedef typename j_sparse_array::const_iterator kji_iterator;

    //! Iterator for looping over i entries
    typedef typename ikj_sparse_array::const_iterator i_iterator;

    //! Iterator for looping over i entries in reverse
    typedef typename ikj_sparse_array::const_reverse_iterator i_reverse_iterator;

    //! Iterator for looping over k entries given i
    typedef typename kj_sparse_array::const_iterator ik_iterator;

    //! Iterator for looping over j entries given i and k
    typedef typename j_sparse_array::const_iterator ikj_iterator;
    
    //! Constructor
    Sparse3Tensor();
    
    //! Destructor
    ~Sparse3Tensor() {}

    //! Add new term for given (i,j,k)
    void add_term(ordinal_type i, ordinal_type j, ordinal_type k,
		  const value_type& c);

    //! Add new term for given (i,j,k) and sum in if already there
    void sum_term(ordinal_type i, ordinal_type j, ordinal_type k,
		  const value_type& c);

    //! Signal all terms have been added
    /*!
     * This method must be called before using any of the accessor methods
     * below.  It sets up a new set of data structures that are much more
     * efficient for iterating through the tensor.
     */
    void fillComplete();

    //! Return whether fillComplete() has been called
    bool fillCompleted() const;

    //! Print tensor
    void print(std::ostream& os) const;

    //! Get Cijk value for a given i, j, k indices
    value_type getValue(ordinal_type i, ordinal_type j, ordinal_type k) const;

    //! Return number of non-zero entries
    ordinal_type num_entries() const;

    /** \name k-based data-structure accessor methods */
    //@{

    //! Number of k entries in C(i,j,k)
    ordinal_type num_k() const;

    //! Number of j entries in C(i,j,k) for given k
    ordinal_type num_j(const k_iterator& k) const;

    //! Number of i entries in C(i,j,k) for given k and j
    ordinal_type num_i(const kj_iterator& j) const;

    //! Return k iterator for given index k
    k_iterator find_k(ordinal_type k) const;

    //! Return j iterator given k iterator and index j
    kj_iterator find_j(const k_iterator& k, ordinal_type j) const;

    //! Return i iterator given j iterator and index i
    kji_iterator find_i(const kj_iterator& j, ordinal_type i) const;

    //! Iterator pointing to first k entry
    k_iterator k_begin() const;

    //! Iterator pointing to last k entry
    k_iterator k_end() const;

    //! Reverse iterator pointing to last k entry
    k_reverse_iterator k_rbegin() const;

    //! Reverse iterator pointing to first k entry
    k_reverse_iterator k_rend() const;

    //! Iterator pointing to first j entry for given k
    kj_iterator j_begin(const k_iterator& k) const;

    //! Iterator pointing to last j entry for given k
    kj_iterator j_end(const k_iterator& k) const;

    //! Iterator pointing to first j entry for given k
    kj_iterator j_begin(const k_reverse_iterator& k) const;

    //! Iterator pointing to last j entry for given k
    kj_iterator j_end(const k_reverse_iterator& k) const;

    //! Iterator pointing to first i entry for given j and k
    kji_iterator i_begin(const kj_iterator& j) const;

    //! Iterator pointing to last i entry for given j and k
    kji_iterator i_end(const kj_iterator& j) const;

    //@}

    /** \name i-based data-structure accessor methods */
    //@{

    //! Number of i entries in C(i,j,k)
    ordinal_type num_i() const;

    //! Number of k entries in C(i,j,k) for given i
    ordinal_type num_k(const i_iterator& i) const;

    //! Number of j entries in C(i,j,k) for given i and k
    ordinal_type num_j(const ik_iterator& k) const;

    //! Return i iterator for given index i
    i_iterator find_i(ordinal_type i) const;

    //! Return k iterator given i iterator and index k
    ik_iterator find_k(const i_iterator& i, ordinal_type k) const;

    //! Return j iterator given k iterator and index j
    ikj_iterator find_j(const ik_iterator& k, ordinal_type j) const;

    //! Iterator pointing to first k entry
    i_iterator i_begin() const;

    //! Iterator pointing to last k entry
    i_iterator i_end() const;

    //! Reverse iterator pointing to last k entry
    i_reverse_iterator i_rbegin() const;

    //! Reverse iterator pointing to first k entry
    i_reverse_iterator i_rend() const;

    //! Iterator pointing to first k entry for given i
    ik_iterator k_begin(const i_iterator& i) const;

    //! Iterator pointing to last k entry for given i
    ik_iterator k_end(const i_iterator& i) const;

    //! Iterator pointing to first k entry for given i
    ik_iterator k_begin(const i_reverse_iterator& i) const;

    //! Iterator pointing to last k entry for given i
    ik_iterator k_end(const i_reverse_iterator& i) const;

    //! Iterator pointing to first j entry for given i and k
    ikj_iterator j_begin(const ik_iterator& k) const;

    //! Iterator pointing to last j entry for given i and k
    ikj_iterator j_end(const ik_iterator& k) const;

    //@}

  private:

    // Prohibit copying
    Sparse3Tensor(const Sparse3Tensor&);

    // Prohibit Assignment
    Sparse3Tensor& operator=(const Sparse3Tensor& b);

  protected:

    //! Indicate whether fillComplete() has been called
    bool fill_completed;

    /** \name k-based structure */
    //@{

    //! kji indices and values in Cijk (data structure for filling)
    kji_map kji_data;

    //! kji indices and values in Cijk (data structure for iterating)
    kji_sparse_array kji_array;

    //@}

    /** \name i-based structure */
    //@{

    //! ikj indices and values in Cijk (data structure for filling)
    ikj_map ikj_data;

    //! kji indices and values in Cijk (data structure for iterating)
    ikj_sparse_array ikj_array;

    //@}

  }; // class Sparse3Tensor

  /*! \relates Sparse3Tensor
   * Print triple product tensor to output stream
   */
  template <typename ordinal_type, typename value_type>
  std::ostream& 
  operator << (std::ostream& os, 
	       const Sparse3Tensor<ordinal_type, value_type>& Cijk) {
    Cijk.print(os);
    return os;
  }

  /*! \relates Sparse3Tensor
   * Return index of a Sparse3Tensor iterator (e.g., i for a given kji_iterator)
   */
  template <typename index_iterator, typename value_iterator>
  typename SparseArrayIterator<index_iterator, value_iterator>::value_type
  index(const SparseArrayIterator<index_iterator, value_iterator>& it) {
    return *it;
  }

  /*! \relates Sparse3Tensor
   * Return index of a Sparse3Tensor reverse iterator
   */
  template <typename index_iterator, typename value_iterator>
  typename SparseArrayReverseIterator<index_iterator, value_iterator>::value_type
  index(const SparseArrayReverseIterator<index_iterator, value_iterator>& it)
  {
    return *it;
  }

  /*! \relates Sparse3Tensor
   * Return value of a Sparse3Tensor iterator (e.g., c = C(i,j,k) for a 
   * given kji_iterator)
   */
  template <typename index_iterator, typename value_iterator>
  typename SparseArrayIterator<index_iterator, value_iterator>::value_reference
  value(const SparseArrayIterator<index_iterator, value_iterator>& it) {
    return it.value();
  }

} // namespace Stokhos

// Include template definitions
#include "Stokhos_Sparse3TensorImp.hpp"

#endif // STOKHOS_SPARSE3TENSOR_HPP
