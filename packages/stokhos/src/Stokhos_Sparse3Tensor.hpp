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

#ifndef STOKHOS_SPARSE3TENSOR_HPP
#define STOKHOS_SPARSE3TENSOR_HPP

#include <ostream>
#include <map>

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

  public:

    //! Iterator for looping over k entries
    typedef typename kji_map::const_iterator k_iterator;

    //! Iterator for looping over k entries in reverse
    typedef typename kji_map::const_reverse_iterator k_reverse_iterator;

    //! Iterator for looping over j entries given k
    typedef typename ji_map::const_iterator kj_iterator;

    //! Iterator for looping over i entries given k and j
    typedef typename j_map::const_iterator kji_iterator;

    //! Iterator for looping over i entries
    typedef typename ikj_map::const_iterator i_iterator;

    //! Iterator for looping over i entries in reverse
    typedef typename ikj_map::const_reverse_iterator i_reverse_iterator;

    //! Iterator for looping over k entries given i
    typedef typename kj_map::const_iterator ik_iterator;

    //! Iterator for looping over j entries given i and k
    typedef typename j_map::const_iterator ikj_iterator;
    
    //! Constructor
    Sparse3Tensor() {}
    
    //! Destructor
    ~Sparse3Tensor() {}

    //! Add new term for given (i,j,k)
    /*!
     * Currently it fills both data structures for storing coefficients, so 
     * all accessor methods should work.
     */
    void add_term(ordinal_type i, ordinal_type j, ordinal_type k,
		  const value_type& c);

    //! Add new term for given (i,j,k) and sum in if already there
    /*!
     * Currently it does not fill the old data structures, so 
     * those accessor methods will not work.
     */
    void sum_term(ordinal_type i, ordinal_type j, ordinal_type k,
		  const value_type& c);

    //! Print tensor
    void print(std::ostream& os) const;

    /** \name k-based data-structure accessor methods */
    //@{

    //! Number of k entries in C(i,j,k)
    ordinal_type num_k() const { return kji_data.size(); }

    //! Number of j entries in C(i,j,k) for given k
    ordinal_type num_j(const k_iterator& k) const { return k->second.size(); }

    //! Number of i entries in C(i,j,k) for given k and j
    ordinal_type num_i(const kj_iterator& j) const { return j->second.size(); }

    //! Return k iterator for given index k
    k_iterator find_k(ordinal_type k) const { return kji_data.find(k); }

    //! Iterator pointing to first k entry
    k_iterator k_begin() const { return kji_data.begin(); }

    //! Iterator pointing to last k entry
    k_iterator k_end() const { return kji_data.end(); }

    //! Reverse iterator pointing to last k entry
    k_reverse_iterator k_rbegin() const { return kji_data.rbegin(); }

    //! Reverse iterator pointing to first k entry
    k_reverse_iterator k_rend() const { return kji_data.rend(); }

    //! Iterator pointing to first j entry for given k
    kj_iterator j_begin(const k_iterator& k) const { return k->second.begin(); }

    //! Iterator pointing to last j entry for given k
    kj_iterator j_end(const k_iterator& k) const { return k->second.end(); }

    //! Iterator pointing to first j entry for given k
    kj_iterator j_begin(const k_reverse_iterator& k) const { 
      return k->second.begin(); }

    //! Iterator pointing to last j entry for given k
    kj_iterator j_end(const k_reverse_iterator& k) const { 
      return k->second.end(); }

    //! Iterator pointing to first i entry for given j and k
    kji_iterator i_begin(const kj_iterator& j) const { 
      return j->second.begin(); }

    //! Iterator pointing to last i entry for given j and k
    kji_iterator i_end(const kj_iterator& j) const { return j->second.end(); }

    //@}

    /** \name i-based data-structure accessor methods */
    //@{

    //! Number of i entries in C(i,j,k)
    ordinal_type num_i() const { return ikj_data.size(); }

    //! Number of k entries in C(i,j,k) for given i
    ordinal_type num_k(const i_iterator& i) const { return i->second.size(); }

    //! Number of j entries in C(i,j,k) for given i and k
    ordinal_type num_j(const ik_iterator& k) const { return k->second.size(); }

    //! Return i iterator for given index i
    i_iterator find_i(ordinal_type i) const { return ikj_data.find(i); }

    //! Iterator pointing to first k entry
    i_iterator i_begin() const { return ikj_data.begin(); }

    //! Iterator pointing to last k entry
    i_iterator i_end() const { return ikj_data.end(); }

    //! Reverse iterator pointing to last k entry
    i_reverse_iterator i_rbegin() const { return ikj_data.rbegin(); }

    //! Reverse iterator pointing to first k entry
    i_reverse_iterator i_rend() const { return ikj_data.rend(); }

    //! Iterator pointing to first k entry for given i
    ik_iterator k_begin(const i_iterator& i) const { return i->second.begin(); }

    //! Iterator pointing to last k entry for given i
    ik_iterator k_end(const i_iterator& i) const { return i->second.end(); }

    //! Iterator pointing to first k entry for given i
    ik_iterator k_begin(const i_reverse_iterator& i) const { 
      return i->second.begin(); }

    //! Iterator pointing to last k entry for given i
    ik_iterator k_end(const i_reverse_iterator& i) const { 
      return i->second.end(); }

    //! Iterator pointing to first j entry for given i and k
    ikj_iterator j_begin(const ik_iterator& k) const { 
      return k->second.begin(); }

    //! Iterator pointing to last j entry for given i and k
    ikj_iterator j_end(const ik_iterator& k) const { return k->second.end(); }

    //@}

  private:

    // Prohibit copying
    Sparse3Tensor(const Sparse3Tensor&);

    // Prohibit Assignment
    Sparse3Tensor& operator=(const Sparse3Tensor& b);

  protected:

    /** \name k-based structure */
    //@{

    //! kji indices and values in Cijk
    kji_map kji_data;

    //@}

    /** \name i-based structure */
    //@{

    //! ikj indices and values in Cijk
    ikj_map ikj_data;

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
  template <typename iterator_type>
  typename std::iterator_traits<iterator_type>::value_type::first_type
  index(const iterator_type& it) { return it->first; }

  /*! \relates Sparse3Tensor
   * Return value of a Sparse3Tensor iterator (e.g., c = C(i,j,k) for a 
   * given kji_iterator)
   */
  template <typename iterator_type>
  typename std::iterator_traits<iterator_type>::value_type::second_type
  value(const iterator_type& it) { return it->second; }

} // namespace Stokhos

// Include template definitions
#include "Stokhos_Sparse3TensorImp.hpp"

#endif // STOKHOS_SPARSE3TENSOR_HPP
