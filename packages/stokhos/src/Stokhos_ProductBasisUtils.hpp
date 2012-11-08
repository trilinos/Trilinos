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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_PRODUCT_BASIS_UTILS_HPP
#define STOKHOS_PRODUCT_BASIS_UTILS_HPP

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Stokhos_SDMUtils.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

namespace Stokhos {

  /*!
   * \brief Compute bionomial coefficient (n ; k) = n!/( k! (n-k)! )
   */
  template <typename ordinal_type>
  ordinal_type n_choose_k(const ordinal_type& n, const ordinal_type& k) {
    // Use formula
    // n!/(k!(n-k)!) = n(n-1)...(k+1) / ( (n-k)(n-k-1)...1 )  ( n-k terms )
    //               = n(n-1)...(n-k+1) / ( k(k-1)...1 )      ( k terms )
    // which ever has fewer terms
    if (k > n)
      return 0;
    ordinal_type num = 1;
    ordinal_type den = 1;
    ordinal_type l = std::min(n-k,k);
    for (ordinal_type i=0; i<l; i++) {
      num *= n-i;
      den *= i+1;
    }
    return num / den;
  }

  //! A multidimensional index
  template <typename ordinal_t>
  class MultiIndex {
  public:

    typedef ordinal_t ordinal_type;
    typedef ordinal_t element_type;

    //! Constructor
    MultiIndex() {}

    //! Constructor
    MultiIndex(ordinal_type dim, ordinal_type v = ordinal_type(0)) : 
      index(dim,v) {}

    //! Destructor
    ~MultiIndex() {}

    //! Dimension
    ordinal_type dimension() const { return index.size(); }

    //! Size
    ordinal_type size() const { return index.size(); }

    //! Term access
    const ordinal_type& operator[] (ordinal_type i) const { return index[i]; }

    //! Term access
    ordinal_type& operator[] (ordinal_type i) { return index[i]; }

     //! Term access
    const Teuchos::Array<element_type>& getTerm() const { return index; }

    //! Term access
    Teuchos::Array<element_type>& getTerm() { return index; }

    //! Initialize
    void init(ordinal_type v) {
      for (ordinal_type i=0; i<dimension(); i++) 
	index[i] = v;
    }

    //! Resize
    void resize(ordinal_type d, ordinal_type v = ordinal_type(0)) {
      index.resize(d,v);
    }

    //! Compute total order of index
    ordinal_type order() const {
      ordinal_type my_order = 0;
      for (ordinal_type i=0; i<dimension(); ++i) my_order += index[i];
      return my_order;
    }

    //! Compare equality
    bool operator==(const MultiIndex& idx) const {
      if (dimension() != idx.dimension())
        return false;
      for (ordinal_type i=0; i<dimension(); i++) {
        if (index[i] != idx.index[i])
	  return false;
       }
       return true;
    }

    //! Compare equality
    bool operator!=(const MultiIndex& idx) const { return !(*this == idx); }

    //! Compare term-wise less-than or equal-to
    bool termWiseLEQ(const MultiIndex& idx) const {
      for (ordinal_type i=0; i<dimension(); i++) {
        if (index[i] > idx.index[i])
	  return false;
       }
       return true;
    }

    //! Print multiindex
    std::ostream& print(std::ostream& os) const {
      os << "[ ";
      for (ordinal_type i=0; i<dimension(); i++)
	os << index[i] << " ";
      os << "]";
      return os;
    }

  protected:

    //! index terms
    Teuchos::Array<ordinal_type> index;

  };

  template <typename ordinal_type>
  std::ostream& operator << (std::ostream& os, 
			     const MultiIndex<ordinal_type>& m) {
    return m.print(os);
  }

  //! An isotropic total order index set
  /*!
   * Represents the set l <= |i| <= u given upper and lower bounds
   * l and u, and |i| = i_1 + ... + i_d where d is the dimension of the
   * index.
   *
   * Currently this class only really provides an input iterator for 
   * iterating over the elements of the index set.  One should not make
   * any assumption on the order of these elements.
   */
  template <typename ordinal_t>
  class TotalOrderIndexSet {
  public:

    // Forward declaration of our iterator
    class Iterator;

    typedef ordinal_t ordinal_type;
    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef Iterator iterator;
    typedef Iterator const_iterator;

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, \c lower_ is the lower
     * bound of the index set, and \c upper_ is the upper bound (inclusive)
     */
    TotalOrderIndexSet(ordinal_type dim_, 
		       ordinal_type lower_, 
		       ordinal_type upper_) :
      dim(dim_), lower(lower_), upper(upper_) {}

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, the lower bound is zero, 
     * and \c upper_ is the upper bound (inclusive)
     */
    TotalOrderIndexSet(ordinal_type dim_, 
		       ordinal_type upper_) :
      dim(dim_), lower(0), upper(upper_) {}

    //! Return dimension
    ordinal_type dimension() const { return dim; }

    //! Return maximum order for each dimension
    multiindex_type max_orders() const {
      return multiindex_type(dim, upper);
    }

    //! Return iterator for first element in the set
    const_iterator begin() const { 
      multiindex_type index(dim);
      index[0] = lower;
      return Iterator(upper, index); 
    }

    //! Return iterator for end of the index set
    const_iterator end() const { 
      multiindex_type index(dim);
      index[dim-1] = upper+1;
      return Iterator(upper, index); 
    }

  protected:

    //! Dimension
    ordinal_type dim;

    //! Lower order of index set
    ordinal_type lower;

    //! Upper order of index set
    ordinal_type upper;

  public:

    //! Iterator class for iterating over elements of the index set
    class Iterator : public std::iterator<std::input_iterator_tag,
					  multiindex_type> {
    public:

      typedef std::iterator<std::input_iterator_tag,multiindex_type> base_type;
      typedef typename base_type::iterator_category iterator_category;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::difference_type difference_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer pointer;

      typedef const multiindex_type& const_reference;
      typedef const multiindex_type* const_pointer;

      //! Constructor
      /*!
       * \c max_order_ is the maximum order of the set (inclusive) and
       * \c index_ is the starting multi-index.
       */
      Iterator(ordinal_type max_order_, const multiindex_type& index_) : 
	max_order(max_order_), index(index_), dim(index.dimension()), 
	orders(dim) {
	orders[dim-1] = max_order;
	for (ordinal_type i=dim-2; i>=0; --i)
	  orders[i] = orders[i+1] - index[i+1];
      }

      //! Compare equality of iterators
      bool operator==(const Iterator& it) const { return index == it.index; }

      //! Compare inequality of iterators
      bool operator!=(const Iterator& it) const { return index != it.index; }
      
      //! Dereference
      const_reference operator*() const { return index; }

      //! Dereference
      const_pointer operator->() const { return &index; }
      
      //! Prefix increment, i.e., ++iterator
      /*!
       * No particular ordering of the indices is guaranteed.  The current 
       * implementation produces multi-indices sorted lexographically 
       * backwards among the elements, e.g., 
       * [0 0], [1 0], [2 0], ... [0 1], [1 1], [2 1], ... 
       * but one shouldn't assume that.  To obtain a specific
       * ordering, one should implement a "less" functional and put the 
       * indices in a sorted container such as std::map<>.
       */
      Iterator& operator++() {
	++index[0];
	ordinal_type i=0;
	while (i<dim-1 && index[i] > orders[i]) {
	  index[i] = 0;
	  ++i;
	  ++index[i];
	}
	for (ordinal_type i=dim-2; i>=0; --i)
	  orders[i] = orders[i+1] - index[i+1];

	return *this;
      }

      //! Postfix increment, i.e., iterator++
      Iterator& operator++(int) {
	Iterator tmp(*this);
	++(*this);
	return tmp;
      }

    protected:

      //! Maximum order of iterator
      ordinal_type max_order;
      
      //! Current value of iterator
      multiindex_type index;

      //! Dimension
      ordinal_type dim;

      //! Maximum orders for each term to determine how to increment
      Teuchos::Array<ordinal_type> orders;
    };
  };

  //! An anisotropic total order index set
  /*!
   * Represents the set l <= |i| <= u and i_j <= u_j given upper and 
   * lower order bounds l and u, upper component bounds u_j, 
   * and |i| = i_1 + ... + i_d where d is the dimension of the index.
   *
   * Currently this class only really provides an input iterator for 
   * iterating over the elements of the index set.  One should not make
   * any assumption on the order of these elements.
   */
  template <typename ordinal_t>
  class AnisotropicTotalOrderIndexSet {
  public:

    // Forward declaration of our iterator
    class Iterator;

    typedef ordinal_t ordinal_type;
    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef Iterator iterator;
    typedef Iterator const_iterator;

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, \c lower_ is the lower
     * bound of the index set, and \c upper_ is the upper bound (inclusive)
     */
    AnisotropicTotalOrderIndexSet(ordinal_type upper_order_,
				  const multiindex_type& lower_,
				  const multiindex_type& upper_) :
      dim(lower_.dimension()), 
      upper_order(upper_order_),
      lower(lower_), 
      upper(upper_) {}

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, the lower bound is zero, 
     * and \c upper_ is the upper bound (inclusive)
     */
    AnisotropicTotalOrderIndexSet(ordinal_type upper_order_,
				  const multiindex_type& upper_) :
      dim(upper_.dimension()), 
      upper_order(upper_order_),
      lower(dim,0),
      upper(upper_) {}

    //! Return dimension
    ordinal_type dimension() const { return dim; }

    //! Return maximum order for each dimension
    multiindex_type max_orders() const { return upper; }

    //! Return iterator for first element in the set
    const_iterator begin() const { 
      return Iterator(upper_order, upper, lower); 
    }

    //! Return iterator for end of the index set
    const_iterator end() const { 
      multiindex_type index(dim);
      index[dim-1] = std::min(upper_order, upper[dim-1]) + 1;
      return Iterator(upper_order, upper, index); 
    }

  protected:

    //! Dimension
    ordinal_type dim;

    //! Lower order of index set
    ordinal_type lower_order;

    //! Upper order of index set
    ordinal_type upper_order;

    //! Component-wise lower bounds
    multiindex_type lower;

    //! Component-wise upper bounds
    multiindex_type upper;

  public:

    //! Iterator class for iterating over elements of the index set
    class Iterator : public std::iterator<std::input_iterator_tag,
					  multiindex_type> {
    public:

      typedef std::iterator<std::input_iterator_tag,multiindex_type> base_type;
      typedef typename base_type::iterator_category iterator_category;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::difference_type difference_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer pointer;

      typedef const multiindex_type& const_reference;
      typedef const multiindex_type* const_pointer;

      //! Constructor
      /*!
       * \c max_order_ is the maximum order of the set (inclusive) and
       * \c index_ is the starting multi-index.
       */
      Iterator(ordinal_type max_order_, 
	       const multiindex_type& component_max_order_,
	       const multiindex_type& index_) : 
	max_order(max_order_), 
	component_max_order(component_max_order_), 
	index(index_), 
	dim(index.dimension()), 
	orders(dim) 
      {
	orders[dim-1] = max_order;
	for (ordinal_type i=dim-2; i>=0; --i)
	  orders[i] = orders[i+1] - index[i+1];
      }

      //! Compare equality of iterators
      bool operator==(const Iterator& it) const { return index == it.index; }

      //! Compare inequality of iterators
      bool operator!=(const Iterator& it) const { return index != it.index; }
      
      //! Dereference
      const_reference operator*() const { return index; }

      //! Dereference
      const_pointer operator->() const { return &index; }
      
      //! Prefix increment, i.e., ++iterator
      /*!
       * No particular ordering of the indices is guaranteed.  The current 
       * implementation produces multi-indices sorted lexographically 
       * backwards among the elements, e.g., 
       * [0 0], [1 0], [2 0], ... [0 1], [1 1], [2 1], ... 
       * but one shouldn't assume that.  To obtain a specific
       * ordering, one should implement a "less" functional and put the 
       * indices in a sorted container such as std::map<>.
       */
      Iterator& operator++() {
	++index[0];
	ordinal_type i=0;
	while (i<dim-1 && (index[i] > orders[i] || index[i] > component_max_order[i])) {
	  index[i] = 0;
	  ++i;
	  ++index[i];
	}
	for (ordinal_type i=dim-2; i>=0; --i)
	  orders[i] = orders[i+1] - index[i+1];

	return *this;
      }

      //! Postfix increment, i.e., iterator++
      Iterator& operator++(int) {
	Iterator tmp(*this);
	++(*this);
	return tmp;
      }

    protected:

      //! Maximum order of iterator
      ordinal_type max_order;

      //! Maximum order for each component
      multiindex_type component_max_order;
      
      //! Current value of iterator
      multiindex_type index;

      //! Dimension
      ordinal_type dim;

      //! Maximum orders for each term to determine how to increment
      Teuchos::Array<ordinal_type> orders;
    };
  };

  //! A tensor product index set
  /*!
   * Represents the set l_j <= i_j <= u_j given upper and lower bounds
   * l_j and u_j, for j = 1,...,d where d is the dimension of the
   * index.
   *
   * Currently this class only really provides an input iterator for 
   * iterating over the elements of the index set.  One should not make
   * any assumption on the order of these elements.
   */
  template <typename ordinal_t>
  class TensorProductIndexSet {
  public:

    // Forward declaration of our iterator
    class Iterator;

    typedef ordinal_t ordinal_type;
    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef Iterator iterator;
    typedef Iterator const_iterator;

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, \c lower_ is the lower
     * bound of the index set, and \c upper_ is the upper bound (inclusive)
     */
    TensorProductIndexSet(ordinal_type dim_, 
			  ordinal_type lower_, 
			  ordinal_type upper_) :
      dim(dim_), lower(dim_,lower_), upper(dim_,upper_) {}

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, the lower bound is zero, 
     * and \c upper_ is the upper bound (inclusive)
     */
    TensorProductIndexSet(ordinal_type dim_, 
			  ordinal_type upper_) :
      dim(dim_), lower(dim_,ordinal_type(0)), upper(dim_,upper_) {}

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, \c lower_ is the lower
     * bound of the index set, and \c upper_ is the upper bound (inclusive)
     */
    TensorProductIndexSet(const multiindex_type& lower_, 
			  const multiindex_type& upper_) :
      dim(lower_.dimension()), lower(lower_), upper(upper_) {}

    //! Constructor
    /*!
     * \c dim_ is the dimension of the index set, the lower bound is zero, 
     * and \c upper_ is the upper bound (inclusive)
     */
    TensorProductIndexSet(const multiindex_type& upper_) :
      dim(upper_.dimension()), lower(dim,ordinal_type(0)), upper(upper_) {}

    //! Return dimension
    ordinal_type dimension() const { return dim; }

    //! Return maximum order for each dimension
    multiindex_type max_orders() const {
      return upper;
    }

    //! Return iterator for first element in the set
    const_iterator begin() const { 
      return Iterator(upper, lower); 
    }

    //! Return iterator for end of the index set
    const_iterator end() const { 
      multiindex_type index(dim);
      index[dim-1] = upper[dim-1]+1;
      return Iterator(upper, index); 
    }

  protected:

    //! Dimension
    ordinal_type dim;

    //! Lower bound of index set
    multiindex_type lower;

    //! Upper bound of index set
    multiindex_type upper;

  public:

    //! Iterator class for iterating over elements of the index set
    class Iterator : public std::iterator<std::input_iterator_tag,
					  multiindex_type> {
    public:

      typedef std::iterator<std::input_iterator_tag,multiindex_type> base_type;
      typedef typename base_type::iterator_category iterator_category;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::difference_type difference_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer pointer;

      typedef const multiindex_type& const_reference;
      typedef const multiindex_type* const_pointer;

      //! Constructor
      /*!
       * \c upper_ is the upper bound of the set (inclusive) and
       * \c index_ is the starting multi-index.
       */
      Iterator(const multiindex_type& upper_, const multiindex_type& index_) : 
	upper(upper_), index(index_), dim(index.dimension()) {}

      //! Compare equality of iterators
      bool operator==(const Iterator& it) const { return index == it.index; }

      //! Compare inequality of iterators
      bool operator!=(const Iterator& it) const { return index != it.index; }
      
      //! Dereference
      const_reference operator*() const { return index; }

      //! Dereference
      const_pointer operator->() const { return &index; }
      
      //! Prefix increment, i.e., ++iterator
      /*!
       * No particular ordering of the indices is guaranteed.  The current 
       * implementation produces multi-indices sorted lexographically 
       * backwards among the elements, e.g., 
       * [0 0], [1 0], [2 0], ... [0 1], [1 1], [2 1], ... 
       * but one shouldn't assume that.  To obtain a specific
       * ordering, one should implement a "less" functional and put the 
       * indices in a sorted container such as std::map<>.
       */
      Iterator& operator++() {
	++index[0];
	ordinal_type i=0;
	while (i<dim-1 && index[i] > upper[i]) {
	  index[i] = 0;
	  ++i;
	  ++index[i];
	}
	return *this;
      }

      //! Postfix increment, i.e., iterator++
      Iterator& operator++(int) {
	Iterator tmp(*this);
	++(*this);
	return tmp;
      }

    protected:

      //! Upper bound of iterator
      multiindex_type upper;
      
      //! Current value of iterator
      multiindex_type index;

      //! Dimension
      ordinal_type dim;

    };
  };

  //! Container storing a term in a generalized tensor product
  template <typename ordinal_t, typename element_t>
  class TensorProductElement {
  public:

    typedef ordinal_t ordinal_type;
    typedef element_t element_type;

    //! Default constructor
    TensorProductElement() {}

    //! Constructor
    TensorProductElement(ordinal_type dim,
			 const element_type& val = element_type(0)) : 
      term(dim,val) {}

    //! Destructor
    ~TensorProductElement() {};

    //! Return dimension
    ordinal_type dimension() const { return term.size(); }

     //! Return dimension
    ordinal_type size() const { return term.size(); }

    //! Term access
    const element_type& operator[] (ordinal_type i) const { return term[i]; }

    //! Term access
    element_type& operator[] (ordinal_type i) { return term[i]; }

    //! Term access
    const Teuchos::Array<element_type>& getTerm() const { return term; }

    //! Term access
    Teuchos::Array<element_type>& getTerm() { return term; }

    //! Convert to ArrayView
    operator Teuchos::ArrayView<element_type>() { return term; }

    //! Convert to ArrayView
    operator Teuchos::ArrayView<const element_type>() const { return term; }

    //! Compute total order of tensor product element
    element_type order() const {
      element_type my_order = 0;
      for (ordinal_type i=0; i<dimension(); ++i) my_order += term[i];
      return my_order;
    }

    //! Print multiindex
    std::ostream& print(std::ostream& os) const {
      os << "[ ";
      for (ordinal_type i=0; i<dimension(); i++)
	os << term[i] << " ";
      os << "]";
      return os;
    }

  protected:

    //! Array storing term elements
    Teuchos::Array<element_type> term;
    
  };

  template <typename ordinal_type, typename element_type>
  std::ostream& operator << (
    std::ostream& os, 
    const TensorProductElement<ordinal_type,element_type>& m) {
    return m.print(os);
  }

  /*! 
   * \brief A comparison functor implementing a strict weak ordering based
   * lexographic ordering
   */
  /*
   * Objects of type \c term_type must implement \c dimension() and 
   * \c operator[] methods, as well as contain ordinal_type and element_type 
   * nested types.
   */
  template <typename term_type, 
	    typename compare_type = std::less<typename term_type::element_type> >
  class LexographicLess {
  public:
    
    typedef term_type product_element_type;
    typedef typename term_type::ordinal_type ordinal_type;
    typedef typename term_type::element_type element_type;

    //! Constructor
    LexographicLess(const compare_type& cmp_ = compare_type()) : cmp(cmp_) {}

    //! Determine if \c a is less than \c b
    bool operator()(const term_type& a, const term_type& b) const {
      ordinal_type i=0;
      while(i < a.dimension() && i < b.dimension() && equal(a[i],b[i])) i++;

      // if a is shorter than b and the first a.dimension() elements agree
      // then a is always less than b
      if (i == a.dimension()) return i != b.dimension();

      // if a is longer than b and the first b.dimension() elements agree
      // then b is always less than a
      if (i == b.dimension()) return false;

      // a and b different at element i, a is less than b if a[i] < b[i]
      return cmp(a[i],b[i]);
    }

  protected:

    //! Element comparison functor
    compare_type cmp;

    //! Determine if two elements \c a and \c b are equal
    bool equal(const element_type& a, const element_type& b) const {
      return !(cmp(a,b) || cmp(b,a));
    }

  };

  /*! 
   * \brief A comparison functor implementing a strict weak ordering based
   * total-order ordering, recursive on the dimension.
   */
  /*
   * Objects of type \c term_type must implement \c dimension(), \c order, and 
   * \c operator[] methods, as well as contain ordinal_type and element_type 
   * nested types.
   */
  template <typename term_type, 
	    typename compare_type = std::less<typename term_type::element_type> >
  class TotalOrderLess {
  public:
    
    typedef term_type product_element_type;
    typedef typename term_type::ordinal_type ordinal_type;
    typedef typename term_type::element_type element_type;

    //! Constructor
    TotalOrderLess(const compare_type& cmp_ = compare_type()) : cmp(cmp_) {}

    bool operator()(const term_type& a, const term_type& b) const {
      element_type a_order = a.order();
      element_type b_order = b.order();
      ordinal_type i=0;
      while (i < a.dimension() && i < b.dimension() && equal(a_order,b_order)) {
	a_order -= a[i];
	b_order -= b[i];
	++i;
      }
      return cmp(a_order,b_order);
    }

  protected:

    //! Element comparison functor
    compare_type cmp;

    //! Determine if two elements \c a and \c b are equal
    bool equal(const element_type& a, const element_type& b) const {
      return !(cmp(a,b) || cmp(b,a));
    }

  };

  //! A functor for comparing floating-point numbers to some tolerance
  /*!
   * The difference between this and standard < is that if |a-b| < tol, 
   * then this always returns false as a and b are treated as equal.
   */
  template <typename value_type>
  class FloatingPointLess {
  public:
    
    //! Constructor
    FloatingPointLess(const value_type& tol_ = 1.0e-12) : tol(tol_) {}

    //! Destructor
    ~FloatingPointLess() {}

    //! Compare if a < b
    bool operator() (const value_type& a, const value_type& b) const {
      return a < b - tol;
    }

  protected:

    //! Tolerance
    value_type tol;

  };

  //! A growth rule that is the identity
  template <typename value_type>
  class IdentityGrowthRule {
  public:
    //! Constructor
    IdentityGrowthRule() {}

    //! Destructor
    ~IdentityGrowthRule() {}

    //! Evaluate growth rule
    value_type operator() (const value_type& x) const { return x; }
  };

  //! Predicate functor for building sparse triple products
  template <typename ordinal_type>
  struct TensorProductPredicate {
    typedef MultiIndex<ordinal_type> coeff_type;
    coeff_type orders;

    TensorProductPredicate(const coeff_type& orders_) : orders(orders_) {}

    bool operator() (const coeff_type& term) const { 
      return term.termWiseLEQ(orders); 
    }

  };

  //! Predicate functor for building sparse triple products based on total order
  template <typename ordinal_type>
  struct TotalOrderPredicate {
    typedef MultiIndex<ordinal_type> coeff_type;
    ordinal_type max_order;
    coeff_type orders;

    TotalOrderPredicate(ordinal_type max_order_, const coeff_type& orders_) 
      : max_order(max_order_), orders(orders_) {}

    bool operator() (const coeff_type& term) const { 
      return term.termWiseLEQ(orders) && term.order() <= max_order; 
    }

  };

  /*! 
   * \brief Utilities for indexing a multi-variate complete polynomial basis 
   */
  /*!
   * This version allows specification of a growth rule for each dimension
   * allowing the coefficient order to be a function of the corresponding
   * index.
   */
  class ProductBasisUtils {
  public:

    /*!
     * \brief Generate a product basis from an index set.
     */
    template <typename index_set_type, 
	      typename growth_rule_type,
	      typename basis_set_type, 
	      typename basis_map_type>
    static void
    buildProductBasis(const index_set_type& index_set,
		      const growth_rule_type& growth_rule,
		      basis_set_type& basis_set,
		      basis_map_type& basis_map) {

      typedef typename index_set_type::ordinal_type ordinal_type;
      typedef typename index_set_type::iterator index_iterator_type;
      typedef typename basis_set_type::iterator basis_iterator_type;
      typedef typename basis_set_type::key_type coeff_type;
      
      ordinal_type dim = index_set.dimension();

      // Iterator over elements of index set
      index_iterator_type index_iterator = index_set.begin();
      index_iterator_type index_iterator_end = index_set.end();
      for (; index_iterator != index_iterator_end; ++index_iterator) {

	// Generate product coefficient
	coeff_type coeff(dim);
	for (ordinal_type i=0; i<dim; i++)
	  coeff[i] = growth_rule[i]((*index_iterator)[i]);

	// Insert coefficient into set
	basis_set[coeff] = ordinal_type(0);
      }

      // Generate linear ordering of basis_set elements
      basis_map.resize(basis_set.size());
      ordinal_type idx = 0;
      basis_iterator_type basis_iterator = basis_set.begin();
      basis_iterator_type basis_iterator_end = basis_set.end();
      for (; basis_iterator != basis_iterator_end; ++basis_iterator) {
	basis_iterator->second = idx;
	basis_map[idx] = basis_iterator->first;
	++idx;
      }
    }

    /*!
     * \brief Generate a product basis from an index set.
     */
    
    template <typename index_set_type, 
	      typename basis_set_type, 
	      typename basis_map_type>
    static void
    buildProductBasis(const index_set_type& index_set,
		      basis_set_type& basis_set,
		      basis_map_type& basis_map) {
      typedef typename index_set_type::ordinal_type ordinal_type;
      ordinal_type dim = index_set.dimension();
      Teuchos::Array< IdentityGrowthRule<ordinal_type> > growth_rule(dim);
      buildProductBasis(index_set, growth_rule, basis_set, basis_map);
    }

    template <typename ordinal_type, 
	      typename value_type,
	      typename basis_set_type, 
	      typename basis_map_type,
	      typename coeff_predicate_type,
	      typename k_coeff_predicate_type>
    static Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> >
    computeTripleProductTensor(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > >& bases,
      const basis_set_type& basis_set,
      const basis_map_type& basis_map,
      const coeff_predicate_type& coeff_pred,
      const k_coeff_predicate_type& k_coeff_pred,
      const value_type sparse_tol = 1.0e-12)
      {
#ifdef STOKHOS_TEUCHOS_TIME_MONITOR
	TEUCHOS_FUNC_TIME_MONITOR("Stokhos: Total Triple-Product Tensor Time");
#endif
	typedef typename basis_map_type::value_type coeff_type;
	ordinal_type d = bases.size();

	// The algorithm for computing Cijk = < \Psi_i \Psi_j \Psi_k > here 
	// works by factoring 
	// < \Psi_i \Psi_j \Psi_k > = 
	//    < \psi^1_{i_1}\psi^1_{j_1}\psi^1_{k_1} >_1 x ... x
	//    < \psi^d_{i_d}\psi^d_{j_d}\psi^d_{k_d} >_d
	// We compute the sparse triple product < \psi^l_i\psi^l_j\psi^l_k >_l
	// for each dimension l, and then compute all non-zero products of these
	// terms.  The complexity arises from iterating through all possible
	// combinations, throwing out terms that aren't in the basis and are 
	// beyond the k-order limit provided
	Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk = 
	  Teuchos::rcp(new Sparse3Tensor<ordinal_type, value_type>);
	
	// Create 1-D triple products
	Teuchos::Array< Teuchos::RCP<Sparse3Tensor<ordinal_type,value_type> > > Cijk_1d(d);
	for (ordinal_type i=0; i<d; i++) {
	  Cijk_1d[i] = 
	    bases[i]->computeSparseTripleProductTensor(bases[i]->order()+1);
	}
	
	// Create i, j, k iterators for each dimension
	// Note:  we have to supply an initializer in the arrays of iterators 
	// to avoid checked-stl errors about singular iterators
	typedef Sparse3Tensor<ordinal_type,value_type> Cijk_type;
	typedef typename Cijk_type::k_iterator k_iterator;
	typedef typename Cijk_type::kj_iterator kj_iterator;
	typedef typename Cijk_type::kji_iterator kji_iterator;
	Teuchos::Array<k_iterator> k_iterators(d, Cijk_1d[0]->k_begin());
	Teuchos::Array<kj_iterator > j_iterators(d, Cijk_1d[0]->j_begin(k_iterators[0]));
	Teuchos::Array<kji_iterator > i_iterators(d, Cijk_1d[0]->i_begin(j_iterators[0]));
	coeff_type terms_i(d), terms_j(d), terms_k(d);
	for (ordinal_type dim=0; dim<d; dim++) {
	  k_iterators[dim] = Cijk_1d[dim]->k_begin();
	  j_iterators[dim] = Cijk_1d[dim]->j_begin(k_iterators[dim]);
	  i_iterators[dim] = Cijk_1d[dim]->i_begin(j_iterators[dim]);
	  terms_i[dim] = index(i_iterators[dim]);
	  terms_j[dim] = index(j_iterators[dim]);
	  terms_k[dim] = index(k_iterators[dim]);
	}

	ordinal_type I = 0;
	ordinal_type J = 0;
	ordinal_type K = 0;
	bool valid_i = coeff_pred(terms_i);
	bool valid_j = coeff_pred(terms_j);
	bool valid_k = k_coeff_pred(terms_k);
	bool inc_i = true;
	bool inc_j = true;
	bool inc_k = true;
	bool stop = false;
	ordinal_type cnt = 0;
	while (!stop) {
	  
	  // Add term if it is in the basis
	  if (valid_i && valid_j && valid_k) {
	    if (inc_k) {
	      typename basis_set_type::const_iterator k = 
		basis_set.find(terms_k);
	      K = k->second;
            }
	    if (inc_i) {
	      typename basis_set_type::const_iterator i = 
		basis_set.find(terms_i);
	      I = i->second;
	    }
	    if (inc_j) {
	      typename basis_set_type::const_iterator j = 
		basis_set.find(terms_j);
	      J = j->second;
	    }
	    value_type c = value_type(1.0);
	    value_type nrm = value_type(1.0);
	    for (ordinal_type dim=0; dim<d; dim++) {
	      c *= value(i_iterators[dim]);
	      nrm *= bases[dim]->norm_squared(terms_i[dim]);
	    }
	    if (std::abs(c/nrm) > sparse_tol)
	      Cijk->add_term(I,J,K,c);
	  }
	  
	  // Increment iterators to the next valid term
	  ordinal_type cdim = 0;
	  bool inc = true;
	  inc_i = false;
	  inc_j = false; 
	  inc_k = false;
	  while (inc && cdim < d) {
	    ordinal_type cur_dim = cdim;
	    ++i_iterators[cdim];
	    inc_i = true;
	    if (i_iterators[cdim] != Cijk_1d[cdim]->i_end(j_iterators[cdim])) {
	      terms_i[cdim] = index(i_iterators[cdim]);
	      valid_i = coeff_pred(terms_i);
	    }
	    if (i_iterators[cdim] == Cijk_1d[cdim]->i_end(j_iterators[cdim]) ||
		!valid_i) {
	      ++j_iterators[cdim];
	      inc_j = true;
	      if (j_iterators[cdim] != Cijk_1d[cdim]->j_end(k_iterators[cdim])) {
		terms_j[cdim] = index(j_iterators[cdim]);
		valid_j = coeff_pred(terms_j);
	      }
	      if (j_iterators[cdim] == Cijk_1d[cdim]->j_end(k_iterators[cdim]) ||
		  !valid_j) {
		++k_iterators[cdim];
		inc_k = true;
		if (k_iterators[cdim] != Cijk_1d[cdim]->k_end()) {
		  terms_k[cdim] = index(k_iterators[cdim]);
		  valid_k = k_coeff_pred(terms_k);
		}
		if (k_iterators[cdim] == Cijk_1d[cdim]->k_end() || !valid_k) {
		  k_iterators[cdim] = Cijk_1d[cdim]->k_begin();
		  ++cdim;
		  terms_k[cur_dim] = index(k_iterators[cur_dim]);
		  valid_k = k_coeff_pred(terms_k);
		}
		else
		  inc = false;
		j_iterators[cur_dim] = 
		  Cijk_1d[cur_dim]->j_begin(k_iterators[cur_dim]);
		terms_j[cur_dim] = index(j_iterators[cur_dim]);
		valid_j = coeff_pred(terms_j);
	      }
	      else
		inc = false;
	      i_iterators[cur_dim] = 
		Cijk_1d[cur_dim]->i_begin(j_iterators[cur_dim]);
	      terms_i[cur_dim] = index(i_iterators[cur_dim]);
	      valid_i = coeff_pred(terms_i);
	    }
	    else
	      inc = false;
	    
	    if (!valid_i || !valid_j || !valid_k)
	      inc = true;
	  }
	  
	  if (cdim == d)
	    stop = true;
	  
	  cnt++;
	}
	
	Cijk->fillComplete();

	return Cijk;
      }
  };

#if 0
  /*!
   * \brief An operator for building pseudo-spectral coefficients in a tensor
   * product basis using tensor-product quadrature.
   */
  template <typename coeff_compare_type,
	    typename point_compare_type>
  class TensorProductPseudoSpectralOperator {
  public:
    typedef typename coeff_compare_type::ordinal_type ordinal_type;
    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef typename coeff_compare_type::product_element_type coeff_type;
    typedef typename point_compare_type::product_element_type point_type;
    typedef typename point_type::element_type value_type;
    typedef OneDOrthogPolyBasis<ordinal_type,value_type> basis_type;
    
    typedef point_type domain;
    typedef coeff_type range;
    typedef point_compare_type domain_compare;
    typedef coeff_compare_type range_compare;
    typedef std::map<domain,std::pair<value_type,ordinal_type>,domain_compare> domain_set_type;
    typedef std::map<range,ordinal_type,range_compare> range_set_type;
    typedef Teuchos::Array<domain> domain_map_type;
    typedef Teuchos::Array<range> range_map_type;

    typedef typename domain_map_type::iterator domain_iterator;
    typedef typename domain_map_type::const_iterator domain_const_iterator;
    typedef typename range_map_type::iterator range_iterator;
    typedef typename range_map_type::const_iterator range_const_iterator;
    typedef typename domain_set_type::iterator domain_set_iterator;
    typedef typename domain_set_type::const_iterator domain_set_const_iterator;
    typedef typename range_set_type::iterator range_set_iterator;
    typedef typename range_set_type::const_iterator range_set_const_iterator;

    //! Constructor
    template <typename coeff_index_set_type>
    TensorProductPseudoSpectralOperator(
      const Teuchos::Array< Teuchos::RCP<const basis_type> >& bases_, 
      const coeff_index_set_type& coeff_index_set,
      const multiindex_type& multiindex,
      const coeff_compare_type& coeff_compare = coeff_compare_type(),
      const point_compare_type& point_compare = point_compare_type()) : 
      dim(multiindex.dimension()),
      points(point_compare),
      coeffs(coeff_compare) {

      // Generate coefficients
      Stokhos::ProductBasisUtils::buildProductBasis(
	coeff_index_set, coeffs, coeff_map);

      // Make sure order of each 1-D basis is large enough for 
      // supplied set of coefficients
      Teuchos::Array< Teuchos::RCP<const basis_type> > bases(bases_);
      multiindex_type max_orders = coeff_index_set.max_orders();
      for (ordinal_type i=0; i<dim; ++i) {
	if (bases[i]->order() < max_orders[i])
	  bases[i] = bases[i]->cloneWithOrder(max_orders[i]);
      }

      // Compute quad points, weights, values
      Teuchos::Array< Teuchos::Array<value_type> > gp(dim);
      Teuchos::Array< Teuchos::Array<value_type> > gw(dim);
      Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > gv(dim);
      multiindex_type n(dim);
      for (ordinal_type i=0; i<dim; i++) {
	bases[i]->getQuadPoints(2*multiindex[i], gp[i], gw[i], gv[i]);
	n[i] = gp[i].size()-1;
      }

      // Generate points and weights
      Stokhos::TensorProductIndexSet<ordinal_type> pointIndexSet(n);
      typedef typename Stokhos::TensorProductIndexSet<ordinal_type>::iterator index_iterator;
      index_iterator point_iterator = pointIndexSet.begin();
      index_iterator point_end = pointIndexSet.end();
      point_type point(dim);
      while (point_iterator != point_end) {
	value_type w = 1.0;
	for (ordinal_type k=0; k<dim; k++) {
	  point[k] = gp[k][(*point_iterator)[k]];
	  w *= gw[k][(*point_iterator)[k]];
	}
	points[point] = std::make_pair(w,ordinal_type(0));
	++point_iterator;
      }

      // Generate linear ordering of points
      ordinal_type nqp = points.size();
      point_map.resize(nqp);      
      ordinal_type idx=0;
      typename domain_set_type::iterator di = points.begin();
      typename domain_set_type::iterator di_end = points.end();
      while (di != di_end) {
	di->second.second = idx;
	point_map[idx] = di->first;
	++idx;
	++di;
      }

      // Generate quadrature operator
      ordinal_type npc = coeffs.size();
      A.reshape(npc,nqp);
      A.putScalar(1.0);
      for (point_iterator = pointIndexSet.begin(); point_iterator != point_end; 
	   ++point_iterator) {
	for (ordinal_type k=0; k<dim; k++)
	  point[k] = gp[k][(*point_iterator)[k]];
	ordinal_type j = points[point].second;
	for (ordinal_type i=0; i<npc; i++) {
	  for (ordinal_type k=0; k<dim; k++) {
	    ordinal_type l = (*point_iterator)[k];
	    ordinal_type m = coeff_map[i][k];
	    A(i,j) *= gw[k][l]*gv[k][l][m] / bases[k]->norm_squared(m);
	  }
	}
      }
      
    }

    ordinal_type domain_size() const { return points.size(); }
    ordinal_type range_size() const { return coeffs.size(); }

    domain_iterator domain_begin() { return point_map.begin(); }
    domain_iterator domain_end() { return point_map.end(); }
    range_iterator range_begin() { return coeff_map.begin(); }
    range_iterator range_end() { return coeff_map.end(); }

    domain_const_iterator domain_begin() const { return point_map.begin(); }
    domain_const_iterator domain_end() const { return point_map.end(); }
    range_const_iterator range_begin() const { return coeff_map.begin(); }
    range_const_iterator range_end() const { return coeff_map.end(); }

    domain_set_iterator domain_set_begin() { return points.begin(); }
    domain_set_iterator domain_set_end() { return points.end(); }
    range_set_iterator range_set_begin() { return coeffs.begin(); }
    range_set_iterator range_set_end() { return coeffs.end(); }

    domain_set_const_iterator domain_set_begin() const { return points.begin(); }
    domain_set_const_iterator domain_set_end() const { return points.end(); }
    range_set_const_iterator range_set_begin() const { return coeffs.begin(); }
    range_set_const_iterator range_set_end() const { return coeffs.end(); }

    ordinal_type getDomainIndex(const domain& term) const { 
      return points[term];
    }
    ordinal_type getRangeIndex(const range& term) const { 
      return coeffs[term];
    }

    const domain& getDomainTerm(ordinal_type n) const {
      return point_map[n];
    }
    const range& getRangeTerm(ordinal_type n) const {
      return coeff_map[n];
    }

    //! Apply tensor product quadrature operator
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     */
    void apply(const value_type& alpha, 
	       const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
	       Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
	       const value_type& beta,
	       bool trans = false) const {
      ordinal_type ret;
      if (trans)
	ret = result.multiply(Teuchos::NO_TRANS, Teuchos::TRANS, alpha, input,
			      A, beta);
      else
	ret = result.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, alpha, A,
			      input, beta);
      TEUCHOS_ASSERT(ret == 0);
    }

  protected:

    //! Dimension
    ordinal_type dim;

    //! Quadrature points
    domain_set_type points;

    //! Coefficients we compute
    range_set_type coeffs;

    //! Map index to domain term
    domain_map_type point_map;

    //! Map index to range term
    range_map_type coeff_map;

    //! Matrix mapping point to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> A;

  };

  /*! 
   * A factory for building TensorProductPseudoSpectralOperator's for a
   * given multi-index.
   */
  template <typename coeff_compare_type,
	    typename point_compare_type>
  class TensorProductPseudoSpectralOperatorFactory {
  public:
    typedef TensorProductPseudoSpectralOperator<coeff_compare_type,
						point_compare_type> operator_type;
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::basis_type basis_type;
    typedef typename operator_type::multiindex_type multiindex_type;

    //! Constructor
    TensorProductPseudoSpectralOperatorFactory(
      const Teuchos::Array< Teuchos::RCP<const basis_type > >& bases_,
      const coeff_compare_type& coeff_compare_ = coeff_compare_type(),
      const point_compare_type& point_compare_ = point_compare_type()) : 
      bases(bases_),
      coeff_compare(coeff_compare_),
      point_compare(point_compare_) {}

    //! Destructor
    ~TensorProductPseudoSpectralOperatorFactory() {}

    //! Generate operator
    template <typename coeff_index_set_type>
    Teuchos::RCP<operator_type> operator() (
      const coeff_index_set_type& coeff_index_set,
      const multiindex_type& multiindex) const {
      return Teuchos::rcp(new operator_type(bases, coeff_index_set, multiindex,
					    coeff_compare, point_compare));
    }

  protected:

    //! 1-D bases comprising tensor product
    Teuchos::Array< Teuchos::RCP<const basis_type > > bases;

    //! Coefficient comparison functor
    coeff_compare_type coeff_compare;

    //! Point comparison functor
    point_compare_type point_compare;

  };

  /*!
   * \brief An operator for building pseudo-spectral coefficients in a tensor
   * product basis using tensor-product quadrature.
   */
  /*!
   * This version uses partial summation in the quadrature.  It is unclear if
   * this can be done for arbitrary coefficient index sets, so this version
   * uses the same index set for points and coefficients.
   *
   * The partial summation relies on the fact that the matrix of multivariate
   * basis functions evaluated at multivariate quadrature points can be written
   * as a Kronecker product A = A_1 \otimes ... \A_d where d is the dimension
   * of the quadrature, and A_i is the corresponding 1-D quadrature operator.
   * Then A*x can be decomposed into a series of 1-D quadratures, which is
   * signficantly more efficient than computing A*x directly.
   *
   * Since the partial summation construction relies on a specific ordering
   * of the points and coefficients, this operator does not allow specification
   * of the point and coefficient comparison functors.
   */
  template <typename ordinal_t, typename value_t>
  class TensorProductPseudoSpectralOperatorPST {
  public:

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef OneDOrthogPolyBasis<ordinal_type,value_type> basis_type;
    typedef TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef TensorProductElement<ordinal_type,value_type> point_type;
    typedef MultiIndex<ordinal_type> multiindex_type;
    typedef Stokhos::TensorProductIndexSet<ordinal_type> coeff_index_set_type;

    typedef point_type domain;
    typedef coeff_type range;
    typedef FloatingPointLess<value_type> point_compare;
    typedef std::less<ordinal_type> coeff_compare;
    typedef LexographicLess<point_type,point_compare> domain_compare;
    typedef LexographicLess<coeff_type,coeff_compare> range_compare;
    typedef std::map<domain,std::pair<value_type,ordinal_type>,domain_compare> domain_set_type;
    typedef std::map<range,ordinal_type,range_compare> range_set_type;
    typedef Teuchos::Array<domain> domain_map_type;
    typedef Teuchos::Array<range> range_map_type;

    typedef typename domain_map_type::iterator domain_iterator;
    typedef typename domain_map_type::const_iterator domain_const_iterator;
    typedef typename range_map_type::iterator range_iterator;
    typedef typename range_map_type::const_iterator range_const_iterator;
    typedef typename domain_set_type::iterator domain_set_iterator;
    typedef typename domain_set_type::const_iterator domain_set_const_iterator;
    typedef typename range_set_type::iterator range_set_iterator;
    typedef typename range_set_type::const_iterator range_set_const_iterator;

    //! Constructor
    TensorProductPseudoSpectralOperatorPST(
      const Teuchos::Array< Teuchos::RCP<const basis_type> >& bases_,
      const coeff_index_set_type& coeff_index_set,
      const multiindex_type& multiindex,
      const value_type& point_tol=1e-12) : 
      dim(multiindex.dimension()),
      points(domain_compare(point_compare(point_tol))),
      coeffs(range_compare(coeff_compare())),
      Ak(dim) {

      // Generate coefficients
      Stokhos::ProductBasisUtils::buildProductBasis(
	coeff_index_set, coeffs, coeff_map);

      // Make sure order of each 1-D basis is large enough for 
      // supplied set of coefficients
      Teuchos::Array< Teuchos::RCP<const basis_type> > bases(bases_);
      multiindex_type max_orders = coeff_index_set.max_orders();
      for (ordinal_type i=0; i<dim; ++i) {
	if (bases[i]->order() != max_orders[i])
	  bases[i] = bases[i]->cloneWithOrder(max_orders[i]);
      }

      // Compute quad points, weights, values
      Teuchos::Array< Teuchos::Array<value_type> > gp(dim);
      Teuchos::Array< Teuchos::Array<value_type> > gw(dim);
      Teuchos::Array< Teuchos::Array< Teuchos::Array<value_type> > > gv(dim);
      multiindex_type n(dim);
      for (ordinal_type k=0; k<dim; k++) {
	bases[k]->getQuadPoints(2*multiindex[k], gp[k], gw[k], gv[k]);
	n[k] = gp[k].size()-1;

	// Generate quadrature operator
	ordinal_type npc = bases[k]->size();
	ordinal_type nqp = gp[k].size();
	Ak[k].reshape(npc,nqp);
	Ak[k].putScalar(1.0);
	for (ordinal_type j=0; j<nqp; j++) {
	  for (ordinal_type i=0; i<npc; i++) {
	      Ak[k](i,j) *= gw[k][j]*gv[k][j][i] / bases[k]->norm_squared(i);
	  }
	}
      }

      // Generate points and weights
      Stokhos::TensorProductIndexSet<ordinal_type> pointIndexSet(n);
      typedef typename Stokhos::TensorProductIndexSet<ordinal_type>::iterator index_iterator;
      index_iterator point_iterator = pointIndexSet.begin();
      index_iterator point_end = pointIndexSet.end();
      point_type point(dim);
      while (point_iterator != point_end) {
	value_type w = 1.0;
	for (ordinal_type k=0; k<dim; k++) {
	  point[k] = gp[k][(*point_iterator)[k]];
	  w *= gw[k][(*point_iterator)[k]];
	}
	points[point] = std::make_pair(w,ordinal_type(0));
	++point_iterator;
      }

      // Generate linear ordering of points
      ordinal_type nqp = points.size();
      point_map.resize(nqp);
      ordinal_type idx=0;
      typename domain_set_type::iterator di = points.begin();
      typename domain_set_type::iterator di_end = points.end();
      while (di != di_end) {
	di->second.second = idx;
	point_map[idx] = di->first;
	++idx;
	++di;
      }
      
    }

    ordinal_type domain_size() const { return points.size(); }
    ordinal_type range_size() const { return coeffs.size(); }

    domain_iterator domain_begin() { return point_map.begin(); }
    domain_iterator domain_end() { return point_map.end(); }
    range_iterator range_begin() { return coeff_map.begin(); }
    range_iterator range_end() { return coeff_map.end(); }

    domain_const_iterator domain_begin() const { return point_map.begin(); }
    domain_const_iterator domain_end() const { return point_map.end(); }
    range_const_iterator range_begin() const { return coeff_map.begin(); }
    range_const_iterator range_end() const { return coeff_map.end(); }

    domain_set_iterator domain_set_begin() { return points.begin(); }
    domain_set_iterator domain_set_end() { return points.end(); }
    range_set_iterator range_set_begin() { return coeffs.begin(); }
    range_set_iterator range_set_end() { return coeffs.end(); }

    domain_set_const_iterator domain_set_begin() const { return points.begin(); }
    domain_set_const_iterator domain_set_end() const { return points.end(); }
    range_set_const_iterator range_set_begin() const { return coeffs.begin(); }
    range_set_const_iterator range_set_end() const { return coeffs.end(); }

    ordinal_type getDomainIndex(const domain& term) const { 
      return points[term];
    }
    ordinal_type getRangeIndex(const range& term) const { 
      return coeffs[term];
    }

    const domain& getDomainTerm(ordinal_type n) const {
      return point_map[n];
    }
    const range& getRangeTerm(ordinal_type n) const {
      return coeff_map[n];
    }

    //! Apply tensor product quadrature operator
    /*!
     * \c input is a vector storing values of a function at the quadrature
     * points, and \c result will contain the resulting polynomial chaos
     * coefficients.  \c input and \c result can have multiple columns for
     * vector-valued functions and set \c trans to true if these (multi-) 
     * vectors are layed out in a transposed fashion.
     *
     * For k=1,...,d, let A_k be the m_k-by-n_k quadrature operator for
     * dimension k and A = A_1 \otimes ... \otimes A_d be m-\by-n where 
     * m = m_1...m_d and n = n_1...n_d.  For any two matrices
     * B and C and vector x, (B \otimes C)x = vec( C*X*B^T ) = 
     * vec( C*(B*X^T)^T ) where x = vec(X), and the vec() operator makes a 
     * vector out of a matrix by stacking columns.  Applying this formula
     * recursively to A yields the simple algorithm for computing y = A*x
     * (x is n-by-1 and y is m-by-1):
     *    X = x;
     *    for k=1:d
     *      n = n / n_k;
     *      X = reshape(X, n, n_k);
     *      X = A_k*X';
     *      n = n * m_k;
     *    end
     *    y = reshape(X, m, 1);
     *
     * When x has p columns, it is somehwat more complicated because the
     * standard transpose above isn't correct.  Instead the transpose
     * needs to be applied to each p block in X.
     *
     * The strategy here for dealing with transposed input and result
     * is to transpose them when copying to/from the temporary buffers
     * used in the algorithm.
     */
    void apply(
      const value_type& alpha, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input,
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result, 
      const value_type& beta,
      bool trans = false) const {

      ordinal_type n, m, p;
      if (trans) {
	TEUCHOS_ASSERT(input.numRows() == result.numRows());
	n = input.numCols();
	p = input.numRows();
	m = result.numCols();
      }
      else {
	TEUCHOS_ASSERT(input.numCols() == result.numCols());
	n = input.numRows();
	p = input.numCols();
	m = result.numRows();
      }    
      TEUCHOS_ASSERT(n == domain_size());
      TEUCHOS_ASSERT(m == range_size());
      ordinal_type sz = std::max(n,m);  

      // Temporary buffers used in algorithm
      Teuchos::Array<value_type> tmp1(sz*p),  tmp2(sz*p);

      // Copy input to tmp1 (transpose if necessary)
      if (trans) {
	for (ordinal_type j=0; j<p; j++)
	  for (ordinal_type i=0; i<n; i++)
	    tmp1[i+j*n] = input(j,i);
      }
      else {
	for (ordinal_type j=0; j<p; j++)
	  for (ordinal_type i=0; i<n; i++)
	    tmp1[i+j*n] = input(i,j);
      }
      
      // Loop over each term in Kronecker product
      for (ordinal_type k=0; k<dim; k++) {
	ordinal_type mk = Ak[k].numRows();
	ordinal_type nk = Ak[k].numCols();
	n = n / nk;
	
	// x = reshape(x, n, n_k)
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(
	  Teuchos::View, &tmp1[0], n, n, nk*p);

	// y = x'
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> y(
	  Teuchos::View, &tmp2[0], nk, nk, n*p);
	for (ordinal_type l=0; l<p; l++)
	  for (ordinal_type j=0; j<n; j++)
	    for (ordinal_type i=0; i<nk; i++)
	      y(i,j+l*n) = x(j,i+l*nk);

	// x = A_k*x
	Teuchos::SerialDenseMatrix<ordinal_type,value_type> z(
	  Teuchos::View, &tmp1[0], mk, mk, n*p);
	ordinal_type ret = 
	  z.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, Ak[k], y, 0.0);
	TEUCHOS_ASSERT(ret == 0);

	n = n * mk;
      }

      // Sum answer into result (transposing if necessary)
      if (trans) {
	for (ordinal_type j=0; j<p; j++)
	  for (ordinal_type i=0; i<m; i++)
	    result(j,i) = beta*result(j,i) + alpha*tmp1[i+j*m];
      }
      else {
	for (ordinal_type j=0; j<p; j++)
	  for (ordinal_type i=0; i<m; i++)
	    result(i,j) = beta*result(i,j) + alpha*tmp1[i+j*m];
      }
    }

  protected:

    //! Dimension
    ordinal_type dim;

    //! Quadrature points
    domain_set_type points;

    //! Coefficients we compute
    range_set_type coeffs;

    //! Map index to domain term
    domain_map_type point_map;

    //! Map index to range term
    range_map_type coeff_map;

    //! Matrix mapping point to coefficients for each dimension
    Teuchos::Array< Teuchos::SerialDenseMatrix<ordinal_type,value_type> > Ak;

    //! Matrix mapping point to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> A;

    Teuchos::BLAS<ordinal_type,value_type> blas;

  };

  /*! 
   * A factory for building TensorProductPseudoSpectralOperator's for a
   * given multi-index.
   */
  template <typename ordinal_t, typename value_type>
  class TensorProductPseudoSpectralOperatorPSTFactory {
  public:
    typedef TensorProductPseudoSpectralOperatorPST<ordinal_t,
						   value_type> operator_type;
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::basis_type basis_type;
    typedef typename operator_type::multiindex_type multiindex_type;
    typedef typename operator_type::coeff_index_set_type coeff_index_set_type;

    //! Constructor
    TensorProductPseudoSpectralOperatorPSTFactory(
      const Teuchos::Array< Teuchos::RCP<const basis_type > >& bases_,
      const value_type& point_tol_=1e-12) : 
      bases(bases_),
      point_tol(point_tol_) {}

    //! Destructor
    ~TensorProductPseudoSpectralOperatorPSTFactory() {}

    //! Generate operator
    Teuchos::RCP<operator_type> operator() (
      const coeff_index_set_type& coeff_index_set,
      const multiindex_type& multiindex) const {
      return Teuchos::rcp(new operator_type(bases, coeff_index_set, multiindex, 
					    point_tol));
    }

  protected:

    //! 1-D bases comprising tensor product
    Teuchos::Array< Teuchos::RCP<const basis_type > > bases;

    //! Point tolerance
    value_type point_tol;

  };

#endif

  /*! 
   * \brief Utilities for indexing a multi-variate complete polynomial basis 
   */
  template <typename ordinal_type, typename value_type>
  class CompletePolynomialBasisUtils {
  public:

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension.
     */
    /*
     * Returns expansion total order.  
     * This version is for an isotropic expansion of total order \c p in
     * \c d dimensions.
     */
    static ordinal_type 
    compute_terms(ordinal_type p, ordinal_type d, 
		  ordinal_type& sz,
		  Teuchos::Array< MultiIndex<ordinal_type> >& terms,
		  Teuchos::Array<ordinal_type>& num_terms);

    /*!
     * \brief Compute the 2-D array of basis terms which maps a basis index
     * into the orders for each basis dimension
     */
    /*
     * Returns expansion total order.  
     * This version allows for anisotropy in the maximum order in each 
     * dimension.
     */
    static ordinal_type 
    compute_terms(const Teuchos::Array<ordinal_type>& basis_orders, 
		  ordinal_type& sz,
		  Teuchos::Array< MultiIndex<ordinal_type> >& terms,
		  Teuchos::Array<ordinal_type>& num_terms);

    /*!
     * \brief Compute basis index given the orders for each basis
     * dimension.
     */
    static ordinal_type 
    compute_index(const MultiIndex<ordinal_type>& term,
		  const Teuchos::Array< MultiIndex<ordinal_type> >& terms,
		  const Teuchos::Array<ordinal_type>& num_terms,
		  ordinal_type max_p);

  };

}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_terms(ordinal_type p, ordinal_type d, 
	      ordinal_type& sz,
	      Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms,
	      Teuchos::Array<ordinal_type>& num_terms)
{
  Teuchos::Array<ordinal_type> basis_orders(d, p);
  return compute_terms(basis_orders, sz, terms, num_terms);
}

template <typename ordinal_type, typename value_type>
ordinal_type
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_terms(const Teuchos::Array<ordinal_type>& basis_orders, 
	      ordinal_type& sz,
	      Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms,
	      Teuchos::Array<ordinal_type>& num_terms)
{
  // The approach here for ordering the terms is inductive on the total
  // order p.  We get the terms of total order p from the terms of total
  // order p-1 by incrementing the orders of the first dimension by 1.
  // We then increment the orders of the second dimension by 1 for all of the
  // terms whose first dimension order is 0.  We then repeat for the third
  // dimension whose first and second dimension orders are 0, and so on.
  // How this is done is most easily illustrated by an example of dimension 3:
  //
  // Order  terms   cnt  Order  terms   cnt
  //   0    0 0 0          4    4 0 0  15 5 1
  //                            3 1 0
  //   1    1 0 0  3 2 1        3 0 1
  //        0 1 0               2 2 0
  //        0 0 1               2 1 1
  //                            2 0 2
  //   2    2 0 0  6 3 1        1 3 0
  //        1 1 0               1 2 1
  //        1 0 1               1 1 2
  //        0 2 0               1 0 3
  //        0 1 1               0 4 0
  //        0 0 2               0 3 1
  //                            0 2 2
  //   3    3 0 0  10 4 1       0 1 3
  //        2 1 0               0 0 4
  //        2 0 1
  //        1 2 0
  //        1 1 1
  //        1 0 2
  //        0 3 0
  //        0 2 1
  //        0 1 2
  //        0 0 3

  // Compute total order
  ordinal_type d = basis_orders.size();
  ordinal_type p = 0;
  for (ordinal_type i=0; i<d; i++) {
    if (basis_orders[i] > p)
      p = basis_orders[i];
  }

  // Temporary array of terms grouped in terms of same order
  Teuchos::Array< Teuchos::Array< MultiIndex<ordinal_type> > > terms_order(p+1);

  // Store number of terms up to each order
  num_terms.resize(p+2, ordinal_type(0));

  // Set order 0
  terms_order[0].resize(1);
  terms_order[0][0].resize(d, ordinal_type(0));
  num_terms[0] = 1;

  // The array "cnt" stores the number of terms we need to increment for each
  // dimension.  
  Teuchos::Array<ordinal_type> cnt(d), cnt_next(d);
  MultiIndex<ordinal_type> term(d);
  for (ordinal_type j=0; j<d; j++) {
    if (basis_orders[j] >= 1)
      cnt[j] = 1;
    else
      cnt[j] = 0;
    cnt_next[j] = 0;
  }

  sz = 1;
  // Loop over orders
  for (ordinal_type k=1; k<=p; k++) {

    num_terms[k] = num_terms[k-1];

    // Stores the index of the term we copying
    ordinal_type prev = 0;

    // Loop over dimensions
    for (ordinal_type j=0; j<d; j++) {

      // Increment orders of cnt[j] terms for dimension j
      for (ordinal_type i=0; i<cnt[j]; i++) {
	if (terms_order[k-1][prev+i][j] < basis_orders[j]) {
	  term = terms_order[k-1][prev+i];
	  ++term[j];
	  terms_order[k].push_back(term);
	  ++sz;
	  num_terms[k]++;
	  for (ordinal_type l=0; l<=j; l++)
	    ++cnt_next[l];
	}
      }

      // Move forward to where all orders for dimension j are 0
      if (j < d-1)
	prev += cnt[j] - cnt[j+1];

    }

    // Update the number of terms we must increment for the new order
    for (ordinal_type j=0; j<d; j++) {
      cnt[j] = cnt_next[j];
      cnt_next[j] = 0;
    }

  }

  num_terms[p+1] = sz;

  // Copy into final terms array
  terms.resize(sz);
  ordinal_type i = 0;
  for (ordinal_type k=0; k<=p; k++) {
    ordinal_type num_k = terms_order[k].size();
    for (ordinal_type j=0; j<num_k; j++)
      terms[i++] = terms_order[k][j];
  }

  return p;
}

template <typename ordinal_type, typename value_type>
ordinal_type 
Stokhos::CompletePolynomialBasisUtils<ordinal_type, value_type>::
compute_index(const Stokhos::MultiIndex<ordinal_type>& term,
	      const Teuchos::Array< Stokhos::MultiIndex<ordinal_type> >& terms,
	      const Teuchos::Array<ordinal_type>& num_terms,
	      ordinal_type max_p)
{
  // The approach here for computing the index is to find the order block
  // corresponding to this term by adding up the component orders.  We then
  // do a linear search through the terms_order array for this order

  // First compute order of term
  ordinal_type d = term.dimension();
  ordinal_type ord = 0;
  for (ordinal_type i=0; i<d; i++)
    ord += term[i];
  TEUCHOS_TEST_FOR_EXCEPTION(ord < 0 || ord > max_p, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Term has invalid order " << ord);

  // Now search through terms of that order to find a match
  ordinal_type k;
  if (ord == 0)
    k = 0;
  else
    k = num_terms[ord-1];
  ordinal_type k_max=num_terms[ord];
  bool found = false;
  while (k < k_max && !found) {
    bool found_term = true;
    for (ordinal_type j=0; j<d; j++) {
      found_term = found_term && (term[j] == terms[k][j]);
      if (!found_term)
	break;
    }
    found = found_term;
    ++k;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(k >= k_max && !found, std::logic_error,
		     "Stokhos::CompletePolynomialBasis::compute_index(): " <<
		     "Could not find specified term.");

  return k-1;
}

#endif
