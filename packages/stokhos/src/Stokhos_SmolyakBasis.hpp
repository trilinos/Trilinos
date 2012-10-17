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

#ifndef STOKHOS_SMOLYAK_BASIS_HPP
#define STOKHOS_SMOLYAK_BASIS_HPP

#include <map>

#include "Teuchos_RCP.hpp"

#include "Stokhos_ProductBasis.hpp"
#include "Stokhos_ProductBasisUtils.hpp"

namespace Stokhos {

  //! A linear growth rule
  template <typename ordinal_type>
  class LinearGrowthRule {
  public:
    //! Constructor
    LinearGrowthRule(const ordinal_type& a_ = ordinal_type(1), 
		     const ordinal_type& b_ = ordinal_type(0)) : 
      a(a_), b(b_) {}

    //! Destructor
    ~LinearGrowthRule() {}

    //! Evaluate growth rule
    ordinal_type operator() (const ordinal_type& x) const { return a*x+b; }

  protected:

    //! Slope
    ordinal_type a;

    //! Offset
    ordinal_type b;
  };

  //! A growth rule that always makes the supplied order even
  /*!
   * When used in conjunction with Gaussian quadrature that generates n+1
   * points for a quadrature of order n, this always results in an odd
   * number of points, and thus includes 0.  This allows some nesting
   * in Gaussian-based sparse grids.
   */
  template <typename ordinal_type>
  class EvenGrowthRule {
  public:
    //! Constructor
    EvenGrowthRule() {}

    //! Destructor
    ~EvenGrowthRule() {}

    //! Evaluate growth rule
    ordinal_type operator() (const ordinal_type& x) const { 
      if (x % 2 == 1) return x+1;
      return x;
    }

  };

  //! An exponential growth rule for Clenshaw-Curtis
  template <typename ordinal_type>
  class ClenshawCurtisExponentialGrowthRule {
  public:
    //! Constructor
    ClenshawCurtisExponentialGrowthRule() {}

    //! Destructor
    ~ClenshawCurtisExponentialGrowthRule() {}

    //! Evaluate growth rule
    ordinal_type operator() (const ordinal_type& x) const { 
      if (x == 0) return 0;
      return std::pow(ordinal_type(2),x-1);
    }

  };

  //! An exponential growth rule for Gauss-Patterson
  template <typename ordinal_type>
  class GaussPattersonExponentialGrowthRule {
  public:
    //! Constructor
    GaussPattersonExponentialGrowthRule() {}

    //! Destructor
    ~GaussPattersonExponentialGrowthRule() {}

    //! Evaluate growth rule
    ordinal_type operator() (const ordinal_type& x) const { 
      // Gauss-Patterson rules have precision 3*2*l-1, which is odd.
      // Since discrete orthogonality requires integrating polynomials of
      // order 2*p, setting p = 3*2*{l-1}-1 will yield the largest p such that
      // 2*p <= 3*2^l-1
      if (x == 0) return 0;
      return 3*std::pow(2,x-1)-1;
    }

  };

  /*!
   * \brief An operator for building pseudo-spectral coefficients using
   * a sparse Smolyak construction.
   */
  template <typename factory_type, 
	    typename coeff_compare_type,
	    typename point_compare_type>
  class SparseGridPseudoSpectralOperator {
  public:

    typedef typename factory_type::operator_type operator_type;
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;
    typedef typename operator_type::basis_type basis_type;
    typedef typename operator_type::coeff_type coeff_type;
    typedef typename operator_type::point_type point_type;
    typedef typename operator_type::multiindex_type multiindex_type;

    typedef typename operator_type::domain domain;
    typedef typename operator_type::range range;
    typedef point_compare_type domain_compare;
    typedef coeff_compare_type range_compare;
    typedef std::map<domain,std::pair<value_type,ordinal_type>,domain_compare> domain_set_type;
    typedef std::map<range,ordinal_type,range_compare> range_set_type;
    typedef Teuchos::Array<domain> domain_map_type;
    typedef Teuchos::Array<range> range_map_type;
    typedef Stokhos::LexographicLess<multiindex_type> index_compare;
    typedef std::map<multiindex_type,ordinal_type,index_compare> index_map_type;
    typedef typename domain_set_type::iterator domain_set_iterator;
    typedef typename range_set_type::iterator range_set_iterator;
    typedef typename domain_set_type::const_iterator domain_set_const_iterator;
    typedef typename range_set_type::const_iterator range_set_const_iterator;
    typedef typename domain_map_type::iterator domain_iterator;
    typedef typename domain_map_type::const_iterator domain_const_iterator;
    typedef typename range_map_type::iterator range_iterator;
    typedef typename range_map_type::const_iterator range_const_iterator;

    typedef typename operator_type::domain_iterator op_domain_iterator_type;
    typedef typename operator_type::range_iterator op_range_iterator_type;
    typedef typename operator_type::domain_set_iterator op_domain_set_iterator_type;
    typedef typename operator_type::range_set_iterator op_range_set_iterator_type;
    typedef Teuchos::Array< Teuchos::RCP<operator_type> > operator_set_type;
    
    //! Constructor
    template <typename index_set_type,
	      typename coeff_growth_rule_type,
	      typename point_growth_rule_type>
    SparseGridPseudoSpectralOperator(
      const Teuchos::Array< Teuchos::RCP<const basis_type > >& bases,
      const index_set_type& index_set,
      const coeff_growth_rule_type& coeff_growth_rule,
      const point_growth_rule_type& point_growth_rule,
      const factory_type& operator_factory,
      const coeff_compare_type& coeff_compare = coeff_compare_type(),
      const point_compare_type& point_compare = point_compare_type()) :
      domain_set(point_compare),
      range_set(coeff_compare) {

      // Generate index set for the final Smolyak operator and 
      // corresponding coefficients
      //
      // The Smolyak operator is given by the formula
      //
      // A = \sum_{k\in\K} \bigotimes_{i=1}^d \Delta^i_{k_i}
      //
      // where \Delta^i_0 = 0, \Delta^i_{k_i} = L^i_{k_i} - L^i_{k_i-1},
      // and K is the supplied index set.  This becomes 
      //
      // A = \sum_{k\in\tilde{K}} c(k) \bigotimes_{i=1}^d L^i_{k_i}
      //
      // for some new index set \tilde{K} and coefficient c(k).  Using the
      // formula (cf. G W Wasilkowski and H Wozniakowski, "Explicit cost bounds 
      // of algorithms for multivariate tensor product problems," 
      // Journal of Complexity (11), 1995)
      // 
      // \bigotimes_{i=1}^d \Delta^i_{k_i} = 
      //    \sum_{\alpha\in\Alpha} (-1)^{|\alpha|} 
      //        \bigotimes_{i=1}^d L^i_{k_i-\alpha_i}
      //
      // where \Alpha = {0,1}^d and |\alpha| = \alpha_1 + ... + \alpha_d, we
      // iterate over K and \Alpha, compute k-\alpha and the corresponding
      // coefficient contribution (-1)^{|\alpha|} and store these in a map.  
      // The keys of of this map with non-zero coefficients define
      // \tilde{K} and c(k).
      typedef Stokhos::TensorProductIndexSet<ordinal_type> alpha_set_type; 
      ordinal_type dim = index_set.dimension();
      alpha_set_type alpha_set(dim, 1);
      typename alpha_set_type::iterator alpha_begin = alpha_set.begin();
      typename alpha_set_type::iterator alpha_end = alpha_set.end();
      typename index_set_type::iterator index_iterator = index_set.begin();
      typename index_set_type::iterator index_end = index_set.end();
      multiindex_type diff(dim);
      index_map_type index_map;
      for (; index_iterator != index_end; ++index_iterator) {
        for (typename alpha_set_type::iterator alpha = alpha_begin; 
	     alpha != alpha_end; ++alpha) {
	  bool valid_index = true;
	  for (ordinal_type i=0; i<dim; ++i) {
	    diff[i] = (*index_iterator)[i] - (*alpha)[i];
	    if (diff[i] < 0) {
	      valid_index = false;
	      break;
	    }
	  }
	  if (valid_index) {
	    ordinal_type alpha_order = alpha->order();
	    ordinal_type val;	  
	    if (alpha_order % 2 == 0)
	      val = 1;
	    else
	      val = -1;
	    typename index_map_type::iterator index_map_iterator =
	      index_map.find(diff);
	    if (index_map_iterator == index_map.end())
	      index_map[diff] = val;
	    else
	      index_map_iterator->second += val;
	  }
	}
      }

      // Generate operators, sparse domain and range
      typename index_map_type::iterator index_map_iterator = index_map.begin();
      typename index_map_type::iterator index_map_end = index_map.end();
      for (; index_map_iterator != index_map_end; ++index_map_iterator) {
	
	// Skip indices with zero coefficient
	if (index_map_iterator->second == 0)
	  continue;

	// Apply growth rule to cofficient multi-index
	multiindex_type coeff_growth_index(dim), point_growth_index(dim);
	for (ordinal_type i=0; i<dim; ++i) {
	  coeff_growth_index[i] = 
	    coeff_growth_rule[i](index_map_iterator->first[i]);
	  point_growth_index[i] = 
	    point_growth_rule[i](coeff_growth_index[i]);
	}

	// Build tensor product operator for given index
	Stokhos::TensorProductIndexSet<ordinal_type> tp_index_set(
	  coeff_growth_index);
	Teuchos::RCP<operator_type> op = 
	  operator_factory(tp_index_set, point_growth_index);

	// Get smolyak cofficient for given index
	value_type c = index_map_iterator->second;
	
	// Include points in union over all sets
	op_domain_set_iterator_type op_domain_iterator = op->domain_set_begin();
	op_domain_set_iterator_type op_domain_end = op->domain_set_end();
	for (; op_domain_iterator != op_domain_end; ++op_domain_iterator) {
	  const domain& point = op_domain_iterator->first;
	  value_type w = op_domain_iterator->second.first;
	  domain_set_iterator dsi = domain_set.find(point);
	  if (dsi == domain_set.end())
	    domain_set[point] = std::make_pair(c*w,ordinal_type(0));
	  else
	    dsi->second.first += c*w;
	}

	// Include coefficients in union over all sets
	op_range_iterator_type op_range_iterator = op->range_begin();
	op_range_iterator_type op_range_end = op->range_end();
	for (; op_range_iterator != op_range_end; ++op_range_iterator)
	  range_set[*op_range_iterator] = ordinal_type(0);
      }

      // Generate linear ordering of points
      ordinal_type idx = 0;
      domain_map.resize(domain_set.size());
      for (domain_set_iterator sdi = domain_set.begin();
	   sdi != domain_set.end(); ++sdi) {
	sdi->second.second = idx;
	domain_map[idx] = sdi->first;
	++idx;
      }
 
      // Generate linear odering of coefficients
      idx = 0;
      range_map.resize(range_set.size());
      for (range_set_iterator sri = range_set.begin(); sri != range_set.end(); 
	   ++sri) {
	sri->second = idx;
	range_map[idx] = sri->first;
	++idx;
      }

      // Compute max coefficient orders
      coeff_type coeff_max(dim);
      for (range_set_iterator sri = range_set.begin(); sri != range_set.end(); 
	   ++sri) {
	for (ordinal_type i=0; i<dim; ++i)
	  if (sri->first[i] > coeff_max[i])
	    coeff_max[i] = sri->first[i];
      }
      Teuchos::Array< Teuchos::RCP<const basis_type > > bases2(bases);
      for (ordinal_type i=0; i<dim; i++)
	if (bases[i]->order() < coeff_max[i])
	  bases2[i] = bases[i]->cloneWithOrder(coeff_max[i]);  

      // Generate quadrature operator
      ordinal_type nqp = domain_set.size();
      ordinal_type npc = range_set.size();
      A.reshape(npc,nqp);
      A.putScalar(1.0);
      Teuchos::Array< Teuchos::Array<value_type> > gv(dim);
      for (ordinal_type k=0; k<dim; ++k)
	  gv[k].resize(coeff_max[k]+1);
      for (domain_set_iterator sdi = domain_set.begin();
	   sdi != domain_set.end(); ++sdi) {
	ordinal_type j = sdi->second.second;
	value_type w = sdi->second.first;
	point_type point = sdi->first;
	for (ordinal_type k=0; k<dim; ++k) 
	  bases2[k]->evaluateBases(point[k], gv[k]);
	for (range_set_iterator sri = range_set.begin();
	     sri != range_set.end(); ++sri) {
	  ordinal_type i = sri->second;
	  coeff_type coeff = sri->first;
	  A(i,j) = w;
	  for (ordinal_type k=0; k<dim; k++) {
	    ordinal_type m = coeff[k];
	    A(i,j) *= gv[k][m] / bases2[k]->norm_squared(m);
	  }
	}
      }
      
    }

    //! Destructor
    ~SparseGridPseudoSpectralOperator() {}

    ordinal_type domain_size() const { return domain_set.size(); }
    ordinal_type range_size() const { return range_set.size(); }

    domain_iterator domain_begin() { return domain_map.begin(); }
    domain_iterator domain_end() { return domain_map.end(); }
    range_iterator range_begin() { return range_map.begin(); }
    range_iterator range_end() { return range_map.end(); }

    domain_const_iterator domain_begin() const { return domain_map.begin(); }
    domain_const_iterator domain_end() const { return domain_map.end(); }
    range_const_iterator range_begin() const { return range_map.begin(); }
    range_const_iterator range_end() const { return range_map.end(); }

    domain_set_iterator domain_set_begin() { return domain_set.begin(); }
    domain_set_iterator domain_set_end() { return domain_set.end(); }
    range_set_iterator range_set_begin() { return range_set.begin(); }
    range_set_iterator range_set_end() { return range_set.end(); }

    domain_set_const_iterator domain_set_begin() const { 
      return domain_set.begin(); }
    domain_set_const_iterator domain_set_end() const { 
      return domain_set.end(); }
    range_set_const_iterator range_set_begin() const { 
      return range_set.begin(); }
    range_set_const_iterator range_set_end() const { 
      return range_set.end(); }

    ordinal_type getDomainIndex(const domain& term) const { 
      return domain_set[term];
    }
    ordinal_type getRangeIndex(const range& term) const { 
      return range_set[term];
    }

    const domain& getDomainTerm(ordinal_type n) const {
      return domain_map[n];
    }
    const range& getRangeTerm(ordinal_type n) const {
      return range_map[n];
    }

    //! Apply Smolyak pseudo-spectral operator
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

    //! Smolyak sparse domain
    domain_set_type domain_set;

    //! Smolyak sparse range
    range_set_type range_set;

    //! Map index to domain term
    domain_map_type domain_map;

    //! Map index to range term
    range_map_type range_map;

    //! Matrix mapping point to coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> A;

  };

  /*!
   * \brief An operator for building pseudo-spectral coefficients using
   * a sparse Smolyak construction.
   */
  template <typename factory_type,
	    typename coeff_compare_type,
	    typename point_compare_type>
  class SmolyakPseudoSpectralOperator {
  public:

    typedef typename factory_type::operator_type operator_type;
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;
    typedef typename operator_type::basis_type basis_type;
    typedef typename operator_type::coeff_type coeff_type;
    typedef typename operator_type::point_type point_type;
    typedef typename operator_type::multiindex_type multiindex_type;

    typedef typename operator_type::domain domain;
    typedef typename operator_type::range range;
    typedef point_compare_type domain_compare;
    typedef coeff_compare_type range_compare;
    typedef std::map<domain,std::pair<value_type,ordinal_type>,domain_compare> domain_set_type;
    typedef std::map<range,ordinal_type,range_compare> range_set_type;
    typedef Teuchos::Array<domain> domain_map_type;
    typedef Teuchos::Array<range> range_map_type;
    typedef Stokhos::LexographicLess<multiindex_type> index_compare;
    typedef std::map<multiindex_type,ordinal_type,index_compare> index_map_type;
    typedef typename domain_set_type::iterator domain_set_iterator;
    typedef typename range_set_type::iterator range_set_iterator;
    typedef typename domain_set_type::const_iterator domain_set_const_iterator;
    typedef typename range_set_type::const_iterator range_set_const_iterator;
    typedef typename domain_map_type::iterator domain_iterator;
    typedef typename domain_map_type::const_iterator domain_const_iterator;
    typedef typename range_map_type::iterator range_iterator;
    typedef typename range_map_type::const_iterator range_const_iterator;

    typedef typename operator_type::domain_iterator op_domain_iterator_type;
    typedef typename operator_type::range_iterator op_range_iterator_type;
    typedef typename operator_type::domain_set_iterator op_domain_set_iterator_type;
    typedef typename operator_type::range_set_iterator op_range_set_iterator_type;
    typedef Teuchos::Array< Teuchos::RCP<operator_type> > operator_set_type;
    
    //! Constructor
     template <typename index_set_type,
	      typename coeff_growth_rule_type,
	      typename point_growth_rule_type>
    SmolyakPseudoSpectralOperator(
      const Teuchos::Array< Teuchos::RCP<const basis_type > >& bases,
      const index_set_type& index_set,
      const coeff_growth_rule_type& coeff_growth_rule,
      const point_growth_rule_type& point_growth_rule,
      const factory_type& operator_factory,
      const coeff_compare_type& coeff_compare = coeff_compare_type(),
      const point_compare_type& point_compare = point_compare_type()) :
      domain_set(point_compare),
      range_set(coeff_compare) {

      // Generate index set for the final Smolyak operator and 
      // corresponding coefficients
      //
      // The Smolyak operator is given by the formula
      //
      // A = \sum_{k\in\K} \bigotimes_{i=1}^d \Delta^i_{k_i}
      //
      // where \Delta^i_0 = 0, \Delta^i_{k_i} = L^i_{k_i} - L^i_{k_i-1},
      // and K is the supplied index set.  This becomes 
      //
      // A = \sum_{k\in\tilde{K}} c(k) \bigotimes_{i=1}^d L^i_{k_i}
      //
      // for some new index set \tilde{K} and coefficient c(k).  Using the
      // formula (cf. G W Wasilkowski and H Wozniakowski, "Explicit cost bounds 
      // of algorithms for multivariate tensor product problems," 
      // Journal of Complexity (11), 1995)
      // 
      // \bigotimes_{i=1}^d \Delta^i_{k_i} = 
      //    \sum_{\alpha\in\Alpha} (-1)^{|\alpha|} 
      //        \bigotimes_{i=1}^d L^i_{k_i-\alpha_i}
      //
      // where \Alpha = {0,1}^d and |\alpha| = \alpha_1 + ... + \alpha_d, we
      // iterate over K and \Alpha, compute k-\alpha and the corresponding
      // coefficient contribution (-1)^{|\alpha|} and store these in a map.  
      // The keys of of this map with non-zero coefficients define
      // \tilde{K} and c(k).
      typedef Stokhos::TensorProductIndexSet<ordinal_type> alpha_set_type; 
      ordinal_type dim = index_set.dimension();
      alpha_set_type alpha_set(dim, 1);
      typename alpha_set_type::iterator alpha_begin = alpha_set.begin();
      typename alpha_set_type::iterator alpha_end = alpha_set.end();
      typename index_set_type::iterator index_iterator = index_set.begin();
      typename index_set_type::iterator index_end = index_set.end();
      multiindex_type diff(dim);
      index_map_type index_map;
      for (; index_iterator != index_end; ++index_iterator) {
        for (typename alpha_set_type::iterator alpha = alpha_begin; 
	     alpha != alpha_end; ++alpha) {
	  bool valid_index = true;
	  for (ordinal_type i=0; i<dim; ++i) {
	    diff[i] = (*index_iterator)[i] - (*alpha)[i];
	    if (diff[i] < 0) {
	      valid_index = false;
	      break;
	    }
	  }
	  if (valid_index) {
	    ordinal_type alpha_order = alpha->order();
	    ordinal_type val;	  
	    if (alpha_order % 2 == 0)
	      val = 1;
	    else
	      val = -1;
	    typename index_map_type::iterator index_map_iterator =
	      index_map.find(diff);
	    if (index_map_iterator == index_map.end())
	      index_map[diff] = val;
	    else
	      index_map_iterator->second += val;
	  }
	}
      }

      // Generate operators, sparse domain and range
      typename index_map_type::iterator index_map_iterator = index_map.begin();
      typename index_map_type::iterator index_map_end = index_map.end();
      for (; index_map_iterator != index_map_end; ++index_map_iterator) {

	// Skip indices with zero coefficient
	if (index_map_iterator->second == 0)
	  continue;

        // Apply growth rule to cofficient multi-index
	multiindex_type coeff_growth_index(dim), point_growth_index(dim);
	for (ordinal_type i=0; i<dim; ++i) {
	  coeff_growth_index[i] = 
	    coeff_growth_rule[i](index_map_iterator->first[i]);
	  point_growth_index[i] = 
	    point_growth_rule[i](coeff_growth_index[i]);
	}

	// Build tensor product operator for given index
	Stokhos::TensorProductIndexSet<ordinal_type> tp_index_set(
	  coeff_growth_index);
	Teuchos::RCP<operator_type> op = 
	  operator_factory(tp_index_set, point_growth_index);
	operators.push_back(op);

	// Get smolyak cofficient for given index
	ordinal_type c = index_map_iterator->second;
	smolyak_coeff.push_back(c);
	
	// Include points in union over all sets
	op_domain_set_iterator_type op_domain_iterator = op->domain_set_begin();
	op_domain_set_iterator_type op_domain_end = op->domain_set_end();
	for (; op_domain_iterator != op_domain_end; ++op_domain_iterator) {
	  const domain& point = op_domain_iterator->first;
	  value_type w = op_domain_iterator->second.first;
	  domain_set_iterator dsi = domain_set.find(point);
	  if (dsi == domain_set.end())
	    domain_set[point] = std::make_pair(c*w,ordinal_type(0));
	  else
	    dsi->second.first += c*w;
	}

	// Include coefficients in union over all sets
	op_range_iterator_type op_range_iterator = op->range_begin();
	op_range_iterator_type op_range_end = op->range_end();
	for (; op_range_iterator != op_range_end; ++op_range_iterator)
	  range_set[*op_range_iterator] = ordinal_type(0);
      }

      // Generate linear ordering of points
      ordinal_type idx = 0;
      domain_map.resize(domain_set.size());
      for (domain_set_iterator sdi = domain_set.begin();
	   sdi != domain_set.end(); ++sdi) {
	sdi->second.second = idx;
	domain_map[idx] = sdi->first;
	++idx;
      }
 
      // Generate linear odering of coefficients
      idx = 0;
      range_map.resize(range_set.size());
      for (range_set_iterator sri = range_set.begin(); sri != range_set.end(); 
	   ++sri) {
	sri->second = idx;
	range_map[idx] = sri->first;
	++idx;
      }

      // Build gather/scatter maps into global domain/range for each operator
      gather_maps.resize(operators.size());
      scatter_maps.resize(operators.size());
      for (ordinal_type i=0; i<operators.size(); i++) {
	Teuchos::RCP<operator_type> op = operators[i];
	
	gather_maps[i].reserve(op->domain_size());
	op_domain_iterator_type op_domain_iterator = op->domain_begin();
	op_domain_iterator_type op_domain_end = op->domain_end();
	for (; op_domain_iterator != op_domain_end; ++op_domain_iterator) {
	  gather_maps[i].push_back(domain_set[*op_domain_iterator].second);
	}

	scatter_maps[i].reserve(op->range_size());
	op_range_iterator_type op_range_iterator = op->range_begin();
	op_range_iterator_type op_range_end = op->range_end();
	for (; op_range_iterator != op_range_end; ++op_range_iterator) {
	  scatter_maps[i].push_back(range_set[*op_range_iterator]);
	}
      }
      
    }

    //! Destructor
    ~SmolyakPseudoSpectralOperator() {}

    ordinal_type domain_size() const { return domain_set.size(); }
    ordinal_type range_size() const { return range_set.size(); }

    domain_iterator domain_begin() { return domain_map.begin(); }
    domain_iterator domain_end() { return domain_map.end(); }
    range_iterator range_begin() { return range_map.begin(); }
    range_iterator range_end() { return range_map.end(); }

    domain_const_iterator domain_begin() const { return domain_map.begin(); }
    domain_const_iterator domain_end() const { return domain_map.end(); }
    range_const_iterator range_begin() const { return range_map.begin(); }
    range_const_iterator range_end() const { return range_map.end(); }

    domain_set_iterator domain_set_begin() { return domain_set.begin(); }
    domain_set_iterator domain_set_end() { return domain_set.end(); }
    range_set_iterator range_set_begin() { return range_set.begin(); }
    range_set_iterator range_set_end() { return range_set.end(); }

    domain_set_const_iterator domain_set_begin() const { 
      return domain_set.begin(); }
    domain_set_const_iterator domain_set_end() const { 
      return domain_set.end(); }
    range_set_const_iterator range_set_begin() const { 
      return range_set.begin(); }
    range_set_const_iterator range_set_end() const { 
      return range_set.end(); }

    ordinal_type getDomainIndex(const domain& term) const { 
      return domain_set[term];
    }
    ordinal_type getRangeIndex(const range& term) const { 
      return range_set[term];
    }

    const domain& getDomainTerm(ordinal_type n) const {
      return domain_map[n];
    }
    const range& getRangeTerm(ordinal_type n) const {
      return range_map[n];
    }

    //! Apply Smolyak pseudo-spectral operator
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
      Teuchos::SerialDenseMatrix<ordinal_type,value_type> op_input, op_result;
      result.scale(beta);
      for (ordinal_type i=0; i<operators.size(); i++) {
	Teuchos::RCP<operator_type> op = operators[i];
	if (trans) {
	  op_input.reshape(input.numRows(), op->domain_size());
	  op_result.reshape(result.numRows(), op->range_size());
	}
	else {
	  op_input.reshape(op->domain_size(), input.numCols());
	  op_result.reshape(op->range_size(), result.numCols());
	}
	gather(gather_maps[i], input, trans, op_input);
	op->apply(smolyak_coeff[i], op_input, op_result, 0.0, trans);
	scatter(scatter_maps[i], op_result, trans, result);
      }
    }

  protected:

    void gather(
      const Teuchos::Array<ordinal_type>& map, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
      bool trans, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const {
      if (trans) {
	for (ordinal_type j=0; j<map.size(); j++)
	  for (ordinal_type i=0; i<input.numRows(); i++)
	    result(i,j) = input(i,map[j]);
      }
      else {
	for (ordinal_type j=0; j<input.numCols(); j++)
	  for (ordinal_type i=0; i<map.size(); i++)
	    result(i,j) = input(map[i],j);
      }
    }

    void scatter(
      const Teuchos::Array<ordinal_type>& map, 
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type>& input, 
      bool trans, 
      Teuchos::SerialDenseMatrix<ordinal_type,value_type>& result) const {
      if (trans) {
	for (ordinal_type j=0; j<map.size(); j++)
	  for (ordinal_type i=0; i<input.numRows(); i++)
	    result(i,map[j]) += input(i,j);
      }
      else {
	for (ordinal_type j=0; j<input.numCols(); j++)
	  for (ordinal_type i=0; i<map.size(); i++)
	    result(map[i],j) += input(i,j);
      }
    }

  protected:

    //! Operators comprising smolyak construction
    operator_set_type operators;

    //! Index set map for building Smolyak coefficients
    index_map_type index_map;

    //! Smolyak sparse domain
    domain_set_type domain_set;

    //! Smolyak sparse range
    range_set_type range_set;

    //! Smolyak coefficients
    Teuchos::Array<value_type> smolyak_coeff;

    //! Map index to domain term
    domain_map_type domain_map;

    //! Map index to range term
    range_map_type range_map;

    //! Gather maps for each operator
    Teuchos::Array< Teuchos::Array<ordinal_type> > gather_maps;

    //! Scatter maps for each operator
    Teuchos::Array< Teuchos::Array<ordinal_type> > scatter_maps;

  };

#if 0

  /*!
   * \brief Multivariate orthogonal polynomial basis generated from a
   * Smolyak sparse grid.
   */
  template <typename ordinal_type, typename value_type>
  class SmolyakBasis : 
    public ProductBasis<ordinal_type,value_type> {
  public:

    //! Constructor
    /*!
     * \param bases array of 1-D coordinate bases
     * \param sparse_tol tolerance used to drop terms in sparse triple-product
     *                   tensors
     * \param use_old_cijk_alg use old algorithm for computing the sparse
     *                         triple product tensor  (significantly slower,
     *                         but simpler)
     * \param deriv_coeffs direction used to define derivatives for
     *                     derivative product tensors.  Defaults to
     *                     all one's if not supplied.
     */
    SmolyakBasis(
      const Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type,
 value_type> > >& bases,
      const value_type& sparse_tol = 1.0e-12);

    //! Destructor
    virtual ~SmolyakBasis();

    //! \name Implementation of Stokhos::OrthogPolyBasis methods
    //@{

    //! Return order of basis
    ordinal_type order() const;

    //! Return dimension of basis
    ordinal_type dimension() const;

    //! Return total size of basis
    virtual ordinal_type size() const;

    //! Return array storing norm-squared of each basis polynomial
    /*!
     * Entry \f$l\f$ of returned array is given by \f$\langle\Psi_l^2\rangle\f$
     * for \f$l=0,\dots,P\f$ where \f$P\f$ is size()-1.
     */
    virtual const Teuchos::Array<value_type>& norm_squared() const;

    //! Return norm squared of basis polynomial \c i.
    virtual const value_type& norm_squared(ordinal_type i) const;

    //! Compute triple product tensor
    /*!
     * The \f$(i,j,k)\f$ entry of the tensor \f$C_{ijk}\f$ is given by
     * \f$C_{ijk} = \langle\Psi_i\Psi_j\Psi_k\rangle\f$ where \f$\Psi_l\f$
     * represents basis polynomial \f$l\f$ and \f$i,j=0,\dots,P\f$ where
     * \f$P\f$ is size()-1 and \f$k=0,\dots,p\f$ where \f$p\f$
     * is the supplied \c order.
     */
    virtual 
    Teuchos::RCP< Stokhos::Sparse3Tensor<ordinal_type, value_type> > 
    computeTripleProductTensor(ordinal_type order) const;

    //! Evaluate basis polynomial \c i at zero
    virtual value_type evaluateZero(ordinal_type i) const;

    //! Evaluate basis polynomials at given point \c point
    /*!
     * Size of returned array is given by size(), and coefficients are
     * ordered from order 0 up to size size()-1.
     */
    virtual void evaluateBases(const Teuchos::Array<value_type>& point,
			       Teuchos::Array<value_type>& basis_vals) const;

    //! Print basis to stream \c os
    virtual void print(std::ostream& os) const;

    //! Return string name of basis
    virtual const std::string& getName() const;

    //@}

    //! \name Implementation of Stokhos::ProductBasis methods
    //@{

    //! Get orders of each coordinate polynomial given an index \c i
    /*!
     * The returned array is of size \f$d\f$, where \f$d\f$ is the dimension of
     * the basis, and entry \f$l\f$ is given by \f$i_l\f$ where
     * \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual Teuchos::Array<ordinal_type> getTerm(ordinal_type i) const;

    //! Get index of the multivariate polynomial given orders of each coordinate
    /*!
     * Given the array \c term storing \f$i_1,\dots,\i_d\f$, returns the index
     * \f$i\f$ such that \f$\Psi_i(x) = \psi_{i_1}(x_1)\dots\psi_{i_d}(x_d)\f$.
     */
    virtual ordinal_type 
    getIndex(const Teuchos::Array<ordinal_type>& term) const;

    //! Return coordinate bases
    /*!
     * Array is of size dimension().
     */
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, 
							   value_type> > > 
    getCoordinateBases() const;

    //@}

  private:

    // Prohibit copying
    SmolyakBasis(const SmolyakBasis&);

    // Prohibit Assignment
    SmolyakBasis& operator=(const SmolyakBasis& b);
    
  protected:

    //! Name of basis
    std::string name;

    //! Total order of basis
    ordinal_type p;

    //! Total dimension of basis
    ordinal_type d;

    //! Total size of basis
    ordinal_type sz;

    //! Array of bases
    Teuchos::Array< Teuchos::RCP<const OneDOrthogPolyBasis<ordinal_type, value_type> > > bases;

    //! Array storing order of each basis
    Teuchos::Array<ordinal_type> basis_orders;

    //! Tolerance for computing sparse Cijk
    value_type sparse_tol;

    //! Norms
    Teuchos::Array<value_type> norms;

    //! Temporary array used in basis evaluation
    mutable Teuchos::Array< Teuchos::Array<value_type> > basis_eval_tmp;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

  }; // class SmolyakBasis

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_SmolyakBasisImp.hpp"

#endif 
}

#endif
