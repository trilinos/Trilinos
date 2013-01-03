// $Id$ 
// $Source$ 
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

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ArrayView.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include <iterator>

namespace ProductBasisUtilsUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType sz,p,d;
    
    UnitTestSetup() {
      rtol = 1e-12;
      atol = 1e-12;
      d = 3;
      p = 5;
    }
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  UnitTestSetup<ordinal_type,value_type> setup;

  // Utility function for computing factorials
  template <typename ordinal_type>
  ordinal_type factorial(const ordinal_type& n) {
    ordinal_type res = 1;
    for (ordinal_type i=1; i<=n; ++i)
      res *= i;
    return res;
  }

  // Function for testing quadratures
  template <typename scalar_type>
  scalar_type quad_func1(const Teuchos::Array<scalar_type>& x) {
    scalar_type val = 0.0;
    for (int i=0; i<x.size(); i++)
      val += x[i];
    return std::exp(val);
  }

  // Function for testing quadratures
  template <typename scalar_type>
  scalar_type quad_func2(const Teuchos::Array<scalar_type>& x) {
    scalar_type val = 0.0;
    for (int i=0; i<x.size(); i++)
      val += x[i];
    return std::sin(val);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, NChooseK ) {
    ordinal_type n, k, v1, v2;

    success = true;

    n = 7; k = 3;  // n-k > k
    v1 = Stokhos::n_choose_k(n,k);
    v2 = factorial(n)/(factorial(k)*factorial(n-k));
    if (v1 != v2) {
      out << "For n  =" << n << ", k = " << k << ", n_choose_k = " << v1 
	  << " != " << v2 << std::endl;
      success = false;
    }

    n = 7; k = 4; // n-k < k
    v1 = Stokhos::n_choose_k(n,k);
    v2 = factorial(n)/(factorial(k)*factorial(n-k));
    if (v1 != v2) {
      out << "For n  =" << n << ", k = " << k << ", n_choose_k = " << v1 
	  << " != " << v2 << std::endl;
      success = false;
    }
    
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderLess ) {
    success = true;

    // Build sorted index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef Stokhos::TotalOrderLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    typedef multiindex_set::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());
    
    // Print sorted index set
    std::ostream_iterator<multiindex_type> out_iterator(out, "\n");
    out << std::endl << "Sorted total order index set (dimension = " << setup.d
	<< ", order = " << setup.p << "):" << std::endl;
    std::copy(sortedIndexSet.begin(), sortedIndexSet.end(), out_iterator);  

    // Ensure orders of each index are increasing
    iterator prev = sortedIndexSet.begin();
    iterator curr = prev; ++curr;
    while (curr != sortedIndexSet.end()) {
      ordinal_type order_prev = prev->order();
      ordinal_type order_curr = curr->order();
      ordinal_type i = 0;
      while (i < setup.d && order_prev == order_curr) {
	order_prev -= (*prev)[i];
	order_curr -= (*curr)[i];
	++i;
      }
      if (order_prev >= order_curr) {
	out << "previous index " << *prev << " and current index "
	    << *curr << " are out of order" << std::endl;
	success = false;
      }
      prev = curr;
      ++curr;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, LexographicLess ) {
    success = true;

    // Build sorted index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef Stokhos::LexographicLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    typedef multiindex_set::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());
    
    // Print sorted index set
    std::ostream_iterator<multiindex_type> out_iterator(out, "\n");
    out << std::endl << "Sorted lexicographic index set (dimension = " 
	<< setup.d << ", order = " << setup.p << "):" << std::endl;
    std::copy(sortedIndexSet.begin(), sortedIndexSet.end(), out_iterator);  

    // Ensure orders of each index are increasing
    iterator prev = sortedIndexSet.begin();
    iterator curr = prev; ++curr;
    while (curr != sortedIndexSet.end()) {
      ordinal_type i = 0;
      while (i < setup.d && (*prev)[i] == (*curr)[i]) ++i;
      if (i == setup.d || (*prev)[i] >= (*curr)[i]) {
	out << "previous index " << *prev << " and current index "
	    << *curr << " are out of order" << std::endl;
	success = false;
      }
      prev = curr;
      ++curr;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, FloatingPointLess ) {
    success = true;

    value_type tol=1e-12;
    Stokhos::FloatingPointLess<value_type> less(tol);

    TEUCHOS_TEST_EQUALITY(less(-0.774597,-0.774597), false, out, success);
    TEUCHOS_TEST_EQUALITY(less(-0.774597+tol/2.0,-0.774597), false, out, success);
    TEUCHOS_TEST_EQUALITY(less(-0.774597-tol/2.0,-0.774597), false, out, success);
    TEUCHOS_TEST_EQUALITY(less(-0.774597,-0.774597+tol/2.0), false, out, success);
    TEUCHOS_TEST_EQUALITY(less(-0.774597,-0.774597-tol/2.0), false, out, success);
    TEUCHOS_TEST_EQUALITY(less(-0.774597,0.0), true, out, success);
    TEUCHOS_TEST_EQUALITY(less(0.0,-0.774597), false, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, LexographicFloatingPointLess ) {
    success = true;
    
    typedef Stokhos::TensorProductElement<ordinal_type,value_type> term_type;
    typedef Stokhos::FloatingPointLess<value_type> comp_type;
    term_type a(2), b(2);
    Stokhos::LexographicLess<term_type,comp_type> less;
    a[0] = -0.774597; a[1] = -0.774597;
    b[0] = -0.774597; b[1] = 0.0;

    TEUCHOS_TEST_EQUALITY(less(a,b), true, out, success);
    TEUCHOS_TEST_EQUALITY(less(b,a), false, out, success);
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderIndexSet ) {
    success = true;

    // Build index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);

    // Print index set
    out << std::endl << "Total order index set (dimension = " << setup.d
	<< ", order = " << setup.p << "):" << std::endl;
    std::ostream_iterator<multiindex_type> out_iterator(out, "\n");
    std::copy(indexSet.begin(), indexSet.end(), out_iterator);

    // Verify each index lies appropriatly in the set
    for (iterator i=indexSet.begin(); i!=indexSet.end(); ++i) {
      if (i->order() < 0 || i->order() > setup.p) {
	out << "index " << *i << " does not lie in total order set! "
	    << std::endl;
	success = false; 
      }
    }
    
    // Put indices in sorted container -- this will ensure there are no
    // duplicates, if we get the right size
    typedef Stokhos::TotalOrderLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());

    out << "sorted index set size = " << sortedIndexSet.size() << std::endl;
    out << "expected index set size = " 
	<< Stokhos::n_choose_k(setup.p+setup.d,setup.d) << std::endl;
    if (static_cast<ordinal_type>(sortedIndexSet.size()) != 
	Stokhos::n_choose_k(setup.p+setup.d,setup.d))
      success = false;
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, 
		     AnisotropicTotalOrderIndexSet ) {
    success = true;

    // Build index set of dimension d and order p
    typedef Stokhos::AnisotropicTotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    multiindex_type upper(setup.d);
    for (ordinal_type i=0; i<setup.d; ++i)
      upper[i] = i+1;
    index_set_type indexSet(setup.p, upper);

    // Print index set
    out << std::endl << "Anisotropic total order index set (dimension = " 
	<< setup.d << ", order = " << setup.p << ", and component orders = "
	<< upper << "):" << std::endl;
    std::ostream_iterator<multiindex_type> out_iterator(out, "\n");
    std::copy(indexSet.begin(), indexSet.end(), out_iterator);

    // Verify each index lies appropriatly in the set
    for (iterator i=indexSet.begin(); i!=indexSet.end(); ++i) {
      if (i->order() < 0 || i->order() > setup.p || !i->termWiseLEQ(upper)) {
	out << "index " << *i << " does not lie in total order set! "
	    << std::endl;
	success = false; 
      }
    }
    
    // Need to figure out how to compute the size of such an index set
    /*
    // Put indices in sorted container -- this will ensure there are no
    // duplicates, if we get the right size
    typedef Stokhos::TotalOrderLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());

    out << "sorted index set size = " << sortedIndexSet.size() << std::endl;
    out << "expected index set size = " 
	<< Stokhos::n_choose_k(setup.p+setup.d,setup.d) << std::endl;
    if (static_cast<ordinal_type>(sortedIndexSet.size()) != 
	Stokhos::n_choose_k(setup.p+setup.d,setup.d))
      success = false;
    */
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderBasis ) {
    success = true;

    // Build index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);

    // Build total-order basis from index set
    typedef Stokhos::TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef Stokhos::TotalOrderLess<coeff_type> less_type;
    typedef std::map<coeff_type, ordinal_type, less_type> basis_set_type;
    typedef basis_set_type::iterator basis_set_iterator;
    typedef Teuchos::Array<coeff_type> basis_map_type;
    basis_set_type basis_set;
    basis_map_type basis_map;
    Stokhos::ProductBasisUtils::buildProductBasis(
      indexSet, basis_set, basis_map);

    // Build total-order basis directly
    ordinal_type sz;
    Teuchos::Array< Stokhos::MultiIndex<ordinal_type> > terms;
    Teuchos::Array<ordinal_type> num_terms;
    Stokhos::CompletePolynomialBasisUtils<ordinal_type,value_type>::
      compute_terms(setup.p, setup.d, sz, terms, num_terms);

    // Check sizes
    TEUCHOS_TEST_EQUALITY(static_cast<ordinal_type>(basis_set.size()), 
			  static_cast<ordinal_type>(basis_map.size()), 
			  out, success);
    TEUCHOS_TEST_EQUALITY(static_cast<ordinal_type>(basis_set.size()), 
			  static_cast<ordinal_type>(terms.size()), 
			  out, success);
    TEUCHOS_TEST_EQUALITY(sz, static_cast<ordinal_type>(terms.size()), 
			  out, success);

    std::ostream_iterator<ordinal_type> out_iterator(out, " ");
    for (ordinal_type i=0; i<sz; i++) {

      // Verify terms match
      out << "term " << basis_map[i] << " == " << terms[i] << " : ";
      bool is_equal = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_equal = is_equal && terms[i][j] == basis_map[i][j];
      if (is_equal)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }

      // Verify global index mapping matches
      TEUCHOS_TEST_EQUALITY(basis_set[basis_map[i]], i, out, success);
    }

  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TensorProductBasis ) {
    success = true;

    // Build index set of dimension d and order p
    typedef Stokhos::TensorProductIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);

    // Build total-order basis from index set
    typedef Stokhos::TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef Stokhos::TotalOrderLess<coeff_type> less_type;
    typedef std::map<coeff_type, ordinal_type, less_type> basis_set_type;
    typedef basis_set_type::iterator basis_set_iterator;
    typedef Teuchos::Array<coeff_type> basis_map_type;
    basis_set_type basis_set;
    basis_map_type basis_map;
    Stokhos::ProductBasisUtils::buildProductBasis(
      indexSet, basis_set, basis_map);

    // Compute expected size
    ordinal_type sz = 1;
    for (ordinal_type i=0; i<setup.d; ++i)
      sz *= setup.p+1;

    // Check sizes
    TEUCHOS_TEST_EQUALITY(static_cast<ordinal_type>(basis_set.size()), 
			  static_cast<ordinal_type>(basis_map.size()), 
			  out, success);
    TEUCHOS_TEST_EQUALITY(sz, static_cast<ordinal_type>(basis_set.size()), 
			  out, success);

    std::ostream_iterator<ordinal_type> out_iterator(out, " ");
    for (ordinal_type i=0; i<sz; i++) {

      // Verify terms match
      out << "term " << basis_map[i] << " <= " << setup.p << " : ";
      bool is_less = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_less = is_less && basis_map[i][j] <= setup.p;
      if (is_less)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }

      // Verify global index mapping matches
      TEUCHOS_TEST_EQUALITY(basis_set[basis_map[i]], i, out, success);
    }

  }

  template <typename ordinal_type>
  struct total_order_predicate {
    ordinal_type dim, order;

    total_order_predicate(ordinal_type dim_, ordinal_type order_) :
      dim(dim_), order(order_) {}

    template <typename term_type>
    bool operator() (const term_type& term) const {
      ordinal_type sum = 0;
      for (ordinal_type i=0; i<dim; ++i)
	sum += term[i];
      return sum <= order;
    }

  };

  template <typename basis_set_type>
  struct general_predicate {
    const basis_set_type& basis_set;

    general_predicate(const basis_set_type& basis_set_) :
      basis_set(basis_set_) {}

    template <typename term_type>
    bool operator() (const term_type& term) const {
      return basis_set.find(term) != basis_set.end();
    }

  };


  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderSparse3Tensor ) {
    success = true;
    ordinal_type dim = setup.d;
    ordinal_type order = setup.p;

    // Build index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    index_set_type indexSet(dim, 0, order);

    // Build total-order basis from index set
    typedef Stokhos::TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef Stokhos::TotalOrderLess<coeff_type> less_type;
    typedef std::map<coeff_type, ordinal_type, less_type> basis_set_type;
    typedef basis_set_type::iterator basis_set_iterator;
    typedef Teuchos::Array<coeff_type> basis_map_type;
    basis_set_type basis_set;
    basis_map_type basis_map;
    Stokhos::ProductBasisUtils::buildProductBasis(
      indexSet, basis_set, basis_map);

    // 1-D bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dim);
    for (ordinal_type i=0; i<dim; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(order, true));

    // Build Cijk tensor
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;
    total_order_predicate<ordinal_type> pred(dim, order);
    //general_predicate<basis_set_type> pred(basis_set);
    Teuchos::RCP<Cijk_type> Cijk =
      Stokhos::ProductBasisUtils::computeTripleProductTensor<ordinal_type,value_type>(bases, basis_set, basis_map, pred, pred);

    // Build Cijk tensor using original approach
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases));
    Teuchos::RCP<Cijk_type> Cijk2 =
      basis->computeTripleProductTensor(basis->size());
    
    // Check sizes
    TEUCHOS_TEST_EQUALITY(Cijk->num_k(), Cijk2->num_k(), out, success);
    TEUCHOS_TEST_EQUALITY(Cijk->num_entries(), Cijk2->num_entries(), out, success);
    
    // Check tensors match
    for (Cijk_type::k_iterator k_it=Cijk2->k_begin(); 
	 k_it!=Cijk2->k_end(); ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = Cijk2->j_begin(k_it); 
	   j_it != Cijk2->j_end(k_it); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk2->i_begin(j_it);
	     i_it != Cijk2->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = Cijk->getValue(i,j,k);
	  double c2 = Stokhos::value(i_it);
	  double tol = setup.atol + c2*setup.rtol;
	  double err = std::abs(c-c2);
	  bool s = err < tol;
	  if (!s) {
	    out << std::endl
		<< "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
		<< " = " << "rel_err( " << c << ", " << c2 << " ) = " << err 
		<< " <= " << tol << " : ";
	    if (s) out << "Passed.";
	    else out << "Failed!";
	    out << std::endl;
	  }
	  success = success && s;
	}
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderSparse3TensorNew ) {
    success = true;
    // ordinal_type dim = setup.d;
    // ordinal_type order = setup.p;
    ordinal_type dim = 2;
    ordinal_type order = 3;

    // Build index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef index_set_type::iterator iterator;
    index_set_type indexSet(dim, 0, order);

    // Build total-order basis from index set
    typedef Stokhos::TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef Stokhos::TotalOrderLess<coeff_type> less_type;
    typedef std::map<coeff_type, ordinal_type, less_type> basis_set_type;
    typedef basis_set_type::iterator basis_set_iterator;
    typedef Teuchos::Array<coeff_type> basis_map_type;
    basis_set_type basis_set;
    basis_map_type basis_map;
    Stokhos::ProductBasisUtils::buildProductBasis(
      indexSet, basis_set, basis_map);

    // 1-D bases
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(dim);
    for (ordinal_type i=0; i<dim; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(order, true));

    // Build Cijk tensor
    typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk_type;
    total_order_predicate<ordinal_type> pred(dim, order);
    //general_predicate<basis_set_type> pred(basis_set);
    Teuchos::RCP<Cijk_type> Cijk =
      Stokhos::ProductBasisUtils::computeTripleProductTensorNew<ordinal_type,value_type>(bases, basis_set, basis_map, pred, pred, false);

    // Build Cijk tensor using original approach
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases));
    Teuchos::RCP<Cijk_type> Cijk2 =
      basis->computeTripleProductTensor(basis->size());

    std::cout << *Cijk <<std::endl;
    std::cout << *Cijk2 <<std::endl;
    
    // Check sizes
    TEUCHOS_TEST_EQUALITY(Cijk->num_k(), Cijk2->num_k(), out, success);
    TEUCHOS_TEST_EQUALITY(Cijk->num_entries(), Cijk2->num_entries(), out, success);
    
    // Check tensors match
    for (Cijk_type::k_iterator k_it=Cijk2->k_begin(); 
	 k_it!=Cijk2->k_end(); ++k_it) {
      int k = Stokhos::index(k_it);
      for (Cijk_type::kj_iterator j_it = Cijk2->j_begin(k_it); 
	   j_it != Cijk2->j_end(k_it); ++j_it) {
	int j = Stokhos::index(j_it);
	for (Cijk_type::kji_iterator i_it = Cijk2->i_begin(j_it);
	     i_it != Cijk2->i_end(j_it); ++i_it) {
	  int i = Stokhos::index(i_it);
	  double c = Cijk->getValue(i,j,k);
	  double c2 = Stokhos::value(i_it);
	  double tol = setup.atol + c2*setup.rtol;
	  double err = std::abs(c-c2);
	  bool s = err < tol;
	  if (!s) {
	    out << std::endl
		<< "Check: rel_err( C(" << i << "," << j << "," << k << ") )"
		<< " = " << "rel_err( " << c << ", " << c2 << " ) = " << err 
		<< " <= " << tol << " : ";
	    if (s) out << "Passed.";
	    else out << "Failed!";
	    out << std::endl;
	  }
	  success = success && s;
	}
      }
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, TotalOrderMapping ) {
    success = true;

    // Build sorted index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef Stokhos::TotalOrderLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    typedef multiindex_set::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());

    // Loop over each index set element and test if mapping 
    // computes proper global index
    iterator i = sortedIndexSet.begin(); 
    ordinal_type idx = 0;
    while (i != sortedIndexSet.end()) {
      ordinal_type idx_mapping = totalOrderMapping(*i);
      out << *i << ":  index = " << idx << " mapped index = " << idx_mapping 
	  << ":  ";
      if (idx == idx_mapping)
	out << "passed";
      else {
	out << "failed";
	success = false;
      }
      out << std::endl;
      ++i;
      ++idx;
    }
  }

  TEUCHOS_UNIT_TEST( Stokhos_ProductBasisUtils, LexographicMapping ) {
    success = true;

    // Build sorted index set of dimension d and order p
    typedef Stokhos::TotalOrderIndexSet<ordinal_type> index_set_type;
    typedef index_set_type::multiindex_type multiindex_type;
    typedef Stokhos::LexographicLess<multiindex_type> less_type;
    typedef std::set<multiindex_type, less_type> multiindex_set;
    typedef multiindex_set::iterator iterator;
    index_set_type indexSet(setup.d, 0, setup.p);
    multiindex_set sortedIndexSet(indexSet.begin(), indexSet.end());

    // Loop over each index set element and test if mapping 
    // computes proper global index
    iterator i = sortedIndexSet.begin(); 
    ordinal_type idx = 0;
    while (i != sortedIndexSet.end()) {
      ordinal_type idx_mapping = lexicographicMapping(*i,setup.p);
      out << *i << ":  index = " << idx << " mapped index = " << idx_mapping 
	  << ":  ";
      if (idx == idx_mapping)
	out << "passed";
      else {
	out << "failed";
	success = false;
      }
      out << std::endl;
      ++i;
      ++idx;
    }
  }
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  Teuchos::TimeMonitor::summarize(std::cout);
  return res;
}
