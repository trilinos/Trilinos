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
#include "Stokhos_SmolyakBasis.hpp"

namespace SmolyakBasisUtilsUnitTest {

  // Common setup for unit tests
  template <typename ordinal_type, typename value_type>
  struct UnitTestSetup {
    typedef Stokhos::TensorProductElement<ordinal_type,ordinal_type> coeff_type;
    typedef Stokhos::TensorProductElement<ordinal_type,value_type> point_type;
    typedef Stokhos::TotalOrderLess<coeff_type> coeff_compare;
    typedef Stokhos::FloatingPointLess<value_type> value_compare;
    typedef Stokhos::LexographicLess<point_type,value_compare> point_compare;
    //typedef Stokhos::TensorProductPseudoSpectralOperatorFactory<coeff_compare,point_compare> factory_type;
    typedef Stokhos::TensorProductPseudoSpectralOperatorPSTFactory<ordinal_type,value_type> factory_type;
    typedef Stokhos::SparseGridPseudoSpectralOperator<factory_type,coeff_compare,point_compare> sparse_grid_operator_type;
    typedef Stokhos::SmolyakPseudoSpectralOperator<factory_type,coeff_compare,point_compare> smolyak_operator_type;
    typedef Stokhos::LegendreBasis<ordinal_type,value_type> basis_type;
    typedef Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type> cc_basis_type;
    typedef Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type> gp_basis_type;

    value_type rtol, atol;
    ordinal_type p,d;
    Stokhos::MultiIndex<ordinal_type> upper;
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set;
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases, bases2, cc_bases, cc_bases2, gp_bases, gp_bases2;
    Teuchos::RCP<factory_type> factory, cc_factory, gp_factory;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<ordinal_type,value_type> > basis, cc_basis, gp_basis;
    
    UnitTestSetup() : 
      rtol(1e-12), 
      atol(1e-12), 
      p(4),
      d(2), 
      upper(d, p),
      coeff_index_set(d, p),
      bases(d),
      bases2(d),
      cc_bases(d),
      cc_bases2(d),
      gp_bases(d),
      gp_bases2(d) {
      for (ordinal_type i=0; i<d; i++) {
	bases[i] = Teuchos::rcp(new basis_type(1, true)); // tests implicit resizing within operator
	bases2[i] = Teuchos::rcp(new basis_type(p, true));
	cc_bases[i] = Teuchos::rcp(new cc_basis_type(1, true)); // tests implicit resizing within operator
	cc_bases2[i] = Teuchos::rcp(new cc_basis_type(p, true));
	gp_bases[i] = Teuchos::rcp(new gp_basis_type(1, true)); // tests implicit resizing within operator
	gp_bases2[i] = Teuchos::rcp(new gp_basis_type(p, true));
      }
      factory = Teuchos::rcp(new factory_type(bases));
      cc_factory = Teuchos::rcp(new factory_type(cc_bases));
      gp_factory = Teuchos::rcp(new factory_type(gp_bases));
      basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(bases2));
      cc_basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(cc_bases2));
      gp_basis = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<ordinal_type,value_type>(gp_bases2));
    
    }
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  typedef UnitTestSetup<ordinal_type,value_type> setup_type;
  setup_type setup;

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

  template <typename operator_type>
  bool testCoefficients(const operator_type& smop, Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    ordinal_type coeff_sz = smop.range_size();
    ordinal_type coeff_sz2 = setup.basis->size();

     // Check sizes
    TEUCHOS_TEST_EQUALITY(coeff_sz, coeff_sz2, out, success);
    if (!success)
      return false;

    std::ostream_iterator<ordinal_type> out_ord_iterator(out, " ");

    // Check coefficients
    ordinal_type idx = 0;
    typename operator_type::range_const_iterator ri = smop.range_begin();
    for (; ri != smop.range_end(); ++ri) {
      out << "coeff " << *ri << " == [ ";
      Teuchos::Array<ordinal_type> term = setup.basis->getTerm(idx);
      std::copy(term.begin(), term.end(), out_ord_iterator);
      out << "] : ";
      bool is_equal = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_equal = is_equal && term[j] == (*ri)[j];
      if (is_equal)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }
      ++idx;
    }

    return success;
  }

  template <typename operator_type,
	    typename quadrature_type>
  bool testApply(const operator_type& smop, 
		 const quadrature_type& quad,
		 Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    ordinal_type coeff_sz = smop.range_size();
    ordinal_type point_sz = smop.domain_size();

    // Evaluate function at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(point_sz,2);
    ordinal_type idx = 0;
    typename operator_type::domain_const_iterator pi = smop.domain_begin();
    for (; pi != smop.domain_end(); ++pi) {
      f(idx,0) = quad_func1(pi->getTerm());
      f(idx,1) = quad_func2(pi->getTerm());
      ++idx;
    }

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(coeff_sz,2);
    smop.apply(1.0, f, x, 0.0);

    // Compute PCE cofficients using original approach
    const Teuchos::Array<value_type>& weights = quad.getQuadWeights();
    const Teuchos::Array< Teuchos::Array<value_type> >& points = 
      quad.getQuadPoints();
    const Teuchos::Array< Teuchos::Array<value_type> > & vals = 
      quad.getBasisAtQuadPoints();
    ordinal_type coeff_sz2 = setup.basis->size();
    ordinal_type point_sz2 = weights.size();

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f2(point_sz2,2);
    for (ordinal_type qp=0; qp<point_sz2; qp++) {
      f2(qp,0) = quad_func1(points[qp]);
      f2(qp,1) = quad_func2(points[qp]);
    }

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x2(coeff_sz2,2);
    for (ordinal_type i=0; i<coeff_sz2; i++) {
      for (ordinal_type j=0; j<point_sz2; j++) {
	x2(i,0) += weights[j]*f2(j,0)*vals[j][i];
	x2(i,1) += weights[j]*f2(j,1)*vals[j][i];
      }
      x2(i,0) /= setup.basis->norm_squared(i);
      x2(i,1) /= setup.basis->norm_squared(i);
    }

    // Compare PCE coefficients
    success =  
      Stokhos::compareSDM(x, "x", x2, "x2", 1e-1, 1e-3, out);

    return success;
  }

  template <typename operator_type,
	    typename quadrature_type>
  bool testApplyTrans(const operator_type& smop, 
		      const quadrature_type& quad,
		      Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    ordinal_type coeff_sz = smop.range_size();
    ordinal_type point_sz = smop.domain_size();

    // Evaluate function at quadrature points
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(2,point_sz);
    ordinal_type idx = 0;
    typename operator_type::domain_const_iterator pi = smop.domain_begin();
    for (; pi != smop.domain_end(); ++pi) {
      f(0,idx) = quad_func1(pi->getTerm());
      f(1,idx) = quad_func2(pi->getTerm());
      ++idx;
    }

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(2,coeff_sz);
    smop.apply(1.0, f, x, 0.0, true);

    // Compute PCE cofficients using original approach
    const Teuchos::Array<value_type>& weights = quad.getQuadWeights();
    const Teuchos::Array< Teuchos::Array<value_type> >& points = 
      quad.getQuadPoints();
    const Teuchos::Array< Teuchos::Array<value_type> > & vals = 
      quad.getBasisAtQuadPoints();
    ordinal_type coeff_sz2 = setup.basis->size();
    ordinal_type point_sz2 = weights.size();

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f2(2,point_sz2);
    for (ordinal_type qp=0; qp<point_sz2; qp++) {
      f2(0,qp) = quad_func1(points[qp]);
      f2(1,qp) = quad_func2(points[qp]);
    }

    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x2(2,coeff_sz2);
    for (ordinal_type i=0; i<coeff_sz2; i++) {
      for (ordinal_type j=0; j<point_sz2; j++) {
	x2(0,i) += weights[j]*f2(0,j)*vals[j][i];
	x2(1,i) += weights[j]*f2(1,j)*vals[j][i];
      }
      x2(0,i) /= setup.basis->norm_squared(i);
      x2(1,i) /= setup.basis->norm_squared(i);
    }

    // Compare PCE coefficients
    success =  
      Stokhos::compareSDM(x, "x", x2, "x2", 1e-1, 1e-3, out);

    return success;
  }

  template <typename operator_type,
	    typename quadrature_type>
  bool testSparseGrid(const operator_type& smop, 
		      const quadrature_type& quad,
		      Teuchos::FancyOStream& out) {
    bool success = true;

    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    ordinal_type point_sz = smop.domain_size();
    ordinal_type point_sz2 = quad.size();
    const Teuchos::Array<value_type>& weights = quad.getQuadWeights();
    const Teuchos::Array< Teuchos::Array<value_type> >& points = 
      quad.getQuadPoints();

     // Check sizes
    TEUCHOS_TEST_EQUALITY(point_sz, point_sz2, out, success);
    if (!success)
      return false;

    std::ostream_iterator<value_type> out_val_iterator(out, " ");

    // Check points & weights
    ordinal_type idx = 0;
    typename operator_type::domain_set_const_iterator di = 
      smop.domain_set_begin();
    for (; di != smop.domain_set_end(); ++di) {
      // Check point
      out << "point " << di->first << " == [ ";
      std::copy(points[idx].begin(), points[idx].end(), out_val_iterator);
      out << "] : ";
      bool is_equal = true;
      for (ordinal_type j=0; j<setup.d; j++)
	is_equal = is_equal && points[idx][j] == di->first[j];
      if (is_equal)
	out << "passed" << std::endl;
      else {
	out << "failed" << std::endl;
	success = false; 
      }

      // Check weight
      success = success && 
	Stokhos::compareValues(di->second.first, "w1", weights[idx], "w2", 
			       setup.rtol, setup.atol, out);

      ++idx;
    }

    return success;
  }

  template <typename operator_type>
  bool testDiscreteOrthogonality(const operator_type& smop,
				 Teuchos::FancyOStream& out) {
    typedef typename operator_type::ordinal_type ordinal_type;
    typedef typename operator_type::value_type value_type;

    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > my_bases(setup.bases);

    // Evaluate full basis at all quadrature points
    ordinal_type coeff_sz = smop.range_size();
    ordinal_type point_sz = smop.domain_size();
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> f(point_sz,coeff_sz);
    ordinal_type jdx = 0;
    typename operator_type::range_const_iterator ri = smop.range_begin();
    for (; ri != smop.range_end(); ++ri) {
      const typename setup_type::coeff_type& coeff = *ri;

      // Make sure 1-D bases are large enough
      for (ordinal_type k=0; k<setup.d; ++k)
	if (my_bases[k]->order() < coeff[k])
	  my_bases[k] = my_bases[k]->cloneWithOrder(coeff[k]);

      ordinal_type idx = 0;
      typename operator_type::domain_const_iterator pi = smop.domain_begin();
      for (; pi != smop.domain_end(); ++pi) {
	const typename setup_type::point_type& point = *pi;

	value_type v = 1.0;
	for (ordinal_type k=0; k<setup.d; ++k)
	  v *= my_bases[k]->evaluate(point[k], coeff[k]);
	f(idx,jdx) = v;

	++idx;
      }

      ++jdx;
    }

    // Compute PCE coefficients
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> x(coeff_sz,coeff_sz);
    smop.apply(1.0, f, x, 0.0);

    // Subtract identity
    for (ordinal_type i=0; i<coeff_sz; ++i)
      x(i,i) -= 1.0;

    // Expected answer, which is zero
    Teuchos::SerialDenseMatrix<ordinal_type,value_type> z(coeff_sz,coeff_sz);
    z.putScalar(0.0);
    
    out << "Discrete orthogonality error = " << x.normInf() << std::endl;

    // Compare PCE coefficients
    bool success = Stokhos::compareSDM(x, "x", z, "zero", 1e-14, 1e-14, out);

    return success;
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Coefficients_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);
    
    success = testCoefficients(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Coefficients_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);
    
    success = testCoefficients(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Apply_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testApply(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Apply_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testApply(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, ApplyTrans_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testApplyTrans(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, ApplyTrans_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testApplyTrans(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Grid_Linear ) {
    // Build sparse grid operator with Gaussian abscissas and linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Grid_Linear ) {
    // Build sparse grid operator with Gaussian abscissas linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Grid_ModerateLinear ) {
    // Build sparse grid operator with Gaussian abscissas and moderate
    // linear growth
    Teuchos::Array< Stokhos::LinearGrowthRule<ordinal_type> > coeff_growth(
      setup.d, Stokhos::LinearGrowthRule<ordinal_type>(2,0));
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::MODERATE_RESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Grid_ModerateLinear ) {
    // Build sparse grid operator with Gaussian abscissas and moderate
    // linear growth
    Teuchos::Array< Stokhos::LinearGrowthRule<ordinal_type> > coeff_growth(
      setup.d, Stokhos::LinearGrowthRule<ordinal_type>(2,0));
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.basis, setup.p, 1e-12, Pecos::MODERATE_RESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Grid_ClenshawCurtis ) {
    // Build sparse grid operator with Clenshaw-Curtis abscissas and
    // exponential growth
    Teuchos::Array< Stokhos::ClenshawCurtisExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.cc_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.cc_factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.cc_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Grid_ClenshawCurtis ) {
    // Build sparse grid operator with Clenshaw-Curtis abscissas and
    // exponential growth
    Teuchos::Array< Stokhos::ClenshawCurtisExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.cc_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.cc_factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.cc_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, Grid_GaussPatterson ) {
    // Build sparse grid operator with Gauss-Patterson abscissas and
    // exponential growth
    Teuchos::Array< Stokhos::GaussPattersonExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.gp_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.gp_factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.gp_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, Grid_GaussPatterson ) {
    // Build sparse grid operator with Gauss-Patterson abscissas and
    // exponential growth
    Teuchos::Array< Stokhos::GaussPattersonExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.gp_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.gp_factory);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      setup.gp_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    
    success = testSparseGrid(op, quad, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridOperator, 
		     DiscreteOrthogonality_Linear ) {
    // Build sparse grid operator with Gaussian abscissas linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::sparse_grid_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    success = testDiscreteOrthogonality(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, 
		     DiscreteOrthogonality_Linear ) {
    // Build sparse grid operator with Gaussian abscissas linear growth
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    success = testDiscreteOrthogonality(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, 
		     DiscreteOrthogonality_ModerateLinear ) {
    // Build sparse grid operator with Gaussian abscissas and moderate
    // linear growth
    Teuchos::Array< Stokhos::LinearGrowthRule<ordinal_type> > coeff_growth(
      setup.d, Stokhos::LinearGrowthRule<ordinal_type>(2,0));
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.factory);

    success = testDiscreteOrthogonality(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, 
		     DiscreteOrthogonality_ClenshawCurtis ) {
    // Build sparse grid operator with Clenshaw-Curtis abscissas and 
    // exponential growth
    Teuchos::Array< Stokhos::ClenshawCurtisExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.cc_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.cc_factory);

    success = testDiscreteOrthogonality(op, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakOperator, DiscreteOrthogonality_GaussPatterson ) {
    // Build sparse grid operator with Gauss-Patterson abscissas and
    // exponential growth
    Teuchos::Array< Stokhos::GaussPattersonExponentialGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    typename setup_type::smolyak_operator_type op(
      setup.gp_bases, setup.coeff_index_set, coeff_growth, point_growth, 
      *setup.gp_factory);

    success = testDiscreteOrthogonality(op, out);
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
