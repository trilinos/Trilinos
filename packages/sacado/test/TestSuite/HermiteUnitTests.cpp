// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "HermiteUnitTests.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION( HermiteUnitTest );

HermiteUnitTest::HermiteUnitTest() :
  urand(0.0, 1.0), tol_a(1.0e-15), tol_r(1.0e-14) {
  std::vector< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(1); 
  bases[0] = Teuchos::rcp(new basis_type(0));
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
  Teuchos::RCP<exp_type> expansion = Teuchos::rcp(new exp_type(basis));
  pce_type::initExpansion(expansion);
}

HermiteUnitTest::HermiteUnitTest(double absolute_tolerance, 
				 double relative_tolerance) :
  urand(0.0, 1.0), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

void HermiteUnitTest::setUp() {
  double val;

  val = urand.number();
  ac = pce_type(val);
  a = val;
  
  val = urand.number();
  bc = pce_type(val);
  b = val;

  cc = pce_type(1.123);
  c = 1.123;
}

void HermiteUnitTest::tearDown() {}

void HermiteUnitTest::comparePCEs(const pce_type& xc, double x) {

  // Compare sizes
  CPPUNIT_ASSERT(xc.size() == 1);
  
  // Compare hasFastAccess
  CPPUNIT_ASSERT(xc.hasFastAccess(0) == true);
  
  // Compare values
  compareDoubles(xc.coeff(0), x);
}

void HermiteUnitTest::compareDoubles(double x, double y) {
  CPPUNIT_ASSERT( fabs(x-y) < tol_a + tol_r*fabs(x) );
}
