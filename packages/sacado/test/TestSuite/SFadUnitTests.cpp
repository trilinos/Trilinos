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

#include "SFadUnitTests.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION( SFadOpsUnitTest );

SFadOpsUnitTest::SFadOpsUnitTest() :
  urand(0.0, 1.0), n(num_comp), tol_a(1.0e-15), tol_r(1.0e-14) {}

SFadOpsUnitTest::SFadOpsUnitTest(int numComponents, double absolute_tolerance, 
				 double relative_tolerance) :
  urand(0.0, 1.0), 
  n(numComponents), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

void SFadOpsUnitTest::setUp() {
  double val;

  val = urand.number();
  a_sfad = SFadType(n,val);
  a_fad = FAD::Fad<double>(n,val);
  
  val = urand.number();
  b_sfad = SFadType(n,val);
  b_fad = FAD::Fad<double>(n,val);

  for (int i=0; i<n; i++) {
    val = urand.number();
    a_sfad.fastAccessDx(i) = val;
    a_fad.fastAccessDx(i) = val;

    val = urand.number();
    b_sfad.fastAccessDx(i) = val;
    b_fad.fastAccessDx(i) = val;
  }
}

void SFadOpsUnitTest::tearDown() {}

void SFadOpsUnitTest::compareFads(const SFadType& x_sfad,
				  const FAD::Fad<double>& x_fad) {

  // Compare sizes
  CPPUNIT_ASSERT(x_sfad.size() == x_fad.size());
  
  // Compare hasFastAccess
  CPPUNIT_ASSERT(x_sfad.hasFastAccess() == x_fad.hasFastAccess());
  
  // Compare values
  compareDoubles(c_sfad.val(), c_fad.val());
  
  for (int i=0; i<x_fad.size(); i++) {
    
    // Compare dx
    compareDoubles(c_sfad.dx(i), c_fad.dx(i));
    
    // Compare fastAccessDx
    compareDoubles(c_sfad.fastAccessDx(i), c_fad.fastAccessDx(i));
  }
}

void SFadOpsUnitTest::compareDoubles(double a, double b) {
  CPPUNIT_ASSERT( fabs(a-b) < tol_a + tol_r*fabs(a) );
}

void SFadOpsUnitTest::testMax() {
  double val;

  SFadType aa_sfad = a_sfad + 1.0;
  c_sfad = max(aa_sfad, a_sfad);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }
  
  c_sfad = max(a_sfad, aa_sfad);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }

  c_sfad = max(a_sfad+1.0, a_sfad);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }
  
  c_sfad = max(a_sfad, a_sfad+1.0);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }
  
  val = a_sfad.val() + 1;
  c_sfad = max(a_sfad, val);
  compareDoubles(c_sfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_sfad.dx(i), 0.0);
  
  val = a_sfad.val() - 1;
  c_sfad = max(a_sfad, val);
  compareDoubles(c_sfad.val(), a_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), a_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), a_sfad.fastAccessDx(i));
  }

  val = b_sfad.val() + 1;
  c_sfad = max(val, b_sfad);
  compareDoubles(c_sfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_sfad.dx(i), 0.0);
  
  val = b_sfad.val() - 1;
  c_sfad = max(val, b_sfad);
  compareDoubles(c_sfad.val(), b_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), b_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), b_sfad.fastAccessDx(i));
  }
}

void SFadOpsUnitTest::testMin() {
  double val;

  SFadType aa_sfad = a_sfad - 1.0;
  c_sfad = min(aa_sfad, a_sfad);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }

  c_sfad = min(a_sfad, aa_sfad);
  compareDoubles(c_sfad.val(), aa_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), aa_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), aa_sfad.fastAccessDx(i));
  }

  val = a_sfad.val() - 1;
  c_sfad = min(a_sfad, val);
  compareDoubles(c_sfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_sfad.dx(i), 0.0);
  
  val = a_sfad.val() + 1;
  c_sfad = min(a_sfad, val);
  compareDoubles(c_sfad.val(), a_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), a_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), a_sfad.fastAccessDx(i));
  }

  val = b_sfad.val() - 1;
  c_sfad = min(val, b_sfad);
  compareDoubles(c_sfad.val(), val);
  for (int i=0; i<n; i++)
    compareDoubles(c_sfad.dx(i), 0.0);
  
  val = b_sfad.val() + 1;
  c_sfad = min(val, b_sfad);
  compareDoubles(c_sfad.val(), b_sfad.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_sfad.dx(i), b_sfad.dx(i));
    compareDoubles(c_sfad.fastAccessDx(i), b_sfad.fastAccessDx(i));
  }
}
