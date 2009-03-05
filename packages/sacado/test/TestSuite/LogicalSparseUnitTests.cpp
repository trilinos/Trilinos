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

#include "LogicalSparseUnitTests.hpp"

CPPUNIT_TEST_SUITE_REGISTRATION( LogicalSparseOpsUnitTest );

LogicalSparseOpsUnitTest::LogicalSparseOpsUnitTest() :
  urand(0.0, 1.0), n(5), tol_a(1.0e-15), tol_r(1.0e-14) {}

LogicalSparseOpsUnitTest::LogicalSparseOpsUnitTest(int numComponents, 
						   double absolute_tolerance, 
						   double relative_tolerance) :
  urand(0.0, 1.0), 
  n(numComponents), 
  tol_a(absolute_tolerance), 
  tol_r(relative_tolerance) {}

void LogicalSparseOpsUnitTest::setUp() {
  double val;

  val = urand.number();
  a_dfad = DFadType(n,val);
  a_ls = LSType(n,val);
  
  val = urand.number();
  b_dfad = DFadType(n,val);
  b_ls = LSType(n,val);

  val = urand.number();
  c_dfad = val;
  c_ls = val;

  for (int i=0; i<n; i++) {
    val = urand.number();
    a_dfad.fastAccessDx(i) = val;
    a_ls.fastAccessDx(i) = 1;

    val = urand.number();
    b_dfad.fastAccessDx(i) = val;
    b_ls.fastAccessDx(i) = 1;
  }
}

void LogicalSparseOpsUnitTest::tearDown() {}

void LogicalSparseOpsUnitTest::compareFads(const DFadType& x_dfad,
					   const LSType& x_ls) {

  // Compare sizes
  CPPUNIT_ASSERT(x_dfad.size() == x_ls.size());
  
  // Compare hasFastAccess
  CPPUNIT_ASSERT(x_dfad.hasFastAccess() == x_ls.hasFastAccess());
  
  // Compare values
  compareDoubles(x_dfad.val(), x_ls.val());
  
  for (int i=0; i<x_ls.size(); i++) {
    
    // Compare dx
    compareDx(x_dfad.dx(i), x_ls.dx(i));
    
    // Compare fastAccessDx
    compareDx(x_dfad.fastAccessDx(i), x_ls.fastAccessDx(i));
  }
}

void LogicalSparseOpsUnitTest::compareDoubles(double a, double b) {
  CPPUNIT_ASSERT( fabs(a-b) < tol_a + tol_r*fabs(a) );
}

void LogicalSparseOpsUnitTest::compareBools(bool a, bool b) {
  CPPUNIT_ASSERT( a == b );
}

void LogicalSparseOpsUnitTest::compareDx(double a, bool b) {
  CPPUNIT_ASSERT( (a && b) || !(a || b) );
}

void LogicalSparseOpsUnitTest::testMax() {
  double val;

  LSType aa_ls = a_ls + 1.0;
  c_ls = max(aa_ls, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  
  c_ls = max(a_ls, aa_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  c_ls = max(a_ls+1.0, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  
  c_ls = max(a_ls, a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (same)
  this->c_ls = max(this->a_ls+1.0, this->a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (different)
  this->c_ls = max(this->a_ls+1.0, this->a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  this->c_ls = max(this->a_ls-1.0, this->a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  
  val = a_ls.val() + 1;
  c_ls = max(a_ls, val);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  
  val = a_ls.val() - 1;
  c_ls = max(a_ls, val);
  compareDoubles(c_ls.val(), a_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), a_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), a_ls.fastAccessDx(i));
  }

  val = b_ls.val() + 1;
  c_ls = max(val, b_ls);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  
  val = b_ls.val() - 1;
  c_ls = max(val, b_ls);
  compareDoubles(c_ls.val(), b_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), b_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), b_ls.fastAccessDx(i));
  }
}

void LogicalSparseOpsUnitTest::testMin() {
  double val;

  LSType aa_ls = a_ls - 1.0;
  c_ls = min(aa_ls, a_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  c_ls = min(a_ls, aa_ls);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), aa_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

    // Expr, Expr (same)
  this->c_ls = min(this->a_ls-1.0, this->a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  // Expr, Expr (different)
  this->c_ls = min(this->a_ls+1.0, this->a_ls-1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }
  this->c_ls = min(this->a_ls-1.0, this->a_ls+1.0);
  compareDoubles(c_ls.val(), aa_ls.val());
  for (int i=0; i<n; i++) {
    compareDoubles(c_ls.dx(i), aa_ls.dx(i));
    compareDoubles(c_ls.fastAccessDx(i), aa_ls.fastAccessDx(i));
  }

  val = a_ls.val() - 1;
  c_ls = min(a_ls, val);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  
  val = a_ls.val() + 1;
  c_ls = min(a_ls, val);
  compareDoubles(c_ls.val(), a_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), a_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), a_ls.fastAccessDx(i));
  }

  val = b_ls.val() - 1;
  c_ls = min(val, b_ls);
  compareDoubles(c_ls.val(), val);
  for (int i=0; i<n; i++)
    compareBools(c_ls.dx(i), 0);
  
  val = b_ls.val() + 1;
  c_ls = min(val, b_ls);
  compareDoubles(c_ls.val(), b_ls.val());
  for (int i=0; i<n; i++) {
    compareBools(c_ls.dx(i), b_ls.dx(i));
    compareBools(c_ls.fastAccessDx(i), b_ls.fastAccessDx(i));
  }
}
