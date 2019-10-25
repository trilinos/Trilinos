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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado_No_Kokkos.hpp"

#if SACADO_ENABLE_NEW_DESIGN
#include "Sacado_Fad_Exp_DMFad.hpp"
#include "Sacado_Fad_Exp_DVFad.hpp"
#include "Sacado_Fad_Exp_ViewFad.hpp"

// Size used for all Fad types
const int global_fad_size = 10;

//
// Move constructor tests
//

TEUCHOS_UNIT_TEST( MoveConstructorTests, SFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::SFad<value_type,global_fad_size> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, SLFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::SLFad<value_type,global_fad_size> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, DFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  value_type *null_ptr = nullptr;
  TEST_EQUALITY_CONST( x.size(), 0 );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY( x.dx(), null_ptr );
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, DMFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DMFad<value_type> ad_type;
  success = true;

  // Setup memory pool
  Sacado::Fad::MemPoolManager<value_type> poolManager(100);
  Sacado::Fad::MemPool *pool = poolManager.getMemoryPool(global_fad_size);
  ad_type::setDefaultPool(pool);

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  value_type *null_ptr = nullptr;
  TEST_EQUALITY_CONST( x.size(), 0 );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY_CONST( x.dx(), null_ptr );
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, DVFad_Owned )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DVFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, DVFad_Unowned )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DVFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  value_type x_val = value_type(1.5);
  std::vector<value_type> x_dx(global_fad_size);
  for (int i=0; i<global_fad_size; ++i)
    x_dx[i] = value_type(2.0+i);
  ad_type x(global_fad_size, &x_val, x_dx.data(), 1, false);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveConstructorTests, ViewFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DFad<value_type> dfad_type;
  typedef Sacado::Fad::Exp::ViewFad<value_type,0,1,dfad_type> ad_type;
  success = true;

  // Initialize AD type
  value_type x_val = value_type(1.5);
  std::vector<value_type> x_dx(global_fad_size);
  for (int i=0; i<global_fad_size; ++i)
    x_dx[i] = value_type(2.0+i);
  ad_type x(x_dx.data(), &x_val, global_fad_size, 1);

  // Move x into y
  ad_type y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  TEST_EQUALITY( y.dx(), x_dx.data() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY( x.dx(), x_dx.data() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

//
// Move assignment tests
//

TEUCHOS_UNIT_TEST( MoveAssignmentTests, SFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::SFad<value_type,global_fad_size> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, SLFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::SLFad<value_type,global_fad_size> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, DFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  value_type *null_ptr = nullptr;
  TEST_EQUALITY_CONST( x.size(), 0 );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY( x.dx(), null_ptr );
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, DMFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DMFad<value_type> ad_type;
  success = true;

  // Setup memory pool
  Sacado::Fad::MemPoolManager<value_type> poolManager(100);
  Sacado::Fad::MemPool *pool = poolManager.getMemoryPool(global_fad_size);
  ad_type::setDefaultPool(pool);

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  value_type *null_ptr = nullptr;
  TEST_EQUALITY_CONST( x.size(), 0 );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY_CONST( x.dx(), null_ptr );
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, DVFad_Owned )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DVFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = value_type(2.0+i);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, DVFad_Unowned )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DVFad<value_type> ad_type;
  success = true;

  // Initialize AD type
  value_type x_val = value_type(1.5);
  std::vector<value_type> x_dx(global_fad_size);
  for (int i=0; i<global_fad_size; ++i)
    x_dx[i] = value_type(2.0+i);
  ad_type x(global_fad_size, &x_val, x_dx.data(), 1, false);

  // Move x into y
  ad_type y(0.5);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_INEQUALITY( x.dx(), y.dx() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

TEUCHOS_UNIT_TEST( MoveAssignmentTests, ViewFad )
{
  typedef double value_type;
  typedef Sacado::Fad::Exp::DFad<value_type> dfad_type;
  typedef Sacado::Fad::Exp::ViewFad<value_type,0,1,dfad_type> ad_type;
  success = true;

  // Initialize AD type
  value_type x_val = value_type(1.5);
  std::vector<value_type> x_dx(global_fad_size);
  for (int i=0; i<global_fad_size; ++i)
    x_dx[i] = value_type(2.0+i);
  ad_type x(x_dx.data(), &x_val, global_fad_size, 1);

  // Move x into y
  value_type y_val = value_type(0.5);
  std::vector<value_type> y_dx(global_fad_size);
  for (int i=0; i<global_fad_size; ++i)
    y_dx[i] = value_type(20.0+i);
  ad_type y(y_dx.data(), &y_val, global_fad_size, 1);
  y = std::move(x);

  // Check y is correct
  TEST_EQUALITY_CONST( y.size(), global_fad_size );
  TEST_EQUALITY_CONST( y.val(), value_type(1.5) );
  TEST_EQUALITY( y.dx(), y_dx.data() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.dx(i), value_type(2.0+i) );
  }

  // Check x is correct
  TEST_EQUALITY_CONST( x.size(), global_fad_size );
  TEST_EQUALITY_CONST( x.val(), value_type(1.5) );
  TEST_EQUALITY( x.dx(), x_dx.data() );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( x.dx(i), value_type(2.0+i) );
  }
}

template <>
Sacado::Fad::MemPool* Sacado::Fad::Exp::MemPoolStorage<double>::defaultPool_ = nullptr;

#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
