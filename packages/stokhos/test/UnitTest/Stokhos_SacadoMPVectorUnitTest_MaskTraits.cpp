// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_MP_Vector_MaskTraits.hpp"

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Create_8)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    std::cout << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << m1 << std::endl;
    TEST_EQUALITY( m1.getSize(), ensemble_size );
    TEST_EQUALITY( m1[0], true );
    TEST_EQUALITY( m1[1], false );
    TEST_EQUALITY( m1[2], true );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( m1[i], false );
    
    TEST_EQUALITY( (double) m1, 2./ensemble_size );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Create_16)
{
    constexpr int ensemble_size = 16;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    std::cout << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << m1 << std::endl;
    TEST_EQUALITY( m1.getSize(), ensemble_size );
    TEST_EQUALITY( m1[0], true );
    TEST_EQUALITY( m1[1], false );
    TEST_EQUALITY( m1[2], true );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( m1[i], false );
    
    TEST_EQUALITY( (double) m1, 2./ensemble_size );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Not_8)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = !m1;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    for (auto i=0; i<ensemble_size; ++i)
        TEST_EQUALITY( m2[i], !m1[i] );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Multiplication_8)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    scalar mul = m1*a;
    
    scalar mul2 = m1*b;
    scalar mul3 = b*m1;
    
    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;
    
    std::cout << mul2 << std::endl;
    std::cout << mul3 << std::endl;
    
    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 0.0 );
    TEST_EQUALITY( mul[2], 2.5 );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( mul[i], 0.0 );
    
    TEST_EQUALITY( mul2, mul3 );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Multiplication_16)
{
    constexpr int ensemble_size = 16;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    scalar mul = m1*a;
    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;
    
    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 0.0 );
    TEST_EQUALITY( mul[2], 2.5 );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( mul[i], 0. );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mul_Add_8)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    scalar mul = m1*a + !m1*b;
    scalar mul2 = a*m1 + !m1*b;
    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;
    std::cout << mul2 << std::endl;
    
    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 2.0 );
    TEST_EQUALITY( mul[2], 2.5 );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( mul[i], 2.0 );
    
    TEST_EQUALITY( mul, mul2 );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_DEFAULT)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    using namespace MaskLogic;
    
    scalar a = (scalar) 1.0;
    a[1] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>(scalar) 0.0;
    auto m3 = a> 0.0;
    auto m4 = 0.0<a;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    std::cout << m3<< std::endl;
    std::cout << m4<< std::endl;
    
    if (m1)
        {TEST_EQUALITY( true, false );}
    else
        {TEST_EQUALITY( true, true );}
    
    TEST_EQUALITY((bool) m1, false );
    TEST_EQUALITY((bool) !m1, true );
    
    if (m2)
        {TEST_EQUALITY( true, true );}
    else
        {TEST_EQUALITY( true, false );}
    
    TEST_EQUALITY((bool) m2, true );
    TEST_EQUALITY((bool) !m2, false );
    
    TEST_EQUALITY( m2, m3 );
    TEST_EQUALITY( m2, m4 );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_AND)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    using namespace MaskLogic;
    
    scalar a = (scalar) 1.0;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.0;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;

    
    TEST_EQUALITY( (AND) true, true );
    TEST_EQUALITY( (AND) false, false );
    TEST_EQUALITY( (AND) m1, false );
    TEST_EQUALITY( (AND) !m1, false );
    TEST_EQUALITY( (AND) m2, true );
    TEST_EQUALITY( (AND) !m2, false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_OR)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    using namespace MaskLogic;
    
    scalar a = (scalar) 1.0;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.0;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    
    
    TEST_EQUALITY( (OR) true, true );
    TEST_EQUALITY( (OR) false, false );
    TEST_EQUALITY( (OR) m1, true );
    TEST_EQUALITY( (OR) !m1, true );
    TEST_EQUALITY( (OR) m2, true );
    TEST_EQUALITY( (OR) !m2, false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_XOR)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    using namespace MaskLogic;
    
    scalar a = (scalar) 1.0;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.0;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    
    
    TEST_EQUALITY( (XOR) true, true );
    TEST_EQUALITY( (XOR) false, false );
    TEST_EQUALITY( (XOR) m1, true );
    TEST_EQUALITY( (XOR) !m1, false );
    TEST_EQUALITY( (XOR) m2, false );
    TEST_EQUALITY( (XOR) !m2, false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_compared_to_double)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.0;
    
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    
    TEST_EQUALITY((double) m1,1./ensemble_size);
    TEST_EQUALITY((double) m2,1.);

    TEST_EQUALITY(m1==1.,false);
    TEST_EQUALITY(m1!=1.,true);
    TEST_EQUALITY(m1==0.,false);
    TEST_EQUALITY(m1!=0.,true);
    
    TEST_EQUALITY(m1>=0.5,false);
    TEST_EQUALITY(m1<=0.5,true);
    TEST_EQUALITY(m1>0.5,false);
    TEST_EQUALITY(m1<0.5,true);
    
    TEST_EQUALITY(m2==1.,true);
    TEST_EQUALITY(m2!=1.,false);
    TEST_EQUALITY(m2==0.,false);
    TEST_EQUALITY(m2!=0.,true);
    
    TEST_EQUALITY(m2>=0.5,true);
    TEST_EQUALITY(m2<=0.5,false);
    TEST_EQUALITY(m2>0.5,true);
    TEST_EQUALITY(m2<0.5,false);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_AND_Mask)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.;
    
    auto m3 = m1 && m2;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    std::cout << m3 << std::endl;
    TEST_EQUALITY(m3,m1);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_OR_Mask)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.0;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.;
    
    auto m3 = m1 || m2;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    std::cout << m3 << std::endl;
    TEST_EQUALITY(m3,m2);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_ADD_Mask)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    
    auto m1 = a>b;
    auto m2 = a>0.;
    //std::cout << m1 << std::endl;
    //std::cout << m2 << std::endl;
    auto m3 = m1 + m2;

    std::cout << m3 << std::endl;
    TEST_EQUALITY(m3,m2);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_SUB_Mask)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    
    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.0;
    
    auto m1 = a>b;
    auto m2 = a>0.;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    auto m3 = (a>0.) - (a>b);
    std::cout << m3 << std::endl;
    TEST_EQUALITY(m3,!m1);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_signbit_v)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    typedef Mask<scalar> mask;
    
    scalar a = (scalar) 1.;
    a[2] = -2.5;
    
    auto m1 = signbit_v(a);
    mask m2;
    m2[2] = true;
    TEST_EQUALITY(m1,m2);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_copysign)
{
    constexpr int ensemble_size = 8;
    
    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;
    typedef Mask<scalar> mask;
    
    scalar a = (scalar) 1.;
    a[2] = -2.5;
    
    scalar b = (scalar) 2.;
    
    using std::copysign;
    
    std::cout << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    b = copysign(b,a);
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    TEST_EQUALITY(b[2],-2.);
}


int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Kokkos::initialize();
//  Kokkos::HostSpace::execution_space::initialize();
//  if (!Kokkos::DefaultExecutionSpace::is_initialized())
//    Kokkos::DefaultExecutionSpace::initialize();

  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);

  Kokkos::finalize();
//  Kokkos::HostSpace::execution_space::finalize();
//  if (Kokkos::DefaultExecutionSpace::is_initialized())
//    Kokkos::DefaultExecutionSpace::finalize();

  return res;
}
