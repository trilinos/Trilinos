// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Create_8)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    std::cout << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << m1 << std::endl;
    TEST_EQUALITY( m1.getSize(), ensemble_size );
    TEST_EQUALITY( m1.get(0), true );
    TEST_EQUALITY( m1.get(1), false );
    TEST_EQUALITY( m1.get(2), true );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( m1.get(i), false );

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
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    std::cout << std::endl;
    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << m1 << std::endl;
    TEST_EQUALITY( m1.getSize(), ensemble_size );
    TEST_EQUALITY( m1.get(0), true );
    TEST_EQUALITY( m1.get(1), false );
    TEST_EQUALITY( m1.get(2), true );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( m1.get(i), false );

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
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = !m1;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;
    for (auto i=0; i<ensemble_size; ++i)
        TEST_EQUALITY( m2.get(i), !m1.get(i) );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Multiplication_8)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    scalar mul = m1*a;

    scalar mul2 = m1*b;
    scalar mul3 = b*m1;

    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;

    std::cout << mul2 << std::endl;
    std::cout << mul3 << std::endl;

    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 0. );
    TEST_EQUALITY( mul[2], 2.5 );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( mul[i], 0. );

    TEST_EQUALITY( mul2, mul3 );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Multiplication_16)
{
    constexpr int ensemble_size = 16;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    scalar mul = m1*a;
    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;

    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 0. );
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

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    scalar mul = m1*a + !m1*b;
    scalar mul2 = a*m1 + !m1*b;
    std::cout << m1 << std::endl;
    std::cout << mul << std::endl;
    std::cout << mul2 << std::endl;

    TEST_EQUALITY( mul[0], 2.5 );
    TEST_EQUALITY( mul[1], 2. );
    TEST_EQUALITY( mul[2], 2.5 );
    for (auto i=3; i<ensemble_size; ++i)
        TEST_EQUALITY( mul[i], 2. );

    TEST_EQUALITY( mul, mul2 );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_DEFAULT)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    using namespace MaskLogic;

    scalar a = (scalar) 1.;
    a[1] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = a>(scalar) 0.;
    auto m3 = a> 0.;
    auto m4 = 0.<a;
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

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = a>0.;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;


    TEST_EQUALITY( AND(true), true );
    TEST_EQUALITY( AND(false), false );
    TEST_EQUALITY( AND(m1), false );
    TEST_EQUALITY( AND(!m1), false );
    TEST_EQUALITY( AND(m2), true );
    TEST_EQUALITY( AND(!m2), false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_OR)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    using namespace MaskLogic;

    scalar a = (scalar) 1.;
    a[0] = 2.5;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = a>0.;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;


    TEST_EQUALITY( OR(true), true );
    TEST_EQUALITY( OR(false), false );
    TEST_EQUALITY( OR(m1), true );
    TEST_EQUALITY( OR(!m1), true );
    TEST_EQUALITY( OR(m2), true );
    TEST_EQUALITY( OR(!m2), false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_XOR)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    using namespace MaskLogic;

    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = a>0.;
    std::cout << m1 << std::endl;
    std::cout << m2 << std::endl;


    TEST_EQUALITY( XOR(true), true );
    TEST_EQUALITY( XOR(false), false );
    TEST_EQUALITY( XOR(m1), true );
    TEST_EQUALITY( XOR(!m1), false );
    TEST_EQUALITY( XOR(m2), false );
    TEST_EQUALITY( XOR(!m2), false );
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_compared_to_double)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

    auto m1 = a>b;
    auto m2 = a>0.;

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

    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

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

    scalar a = (scalar) 1.;
    a[2] = 2.5;
    scalar b = (scalar) 2.;

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
    scalar b = (scalar) 2.;

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
    scalar b = (scalar) 2.;

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
    m2.set(2,true);
    TEST_EQUALITY(m1,m2);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_copysign)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

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

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_assign)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[2] = -2.5;

    mask_assign(a<=0.,a) = {0.,a};

    TEST_EQUALITY(a[1],1.);
    TEST_EQUALITY(a[2],0.);

    double b = 1.;

    mask_assign(b>0.5 && b<2.,b) = {2.*b,-1.};

    TEST_EQUALITY(b,2.);

    mask_assign(b>0.5 && b<2.,b) = {2.*b,-1.};
    TEST_EQUALITY(b,-1.);

}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_pointer_assign)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;

    scalar a = (scalar) 1.;
    a[2] = -2.5;
    scalar *p = &a;

    mask_assign(a<=0.,*p) = {0.,a};

    TEST_EQUALITY(a[1],1.);
    TEST_EQUALITY(a[2],0.);
}

TEUCHOS_UNIT_TEST( MP_Vector_MaskTraits, Mask_div)
{
    constexpr int ensemble_size = 8;

    typedef Kokkos::DefaultExecutionSpace execution_space;
    typedef Stokhos::StaticFixedStorage<int,double,ensemble_size,execution_space> storage_type;
    typedef Sacado::MP::Vector<storage_type> scalar;    

    scalar a2 = {0.,2.,2.,2.,2.,2.,2.,2.};
    std::cout << a2 << std::endl;

    scalar a = (scalar) 1.;
    a[2] = -2.5;
    auto m = (a>(scalar) 0.);
    std::cout << "m is computed" << std::endl;
    std::cout << m << std::endl;
    m = a>0.;
    std::cout << "m is computed" << std::endl;
    std::cout << m << std::endl;

    std::cout << a << std::endl;
    std::cout << m << std::endl;
    std::cout << (a>=(scalar) 0. )<< std::endl;
    std::cout << (a> 0. )<< std::endl;
    std::cout << (a>= 0.) << std::endl;
    std::cout << (0.<a )<< std::endl;
    std::cout << (0.<=a) << std::endl;

    mask_assign<scalar>(m,a) /= {a, 2.,-1.};
    TEST_EQUALITY(a[1],0.5);
    TEST_EQUALITY(a[2],-1.);

    /*
     This test is working only if c++ 14 is allowed due to the fact
     that copy-list-initialization in the constructor of tuple is not allowed before
     as it is an explicit one.

    mask_assign<scalar>(m,a) /= {(scalar) 4.,2.,-1.};
     */
    //std::tuple<scalar,scalar,scalar> ts {4.,2.,-1.};
    mask_assign<scalar>(m,a) /= {4.,2.,-1.};
    TEST_EQUALITY(a[1],2.);
    TEST_EQUALITY(a[2],-1.);


    double b = 1.;
    mask_assign(b>0.5,b) /= {b, 2.,-1.};
    TEST_EQUALITY(b,0.5);
    mask_assign(b>0.5,b) /= {b, 2.,-1.};
    TEST_EQUALITY(b,-1.);

}
