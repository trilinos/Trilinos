// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#define SHARDS_ARRAY_BOUNDS_CHECKING

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cstddef>
#include <Shards_Array.hpp>
#include <Shards_ArrayVector.hpp>

template< class T1 , class T2 > struct AssertSameType ;
template< class T > struct AssertSameType<T,T> {};

namespace {

//----------------------------------------------------------------------

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( TagA )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( TagA )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( TagB )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( TagB )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( TagC )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( TagC )

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION( TagD )
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( TagD )

//----------------------------------------------------------------------

using namespace shards ;

void myfortranfunc( const Array<double,FortranOrder> & xf )
{
  std::cout << "myfortranfunc( Array<double,FortranOrder" ;

  if ( xf.rank() && NULL != xf.tag(0) ) {
    for ( int i = 0 ; i < xf.rank() ; ++i ) {
      std::cout << "," << xf.tag(i)->name();
    }
  }

  std::cout << ">( " ;
  std::cout << (void*) xf.contiguous_data();
  for ( int i = 0 ; i < xf.rank() ; ++i ) {
    std::cout << " , " << xf.dimension(i);
  }
  std::cout << " ) )" << std::endl ;
}

void mynaturalfunc( const Array<double,NaturalOrder> & xf )
{
  std::cout << "mynaturalfunc( Array<double,NaturalOrder" ;

  if ( xf.rank() && NULL != xf.tag(0) ) {
    for ( int i = 0 ; i < xf.rank() ; ++i ) {
      std::cout << "," << xf.tag(i)->name();
    }
  }

  std::cout << ">( " ;
  std::cout << (void*) xf.contiguous_data();
  for ( int i = 0 ; i < xf.rank() ; ++i ) {
    std::cout << " , " << xf.dimension(i);
  }
  std::cout << " ) )" << std::endl ;
}

//----------------------------------------------------------------------

void myfortranA( const Array<double,FortranOrder,TagA> ) {}

void myfortranAB( const Array<double,FortranOrder,TagA,TagB> ) {}

void myfortranABC( const Array<double,FortranOrder,TagA,TagB,TagC> ) {}

void myfortranABCD( const Array<double,FortranOrder,TagA,TagB,TagC,TagD> ) {}

//----------------------------------------------------------------------

typedef Array<double,FortranOrder,TagA> AF1 ;
typedef Array<double,FortranOrder,TagA,TagB> AF2 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC> AF3 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC,TagD> AF4 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA> AF5 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB> AF6 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC> AF7 ;
typedef Array<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC,TagD> AF8 ;

typedef Array<double,NaturalOrder,TagA> AN1 ;
typedef Array<double,NaturalOrder,TagA,TagB> AN2 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC> AN3 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC,TagD> AN4 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA> AN5 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB> AN6 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC> AN7 ;
typedef Array<double,NaturalOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC,TagD> AN8 ;

void local_test_array()
{
  AssertSameType< AF2::TruncateType , AF1 >();
  AssertSameType< AF3::TruncateType , AF2 >();
  AssertSameType< AF4::TruncateType , AF3 >();
  AssertSameType< AF5::TruncateType , AF4 >();
  AssertSameType< AF6::TruncateType , AF5 >();
  AssertSameType< AF7::TruncateType , AF6 >();
  AssertSameType< AF8::TruncateType , AF7 >();

  AssertSameType< ArrayAppend< AF1 , TagB >::type , AF2 >();
  AssertSameType< ArrayAppend< AF2 , TagC >::type , AF3 >();
  AssertSameType< ArrayAppend< AF3 , TagD >::type , AF4 >();
  AssertSameType< ArrayAppend< AF4 , TagA >::type , AF5 >();
  AssertSameType< ArrayAppend< AF5 , TagB >::type , AF6 >();
  AssertSameType< ArrayAppend< AF6 , TagC >::type , AF7 >();
  AssertSameType< ArrayAppend< AF7 , TagD >::type , AF8 >();

  AssertSameType< AF8 , AF8::ReverseType::ReverseType >();
  AssertSameType< AF7 , AF7::ReverseType::ReverseType >();
  AssertSameType< AF6 , AF6::ReverseType::ReverseType >();
  AssertSameType< AF5 , AF5::ReverseType::ReverseType >();
  AssertSameType< AF4 , AF4::ReverseType::ReverseType >();
  AssertSameType< AF3 , AF3::ReverseType::ReverseType >();
  AssertSameType< AF2 , AF2::ReverseType::ReverseType >();
  AssertSameType< AF1 , AF1::ReverseType::ReverseType >();

  // AssertSameType< AF8 , AN8 >(); // Correctly fails to compile

  double storage[100000];

  AF1 af1( storage , 2 );
  AF2 af2( storage , 2 , 3 );
  AF3 af3( storage , 2 , 3 , 4 );
  AF4 af4( storage , 2 , 3 , 4 , 5 );
  AF5 af5( storage , 2 , 3 , 4 , 5 , 6 );
  AF6 af6( storage , 2 , 3 , 4 , 5 , 6 , 7 );
  AF7 af7( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 );
  AF8 af8( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 );

  AN1 an1( storage , 2 );
  AN2 an2( storage , 2 , 3 );
  AN3 an3( storage , 2 , 3 , 4 );
  AN4 an4( storage , 2 , 3 , 4 , 5 );
  AN5 an5( storage , 2 , 3 , 4 , 5 , 6 );
  AN6 an6( storage , 2 , 3 , 4 , 5 , 6 , 7 );
  AN7 an7( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 );
  AN8 an8( storage , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 );

  //--------------------------------

  std::vector<int> dim ;
  af6.dimensions( dim );

  if ( dim.size() != 6u || dim[0] != 2 ||
                           dim[1] != 3 ||
                           dim[2] != 4 ||
                           dim[3] != 5 ||
                           dim[4] != 6 ||
                           dim[5] != 7 ) {
    std::ostringstream msg ;
    msg << "Array<Fortran> dimension(vector) test failed" << std::endl
        << "  dim.size() = " << dim.size() << " vs. 6" << std::endl
        << "  dim[0]     = " << dim[0] << " vs. 2" << std::endl
        << "  dim[1]     = " << dim[1] << " vs. 3" << std::endl
        << "  dim[2]     = " << dim[2] << " vs. 4" << std::endl
        << "  dim[3]     = " << dim[3] << " vs. 5" << std::endl
        << "  dim[4]     = " << dim[4] << " vs. 6" << std::endl
        << "  dim[5]     = " << dim[5] << " vs. 7" << std::endl ;
    throw std::runtime_error( msg.str() );
  }

  an6.dimensions( dim );
  if ( dim.size() != 6u || dim[0] != 2 ||
                           dim[1] != 3 ||
                           dim[2] != 4 ||
                           dim[3] != 5 ||
                           dim[4] != 6 ||
                           dim[5] != 7 ) {
    std::ostringstream msg ;
    msg << "Array<Natural> dimension(vector) test failed" << std::endl
        << "  dim.size() = " << dim.size() << " vs. 6" << std::endl
        << "  dim[0]     = " << dim[0] << " vs. 2" << std::endl
        << "  dim[1]     = " << dim[1] << " vs. 3" << std::endl
        << "  dim[2]     = " << dim[2] << " vs. 4" << std::endl
        << "  dim[3]     = " << dim[3] << " vs. 5" << std::endl
        << "  dim[4]     = " << dim[4] << " vs. 6" << std::endl
        << "  dim[5]     = " << dim[5] << " vs. 7" << std::endl ;
    throw std::runtime_error( msg.str() );
  }

  //--------------------------------

  ArrayAppend< Array<double,FortranOrder,TagA,TagB,TagC> , TagD >::type
    aaf4( af4 );

  Array<double,FortranOrder,TagA,TagB,TagC,TagD> atf4( af5.truncate(0) );

  //--------------------------------

  std::cout << std::endl << "FORTRAN ARRAYS:" << std::endl ;

  myfortranfunc( af1 );
  myfortranfunc( af2 );
  myfortranfunc( af3 );
  myfortranfunc( af4 );
  myfortranfunc( af5 );
  myfortranfunc( af6 );
  myfortranfunc( af7 );
  myfortranfunc( af8 );

  mynaturalfunc( af1 );
  mynaturalfunc( af2 );
  mynaturalfunc( af3 );
  mynaturalfunc( af4 );
  mynaturalfunc( af5 );
  mynaturalfunc( af6 );
  mynaturalfunc( af7 );
  mynaturalfunc( af8 );

  myfortranfunc( af8.truncate(0) );

  std::cout << std::endl << "NATURAL ARRAYS:" << std::endl ;

  mynaturalfunc( an1 );
  mynaturalfunc( an2 );
  mynaturalfunc( an3 );
  mynaturalfunc( an4 );
  mynaturalfunc( an5 );
  mynaturalfunc( an6 );
  mynaturalfunc( an7 );
  mynaturalfunc( an8 );

  myfortranfunc( an1 );
  myfortranfunc( an2 );
  myfortranfunc( an3 );
  myfortranfunc( an4 );
  myfortranfunc( an5 );
  myfortranfunc( an6 );
  myfortranfunc( an7 );
  myfortranfunc( an8 );

  mynaturalfunc( an8.truncate(0) );

  //------------------------------

  myfortranA( af1 );
  myfortranA( an1 ); // Implicit conversion-construction is good

  // myfortranA( af2 ); // Compile error catches correctly
  // myfortranA( af3 ); // Compile error catches correctly
  // myfortranA( af4 ); // Compile error catches correctly

  myfortranAB( af2 );

  // myfortranAB( an2 ); // Compile error catches correctly
  // myfortranAB( af3 ); // Compile error catches correctly
  // myfortranAB( af4 ); // Compile error catches correctly

  myfortranABC( af3 );
  myfortranABCD( af4 );

  //------------------------------

  {
    {
      bool caught_it = false ;
      try { af8( af8.dimension<0>() ,0,0,0,0,0,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0, af8.dimension<1>() ,0,0,0,0,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0, af8.dimension<2>() ,0,0,0,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0,0, af8.dimension<3>() ,0,0,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0,0,0, af8.dimension<4>() ,0,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0,0,0,0, af8.dimension<5>() ,0,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0,0,0,0,0, af8.dimension<6>() ,0); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
    {
      bool caught_it = false ;
      try { af8(0,0,0,0,0,0,0, af8.dimension<7>() ); }
      catch( ... ) { caught_it = true ; }
      if ( ! caught_it ) {
        throw std::runtime_error(
          std::string("Array index bounds check failed") );
      }
    }
  }

  //------------------------------
  // Test flexibility of indexing type.

  {
    int i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }
  {
    long int i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }
  {
    unsigned i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }
  {
    long unsigned i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }
  {
    std::size_t i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }
  {
    std::ptrdiff_t i = 1 ;
    af8( i , i , i , i , i , i , i , i );
  }

  //------------------------------

  {
  AF8::ReverseType rf8( af8 );

  int count = 0 ;

  for ( int i7 = 0 ; i7 < af8.dimension<7>() ; ++i7 ) {
    for ( int i6 = 0 ; i6 < af8.dimension<6>() ; ++i6 ) {
      for ( int i5 = 0 ; i5 < af8.dimension<5>() ; ++i5 ) {
        for ( int i4 = 0 ; i4 < af8.dimension<4>() ; ++i4 ) {
          for ( int i3 = 0 ; i3 < af8.dimension<3>() ; ++i3 ) {
            for ( int i2 = 0 ; i2 < af8.dimension<2>() ; ++i2 ) {
              for ( int i1 = 0 ; i1 < af8.dimension<1>() ; ++i1 ) {
                for ( int i0 = 0 ; i0 < af8.dimension<0>() ; ++i0 ) {
                  ++count ;
                  if ( & af8(i0,i1,i2,i3,i4,i5,i6,i7) !=
                       & rf8(i7,i6,i5,i4,i3,i2,i1,i0) ) {
                    throw std::runtime_error(std::string("Failed index check"));
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if ( af8.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF7::ReverseType rf7( af7 );

  int count = 0 ;

  for ( int i6 = 0 ; i6 < af7.dimension<6>() ; ++i6 ) {
    for ( int i5 = 0 ; i5 < af7.dimension<5>() ; ++i5 ) {
      for ( int i4 = 0 ; i4 < af7.dimension<4>() ; ++i4 ) {
        for ( int i3 = 0 ; i3 < af7.dimension<3>() ; ++i3 ) {
          for ( int i2 = 0 ; i2 < af7.dimension<2>() ; ++i2 ) {
            for ( int i1 = 0 ; i1 < af7.dimension<1>() ; ++i1 ) {
              for ( int i0 = 0 ; i0 < af7.dimension<0>() ; ++i0 ) {
                ++count ;
                if ( & af7(i0,i1,i2,i3,i4,i5,i6) !=
                     & rf7(i6,i5,i4,i3,i2,i1,i0) ) {
                  throw std::runtime_error(std::string("Failed index check"));
                }
              }
            }
          }
        }
      }
    }
  }

  if ( af7.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF6::ReverseType rf6( af6 );

  int count = 0 ;

  for ( int i5 = 0 ; i5 < af6.dimension<5>() ; ++i5 ) {
    for ( int i4 = 0 ; i4 < af6.dimension<4>() ; ++i4 ) {
      for ( int i3 = 0 ; i3 < af6.dimension<3>() ; ++i3 ) {
        for ( int i2 = 0 ; i2 < af6.dimension<2>() ; ++i2 ) {
          for ( int i1 = 0 ; i1 < af6.dimension<1>() ; ++i1 ) {
            for ( int i0 = 0 ; i0 < af6.dimension<0>() ; ++i0 ) {
              ++count ;
              if ( & af6(i0,i1,i2,i3,i4,i5) !=
                   & rf6(i5,i4,i3,i2,i1,i0) ) {
                throw std::runtime_error(std::string("Failed index check"));
              }
            }
          }
        }
      }
    }
  }

  if ( af6.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF5::ReverseType rf5( af5 );

  int count = 0 ;

  for ( int i4 = 0 ; i4 < af5.dimension<4>() ; ++i4 ) {
    for ( int i3 = 0 ; i3 < af5.dimension<3>() ; ++i3 ) {
      for ( int i2 = 0 ; i2 < af5.dimension<2>() ; ++i2 ) {
        for ( int i1 = 0 ; i1 < af5.dimension<1>() ; ++i1 ) {
          for ( int i0 = 0 ; i0 < af5.dimension<0>() ; ++i0 ) {
            ++count ;
            if ( & af5(i0,i1,i2,i3,i4) !=
                 & rf5(i4,i3,i2,i1,i0) ) {
              throw std::runtime_error(std::string("Failed index check"));
            }
          }
        }
      }
    }
  }

  if ( af5.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF4::ReverseType rf4( af4 );

  int count = 0 ;

  for ( int i3 = 0 ; i3 < af4.dimension<3>() ; ++i3 ) {
    for ( int i2 = 0 ; i2 < af4.dimension<2>() ; ++i2 ) {
      for ( int i1 = 0 ; i1 < af4.dimension<1>() ; ++i1 ) {
        for ( int i0 = 0 ; i0 < af4.dimension<0>() ; ++i0 ) {
          ++count ;
          if ( & af4(i0,i1,i2,i3) !=
               & rf4(i3,i2,i1,i0) ) {
            throw std::runtime_error(std::string("Failed index check"));
          }
        }
      }
    }
  }

  if ( af4.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF3::ReverseType rf3( af3 );

  int count = 0 ;

  for ( int i2 = 0 ; i2 < af3.dimension<2>() ; ++i2 ) {
    for ( int i1 = 0 ; i1 < af3.dimension<1>() ; ++i1 ) {
      for ( int i0 = 0 ; i0 < af3.dimension<0>() ; ++i0 ) {
        ++count ;
        if ( & af3(i0,i1,i2) !=
             & rf3(i2,i1,i0) ) {
          throw std::runtime_error(std::string("Failed index check"));
        }
      }
    }
  }

  if ( af3.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF2::ReverseType rf2( af2 );

  int count = 0 ;

  for ( int i1 = 0 ; i1 < af2.dimension<1>() ; ++i1 ) {
    for ( int i0 = 0 ; i0 < af2.dimension<0>() ; ++i0 ) {
      ++count ;
      if ( & af2(i0,i1) !=
           & rf2(i1,i0) ) {
        throw std::runtime_error(std::string("Failed index check"));
      }
    }
  }

  if ( af2.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

  {
  AF1::ReverseType rf1( af1 );

  int count = 0 ;

  for ( int i0 = 0 ; i0 < af1.dimension<0>() ; ++i0 ) {
    ++count ;
    if ( & af1(i0) != & rf1(i0) ) {
      throw std::runtime_error(std::string("Failed index check"));
    }
  }

  if ( af1.size() != count ) {
    throw std::runtime_error( std::string("Failed loop check" ) );
  }
  }

  //------------------------------

}

//----------------------------------------------------------------------

void local_test_array_vector()
{
  typedef ArrayVector<double,FortranOrder,TagA> AVF1 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB> AVF2 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC> AVF3 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC,TagD> AVF4 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC,TagD,TagA> AVF5 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB> AVF6 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC> AVF7 ;
  typedef ArrayVector<double,FortranOrder,TagA,TagB,TagC,TagD,TagA,TagB,TagC,TagD> AVF8 ;

  AVF1 avf1( 2 );
  AVF2 avf2( 2 , 3 );
  AVF3 avf3( 2 , 3 , 4 );
  AVF4 avf4( 2 , 3 , 4 , 5 );
  AVF5 avf5( 2 , 3 , 4 , 5 , 6 );
  AVF6 avf6( 2 , 3 , 4 , 5 , 6 , 7 );
  AVF7 avf7( 2 , 3 , 4 , 5 , 6 , 7 , 8 );
  AVF8 avf8( 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 );
 
  std::cout << std::endl << "FORTRAN ARRAY-VECTORS:" << std::endl ;

  myfortranfunc( avf1 );
  myfortranfunc( avf2 );
  myfortranfunc( avf3 );
  myfortranfunc( avf4 );
  myfortranfunc( avf5 );
  myfortranfunc( avf6 );
  myfortranfunc( avf7 );
  myfortranfunc( avf8 );

  ArrayVector<double,FortranOrder> av_dynamic ;

  av_dynamic.resize<TagA>( 2 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB>( 2 , 3 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB,TagC>( 2 , 3 , 4 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB,TagC,TagD>( 2 , 3 , 4 , 5 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB,TagC,TagD,TagA>( 2 , 3 , 4 , 5 , 6 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB,TagC,TagD,TagA,TagB>( 2 , 3 , 4 , 5 , 6 , 7 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB,TagC,TagD,TagA,TagB>( 2 , 3 , 4 , 5 , 6 , 7 );
  myfortranfunc( av_dynamic );

  av_dynamic.resize<TagA,TagB>( 3 , 2 );
  myfortranfunc( av_dynamic );
}

//----------------------------------------------------------------------

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Point)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Point)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dim)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Dim)

void local_test_array_truncate()
{
  const int num_cells = 5;
  const int num_points = 4;
  const int num_dim = 3;
  const int offset_cell = 2 ;
  const int size = num_cells*num_points*num_dim;
  int stride[3] , dims[3] , indices[3] ;
  double* mem = new double[size];

  for (int i=0; i < size; ++i) {
    mem[i] = static_cast<double>(i);
  }

  Array<double,NaturalOrder,Cell,Point,Dim> a(mem,num_cells,num_points,num_dim);
  Array<double,NaturalOrder,     Point,Dim> b = a.truncate(offset_cell);

  dims[0] = num_cells ;
  dims[1] = num_points ;
  dims[2] = num_dim ;
  array_traits::stride_from_natural_dimensions( 3 , stride , dims );
  dims[0] = 0 ;
  dims[1] = 0 ;
  dims[2] = 0 ;
  array_traits::stride_to_natural_dimensions( 3 , stride , dims );

  if ( num_cells != dims[0] || num_points != dims[1] || num_dim != dims[2] ) {
    throw std::runtime_error( std::string("Array stride test failed") );
  }

  for (int c=0; c < num_cells; ++c) {
    for (int pt=0; pt < num_points; ++pt) {
      for (int dim=0; dim < num_dim; ++ dim) {
        const int offset = dim + num_dim * ( pt + num_points * c );
        const int value  = static_cast<int>( a(c,pt,dim) );

        array_traits::stride_to_natural_indices( 3, stride, offset, indices );

        if ( c != indices[0] || pt != indices[1] || dim != indices[2] ) {
          throw std::runtime_error( std::string("Array indices test failed") );
        }

        if ( a.contiguous_data() + offset != & a(c,pt,dim) ||
             offset != value ) {
          std::ostringstream msg ;
          msg << "Array offset test failed: required = " << offset ;
          msg << " != " << value << " = value" ;
          throw std::runtime_error( msg.str() );
        }
      }
    }
  }
   
  const ptrdiff_t truncate_offset = b.contiguous_data() - a.contiguous_data();
  const ptrdiff_t required_offset = num_dim * num_points * offset_cell ;

  if ( truncate_offset != required_offset ) {
    std::ostringstream msg ;
    msg << "Array Truncate test failed: required = " << required_offset ;
    msg << " != " << truncate_offset << " = truncate_offset" ;
    throw std::runtime_error( msg.str() );
  }

  delete [] mem;
}

//----------------------------------------------------------------------

}

void test_shards_array()
{
  static const char method[] = "test_shards_array" ;

  try {
    local_test_array();
    local_test_array_truncate();
    local_test_array_vector();
    std::cout << method << "\n" << "End Result: TEST PASSED" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << method << "\n" << "End Result: TEST FAILED: " << x.what() << std::endl ;
    throw x ;
  }
  catch( ... ) {
    std::cout << method << "\n" << "End Result: TEST FAILED: <unknown>" << std::endl ;
    throw ;
  }
}

