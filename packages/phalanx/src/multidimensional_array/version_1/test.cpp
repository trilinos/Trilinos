
#include <vector>
#include <iostream>

#include <Dimension.hpp>
#include <Array.hpp>
#include <ArrayVector.hpp>
#include <ArrayScratch.hpp>

using namespace phdmesh ;

struct A : public DimTag {
  const char * name() const ;
  static const A & descriptor();
};

struct B : public DimTag {
  const char * name() const ;
  static const B & descriptor();
};

struct C : public DimTag {
  const char * name() const ;
  static const C & descriptor();
};

struct D : public DimTag {
  const char * name() const ;
  static const D & descriptor();
};

//----------------------------------------------------------------------

void myfunc_natural( const DimNatural<C,B,A> & d );
void myfunc_fortran( const DimFortran<A,B,C> & d );

// Enforce the convention:
// template< typename Type >
// void myfunc( const DimNatural<C,B,A> & d , const ArrayRCP<Type> & );

// Assume people follow the commented convention:
// template< typename Type >
// void myfunc( const ArrayRCP<Type> & );

int main()
{
  //----------------------------------
  // Rank 1 fortran and natural arrays,
  // automatic translation:

  DimFortran<A> d1_for( 10 );
  DimNatural<A> d1_nat( d1_for );

  //----------------------------------
  // Rank 2 fortran and natural arrays,
  // automatic translation:

  DimFortran<A,B> d2_for( 10 , 20 );
  DimNatural<B,A> d2_nat( d2_for );

  // DimNatural<A,B> d2_nat_compile_error( d2_for );

  //----------------------------------
  // Rank 3 fortran and natural arrays,
  // automatic translation:

  DimFortran<A,B,C> d3_for( 10 , 20 , 30 );
  DimFortran<C,A,B> d3_for_other( 10 , 20 , 30 );
  DimNatural<C,B,A> d3_nat( d3_for );
  DimNatural<B,A,C> d3_nat_other( d3_for_other );

  std::cout << d3_nat << std::endl ;
  std::cout << d3_for << std::endl ;

  // Calling functions with specific dimension conventions
  // to test the automatic translation between natural and fortran
  // conventions.

  myfunc_fortran( d3_nat );
  myfunc_fortran( d3_for );

  myfunc_natural( d3_nat );
  myfunc_natural( d3_for );

  // myfunc_fortran( d3_for_other ); // Compile error for incompatibility
  // myfunc_fortran( d3_nat_other ); // Compile error for incompatibility
  // myfunc_natural( d3_for_other ); // Compile error for incompatibility
  // myfunc_natural( d3_nat_other ); // Compile error for incompatibility

  //----------------------------------
  // Array that is a view to some storage:
  // There is no bounds-check possible with this
  // trust-based interface.

  double a[ 2000 ];

  Array< DimFortran<A,B> , double * > a2_for( 50 , 40 , a );
  Array< DimNatural<B,A> , double * > a2_nat( a2_for );

  std::cout << a2_for.dimension() << std::endl ;
  std::cout << a2_nat.dimension() << std::endl ;

  a2_for(1,2) = 5 ;

  // Test truncation:

  Array< DimNatural<A> , double * > a1_nat = a2_nat[5] ;
  Array< DimFortran<A> , double * > a1_for = a2_for[5] ;

  if ( a1_nat.data() != a1_for.data() ||
       a1_nat.data() != & a2_nat(5,0) ) {
    std::cout << "  Pointer-specialized Truncation failed" << std::endl ;
  }
  else {
    std::cout << "  Pointer-specialized Truncation succeded" << std::endl ;
  }

  //----------------------------------
  // Array that uses a vector for storage:

  Array< DimFortran<A,B> , std::vector<double> > vec2_for( 50 , 40 );

  std::cout << vec2_for.dimension() << std::endl ;

  vec2_for(1,2) = 5 ;

  // Test truncation:

  a1_nat = vec2_for[6] ;
  a1_for = vec2_for[6] ;

  if ( a1_nat.data() != a1_for.data() ||
       a1_nat.data() != & vec2_for(0,6) ) {
    std::cout << "  Vector-specialized Truncation failed" << std::endl ;
  }
  else {
    std::cout << "  Vector-specialized Truncation succeded" << std::endl ;
  }

  //----------------------------------
  // Array that is local, not-allocated, scratch space.
  // The dimensions are interpreted according to the
  // convention given by the second argument.

  Array< DimFortran<A,B> , double[50][40] > scr_for ;
  // Array< DimNatural<A,B> , double * > ptr_compile_error( scr_for );
  Array< DimNatural<B,A> , double * > ptr( scr_for );

  // Array< DimNatural<B,A> , ArrayRCP<double> > ptr( scr_for );

  std::cout << scr_for.dimension() << std::endl ;
  std::cout << ptr.dimension() << std::endl ;

  scr_for(1,2) = 5 ;

  // Test truncation:

  a1_nat = scr_for[6] ;
  a1_for = scr_for[6] ;

  if ( a1_nat.data() != a1_for.data() ||
       a1_nat.data() != & scr_for(0,6) ) {
    std::cout << "  Scratch-specialized Truncation failed" << std::endl ;
  }
  else {
    std::cout << "  Scratch-specialized Truncation succeded" << std::endl ;
  }
  //----------------------------------

  return 0 ;
}

void myfunc_natural( const DimNatural<C,B,A> & d )
{
  DimFortran<A,B,C> myd( d );
  std::cout << "myfunc_natural( " << d << " )" << std::endl ;
}

void myfunc_fortran( const DimFortran<A,B,C> & d )
{
  DimNatural<C,B,A> myd( d );
  std::cout << "myfunc_fortran( " << d << " )" << std::endl ;
}

//----------------------------------------------------------------------

const char * A::name() const { static const char n[] = "A" ; return n ; }
const char * B::name() const { static const char n[] = "B" ; return n ; }
const char * C::name() const { static const char n[] = "C" ; return n ; }
const char * D::name() const { static const char n[] = "D" ; return n ; }

const A & A::descriptor() { static const A myself ; return myself ; }
const B & B::descriptor() { static const B myself ; return myself ; }
const C & C::descriptor() { static const C myself ; return myself ; }
const D & D::descriptor() { static const D myself ; return myself ; }

