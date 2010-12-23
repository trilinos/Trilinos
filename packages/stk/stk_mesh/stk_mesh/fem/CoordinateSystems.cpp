/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdlib.h>

#include <sstream>
#include <stdexcept>

#include <stk_util/util/string_case_compare.hpp>
#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

namespace stk {
namespace mesh {

//----------------------------------------------------------------------

namespace {

unsigned get_index( const char * const func ,
                    const unsigned number_names ,
                    const char * const * names ,
                    const unsigned size ,
                    const char * const select )
{
  unsigned index = size <= number_names ? 0 : size ;

  for ( ; index < size && not_equal_case(select,names[index]) ; ++index );

  ThrowErrorMsgIf( index == size,
                   func << ", size = " << size << " label = " << select );
  return index ;
}

const char * get_string( const char * const func ,
                         const unsigned number_names ,
                         const char * const * names ,
                         const unsigned size ,
                         const unsigned index )
{
  ThrowErrorMsgIf( size < number_names || size <= index,
                   func << ", size = " << size << " index = " << index );

  return names[index];
}

}

//----------------------------------------------------------------------

const Cartesian & Cartesian::tag()
{ static const Cartesian self ; return self ; }

const char * Cartesian::name() const
{ static const char n[] = "Cartesian" ; return n ; }

std::string Cartesian::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * label[] = { x , y , z };

  return std::string( get_string( Cartesian::tag().name() ,
                                  3 , label , size , index ) );
}

shards::ArrayDimTag::size_type Cartesian::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * label[] = { x , y , z };

  return get_index( Cartesian::tag().name() ,
                    3 , label , size , arg.c_str() );
}

//----------------------------------------------------------------------

const Cylindrical & Cylindrical::tag()
{ static const Cylindrical self ; return self ; }

const char * Cylindrical::name() const
{ static const char n[] = "Cylindrical" ; return n ; }

std::string Cylindrical::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * label[] = { r , a , z };

  return std::string( get_string( Cylindrical::tag().name() ,
                                  3 , label , size , index ) );
}

shards::ArrayDimTag::size_type Cylindrical::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * label[] = { r , a , z };

  return get_index( Cylindrical::tag().name() ,
                    3 , label , size , arg.c_str() );
}

//----------------------------------------------------------------------

const FullTensor & FullTensor::tag()
{ static const FullTensor self ; return self ; }

const char * FullTensor::name() const
{ static const char n[] = "FullTensor" ; return n ; }

std::string FullTensor::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  static const char xx[] = "xx" ;
  static const char yx[] = "yx" ;
  static const char zx[] = "zx" ;
  static const char xy[] = "xy" ;
  static const char yy[] = "yy" ;
  static const char zy[] = "zy" ;
  static const char xz[] = "xz" ;
  static const char yz[] = "yz" ;
  static const char zz[] = "zz" ;
  static const char * label[] = { xx , yx , zx , xy , yy , zy , xz , yz , zz };

  return std::string( get_string( FullTensor::tag().name() ,
                                  9 , label , size , index ) );
}

shards::ArrayDimTag::size_type FullTensor::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  static const char xx[] = "xx" ;
  static const char yx[] = "yx" ;
  static const char zx[] = "zx" ;
  static const char xy[] = "xy" ;
  static const char yy[] = "yy" ;
  static const char zy[] = "zy" ;
  static const char xz[] = "xz" ;
  static const char yz[] = "yz" ;
  static const char zz[] = "zz" ;
  static const char * label[] = { xx , yx , zx , xy , yy , zy , xz , yz , zz };

  return get_index( FullTensor::tag().name() ,
                    9 , label , size , arg.c_str() );
}

//----------------------------------------------------------------------

const SymmetricTensor & SymmetricTensor::tag()
{ static const SymmetricTensor self ; return self ; }

const char * SymmetricTensor::name() const
{ static const char n[] = "SymmetricTensor" ; return n ; }

std::string SymmetricTensor::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  static const char xx[] = "xx" ;
  static const char yx[] = "yx" ;
  static const char zx[] = "zx" ;
  static const char xy[] = "xy" ;
  static const char yy[] = "yy" ;
  static const char zy[] = "zy" ;
  static const char xz[] = "xz" ;
  static const char yz[] = "yz" ;
  static const char zz[] = "zz" ;
  static const char * label[] = { xx , yx , zx , xy , yy , zy , xz , yz , zz };

  return std::string( get_string( SymmetricTensor::tag().name() ,
                                  9 , label , size , index ) );
}

shards::ArrayDimTag::size_type SymmetricTensor::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  static const char xx[] = "xx" ;
  static const char yy[] = "yy" ;
  static const char zz[] = "zz" ;

  static const char xy[] = "xy" ;
  static const char yz[] = "yz" ;
  static const char xz[] = "xz" ;

  static const char * label[] = { xx , yy , zz , xy , yz , xz };

  return get_index( SymmetricTensor::tag().name() ,
                    6 , label , size , arg.c_str() );
}

//----------------------------------------------------------------------

}//namespace mesh
}//namespace stk

