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

const Cartesian2d & Cartesian2d::tag()
{ static const Cartesian2d self ; return self ; }

const char * Cartesian2d::name() const
{ static const char n[] = "Cartesian2d" ; return n ; }

namespace {
const char * const * Cartesian2d_label() {
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char * label[] = { x , y };
  return label;
}
}
std::string Cartesian2d::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( Cartesian2d::tag().name() ,
                                  2 , Cartesian2d_label() , size , index ) );
}

shards::ArrayDimTag::size_type Cartesian2d::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( Cartesian2d::tag().name() ,
                    2 , Cartesian2d_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const Cartesian3d & Cartesian3d::tag()
{ static const Cartesian3d self ; return self ; }

const char * Cartesian3d::name() const
{ static const char n[] = "Cartesian3d" ; return n ; }

namespace {
const char * const * Cartesian3d_label() {
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * label[] = { x , y , z };
  return label;
}
}
std::string Cartesian3d::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( Cartesian3d::tag().name() ,
                                  3 , Cartesian3d_label() , size , index ) );
}

shards::ArrayDimTag::size_type Cartesian::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( Cartesian3d::tag().name() ,
                    3 , Cartesian3d_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const Cylindrical & Cylindrical::tag()
{ static const Cylindrical self ; return self ; }

const char * Cylindrical::name() const
{ static const char n[] = "Cylindrical" ; return n ; }

namespace {
const char * const * Cylindrical_label() {
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * label[] = { r , a , z };
  return label;
}
}
std::string Cylindrical::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( Cylindrical::tag().name() ,
                                  3 , Cylindrical_label() , size , index ) );
}

shards::ArrayDimTag::size_type Cylindrical::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( Cylindrical::tag().name() ,
                    3 , Cylindrical_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const FullTensor & FullTensor::tag()
{ static const FullTensor self ; return self ; }

const char * FullTensor::name() const
{ static const char n[] = "FullTensor" ; return n ; }

namespace {
const char * const * FullTensor36_label() {
  static const char xx[] = "xx" ;
  static const char yx[] = "yx" ;
  static const char zx[] = "zx" ;
  static const char xy[] = "xy" ;
  static const char yy[] = "yy" ;
  static const char zy[] = "zy" ;
  static const char xz[] = "xz" ;
  static const char yz[] = "yz" ;
  static const char zz[] = "zz" ;
  static const char * label[] = { xx , yy , zz , 
                                  xy , yz , zx , 
                                  yx , zy , xz };
  return label;
}
}
std::string FullTensor36::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( FullTensor36::tag().name() ,
                                  9 , FullTensor36_label() , size , index ) );
}

shards::ArrayDimTag::size_type FullTensor::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( FullTensor36::tag().name() ,
                    9 , FullTensor36_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const FullTensor22 & FullTensor22::tag()
{ static const FullTensor22 self ; return self ; }

const char * FullTensor22::name() const
{ static const char n[] = "FullTensor22" ; return n ; }

namespace {
const char * const * FullTensor22_label() {
  static const char xx[] = "xx" ;
  static const char yy[] = "yy" ;
  static const char xy[] = "xy" ;
  static const char yx[] = "yx" ;
  static const char * label[] = { xx, yy, xy, yx };
  return label;
}
}
std::string FullTensor22::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( FullTensor22::tag().name() ,
                                  4 , FullTensor22_label() , size , index ) );
}

shards::ArrayDimTag::size_type FullTensor22::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( FullTensor22::tag().name() ,
                    4 , FullTensor22_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const SymmetricTensor33 & SymmetricTensor33::tag()
{ static const SymmetricTensor33 self ; return self ; }

const char * SymmetricTensor::name() const
{ static const char n[] = "SymmetricTensor" ; return n ; }

namespace {
const char * const * SymmetricTensor33_label() {
  static const char xx[] = "xx" ;
  static const char yy[] = "yy" ;
  static const char zz[] = "zz" ;
  static const char xy[] = "xy" ;
  static const char yz[] = "yz" ;
  static const char xz[] = "xz" ;
  static const char * label[] = { xx , yy , zz , xy , yz , xz };
  return label;
}
}
std::string SymmetricTensor33::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( SymmetricTensor33::tag().name() ,
                                  6 , SymmetricTensor33_label() , size , index ) );
}

shards::ArrayDimTag::size_type SymmetricTensor33::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( SymmetricTensor33::tag().name() ,
                    6 , SymmetricTensor33_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const SymmetricTensor31 & SymmetricTensor31::tag()
{ static const SymmetricTensor31 self ; return self ; }

const char * SymmetricTensor31::name() const
{ static const char n[] = "SymmetricTensor31" ; return n ; }

namespace {
const char * const * SymmetricTensor31_label() {
  static const char rr[] = "rr" ;
  static const char zz[] = "zz" ;
  static const char rz[] = "rz" ;
  static const char zr[] = "zr" ;
  static const char * label[] = { rr, zz, rz, zr };
  return label;
}
}
std::string SymmetricTensor31::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( SymmetricTensor31::tag().name() ,
                                  4 , SymmetricTensor31_label() , size , index ) );
}

shards::ArrayDimTag::size_type SymmetricTensor31::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( SymmetricTensor31::tag().name() ,
                    4 , SymmetricTensor31_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const SymmetricTensor21 & SymmetricTensor21::tag()
{ static const SymmetricTensor21 self ; return self ; }

const char * SymmetricTensor21::name() const
{ static const char n[] = "SymmetricTensor21" ; return n ; }

namespace {
const char * const * SymmetricTensor21_label() {
  static const char xx[] = "xx" ;
  static const char yy[] = "yy" ;
  static const char xy[] = "xy" ;
  static const char * label[] = { xx, yy, xy };
  return label;
}
}
std::string SymmetricTensor21::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( SymmetricTensor21::tag().name() ,
                                  3 , SymmetricTensor21_label() , size , index ) );
}

shards::ArrayDimTag::size_type SymmetricTensor21::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( SymmetricTensor21::tag().name() ,
                    3 , SymmetricTensor21_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const AsymmetricTensor03 & AsymmetricTensor03::tag()
{ static const AsymmetricTensor03 self ; return self ; }

const char * AsymmetricTensor03::name() const
{ static const char n[] = "AsymmetricTensor03" ; return n ; }

namespace {
const char * const * AsymmetricTensor03_label() {
  static const char yz[] = "yz" ;
  static const char xz[] = "xz" ;
  static const char xy[] = "xy" ;
  static const char * label[] = { xy, yz, xz };
  return label;
}
}
std::string AsymmetricTensor03::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( AsymmetricTensor03::tag().name() ,
                                  3 , AsymmetricTensor03_label() , size , index ) );
}

shards::ArrayDimTag::size_type AsymmetricTensor03::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( AsymmetricTensor03::tag().name() ,
                    3 , AsymmetricTensor03_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const Matrix22 & Matrix22::tag()
{ static const Matrix22 self ; return self ; }

const char * Matrix22::name() const
{ static const char n[] = "Matrix22" ; return n ; }

namespace {
const char * const * Matrix22_label() {
  static const char xx[] = "xx" ;
  static const char yx[] = "yx" ;
  static const char xy[] = "xy" ;
  static const char yy[] = "yy" ;
  static const char * label[] = { xx , yx, xy, yy };
  return label;
}
}
std::string Matrix22::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( Matrix22::tag().name() ,
                                  4 , Matrix22_label() , size , index ) );
}

shards::ArrayDimTag::size_type Matrix22::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( Matrix22::tag().name() ,
                    4 , Matrix22_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

const Matrix33 & Matrix33::tag()
{ static const Matrix33 self ; return self ; }

const char * Matrix33::name() const
{ static const char n[] = "Matrix33" ; return n ; }

namespace {
const char * const * Matrix33_label() {
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
  return label;
}
}
std::string Matrix33::to_string( shards::ArrayDimTag::size_type size , shards::ArrayDimTag::size_type index ) const
{
  return std::string( get_string( Matrix33::tag().name() ,
                                  9 , Matrix33_label() , size , index ) );
}

shards::ArrayDimTag::size_type Matrix33::to_index( shards::ArrayDimTag::size_type size , const std::string & arg ) const
{
  return get_index( Matrix33::tag().name() ,
                    9 , Matrix33_label() , size , arg.c_str() );
}

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

