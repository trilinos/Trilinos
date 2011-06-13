/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef STK_MESH_DATA_TRAITS_HPP
#define STK_MESH_DATA_TRAITS_HPP

#include <cstddef>
#include <iosfwd>
#include <typeinfo>
#include <vector>
#include <string>
#include <complex>

#include <stk_util/parallel/ParallelComm.hpp>

namespace stk {
namespace mesh {

class DataTraits ;

//----------------------------------------------------------------------
/** \brief  Query singleton for data traits of a given data type. */
template< typename T > const DataTraits & data_traits();
 
/** \brief  Query DataTraits for a given data value. */
template< typename T >
inline
const DataTraits & data_traits( const T & ) { return data_traits<T>(); }

//----------------------------------------------------------------------
// Data traits for the fundamental computational data types:

template<> const DataTraits & data_traits< void >();
template<> const DataTraits & data_traits< signed   char >();
template<> const DataTraits & data_traits< unsigned char >();
template<> const DataTraits & data_traits< signed   short >();
template<> const DataTraits & data_traits< unsigned short >();
template<> const DataTraits & data_traits< signed   int >();
template<> const DataTraits & data_traits< unsigned int >();
template<> const DataTraits & data_traits< signed   long >();
template<> const DataTraits & data_traits< unsigned long >();
template<> const DataTraits & data_traits< float >();
template<> const DataTraits & data_traits< double >();
template<> const DataTraits & data_traits< std::complex<float> >();
template<> const DataTraits & data_traits< std::complex<double> >();

template<> const DataTraits & data_traits< void * >();
template<> const DataTraits & data_traits< signed   char * >();
template<> const DataTraits & data_traits< unsigned char * >();
template<> const DataTraits & data_traits< signed   short * >();
template<> const DataTraits & data_traits< unsigned short * >();
template<> const DataTraits & data_traits< signed   int * >();
template<> const DataTraits & data_traits< unsigned int * >();
template<> const DataTraits & data_traits< signed   long * >();
template<> const DataTraits & data_traits< unsigned long * >();
template<> const DataTraits & data_traits< float * >();
template<> const DataTraits & data_traits< double * >();
template<> const DataTraits & data_traits< std::complex<float> * >();
template<> const DataTraits & data_traits< std::complex<double> * >();

//----------------------------------------------------------------------

class DataTraits {
public:
  //------------------------------
  // Standard properties:
  const std::type_info & type_info ;
  std::size_t            size_of ;
 
  //------------------------------
  // TR1 primary type categories:
  bool         is_void ;
  bool         is_integral ;
  bool         is_floating_point ;
  bool         is_array ;
  bool         is_pointer ;
  bool         is_enum ;
  bool         is_class ;

  // TR1 type properties:
  bool         is_pod ;
  bool         is_signed ;     // only if 'is_integral'
  bool         is_unsigned ;   // only if 'is_integral'
  std::size_t  alignment_of ;

  // For memory management contiguous arrays of data:
  // Array must start aligned with 'alignment_of' and
  // stride by 'stride_of'.
  std::size_t  stride_of ;

  // TR1 type manipulators:
  const DataTraits * remove_pointer ; // if 'is_pointer'
 
  //------------------------------
  /** \brief  Namespace-qualified text name as it appears in source code */
  std::string  name ;
 
  //------------------------------
  // Only If 'is_enum'
  struct EnumMember {
    std::string  name ;
    long         value ;
  };
  std::vector< EnumMember > enum_info ;
 
  //------------------------------
  // Only If 'is_class':
  struct ClassMember {
    std::string        name ;
    const DataTraits * traits ;
    std::size_t        offset ;
  };
  std::vector< ClassMember > class_info ;
 
  //------------------------------
  // Functions required for all field data:

  virtual void construct( void * , std::size_t ) const  = 0 ;
  virtual void destroy(   void * , std::size_t ) const  = 0 ;
  virtual void copy( void * , const void * , std::size_t ) const  = 0 ;
  virtual void pack(   CommBuffer & , const void * , std::size_t ) const  = 0 ;
  virtual void unpack( CommBuffer & ,       void * , std::size_t ) const  = 0 ;
  virtual void print( std::ostream & , const void * , std::size_t ) const  = 0 ;

  //------------------------------
  // Commutative and associative ops
  // required for is_integral and is_floating_point data.
  // In-place reduction: x[0..(n-1)] op= y[0..(n-1)]
  virtual void sum( void * x , const void * y , std::size_t n ) const = 0 ;
  virtual void max( void * x , const void * y , std::size_t n ) const = 0 ;
  virtual void min( void * x , const void * y , std::size_t n ) const = 0 ;

  // Commutative and associative ops
  // required for is_integral data.
  // In-place reduction: x[0..(n-1)] op= y[0..(n-1)]
  virtual void bit_and( void * x , const void * y, std::size_t n ) const = 0 ;
  virtual void bit_or(  void * x , const void * y, std::size_t n ) const = 0 ;
  virtual void bit_xor( void * x , const void * y, std::size_t n ) const = 0 ;

  //------------------------------

protected:

  //------------------------------
  /** \brief  CTOR for most types. */
  DataTraits( const std::type_info & arg_type ,
              const char * const     arg_name ,
              const std::size_t      arg_size ,
              const std::size_t      arg_align );
 
  /** \brief  CTOR for pointer type */
  DataTraits( const std::type_info & arg_type , const DataTraits & );

  virtual ~DataTraits();
private:
  DataTraits();
  DataTraits( const DataTraits & );
  DataTraits & operator = ( const DataTraits & );
};

} // namespace mesh
} // namespace stk

#endif

