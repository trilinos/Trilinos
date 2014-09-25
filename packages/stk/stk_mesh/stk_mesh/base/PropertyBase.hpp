// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_mesh_base_PropertyBase_hpp
#define stk_mesh_base_PropertyBase_hpp

#include <iosfwd>
#include <string>
#include <map>
#include <vector>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace mesh {

/** \ingroup stk_mesh_module
 *  \brief  Property base class with an anonymous data type and
 *          anonymous multi-dimension.
 */
template<>
class Property< void > {
public:

  MetaData & mesh_meta_data() const { return m_meta_data ; }

  /** \brief  Internally generated ordinal of this property that is unique
   *          within the owning \ref stk::mesh::MetaData "meta data manager".
   */
  unsigned mesh_meta_data_ordinal() const { return m_meta_data_ordinal ; }

  /** \brief  Application-defined text name of this property */
  const std::string & name() const { return m_name ; }

  /** \brief  Query if the scalar type is DataType */
  template<class DataType> bool type_is() const
    { return m_type == typeid(DataType); }

  /** \brief  Query number of members, if an array */
  unsigned size() const { return m_size ; }

  /** \brief  Type-checked cast to property with members of the given type */
  template< typename DataType >
  Property< DataType > * property()
    {
      Property< DataType > * p = NULL ;
      if ( m_type == typeid(DataType) ) {
        p = static_cast< Property< DataType > * >( this );
      }
      return p ;
    }

  /** \brief  Type-checked cast to property with members of the given type */
  template< typename DataType >
  const Property< DataType > * property() const
    {
      const Property< DataType > * p = NULL ;
      if ( m_type == typeid(DataType) ) {
        p = static_cast< const Property< DataType > * >( this );
      }
      return p ;
    }

  //----------------------------------------

#ifndef DOXYGEN_COMPILE

protected:

  Property( MetaData & my_meta_data ,
            unsigned meta_data_ordinal ,
            const std::string & input_name ,
            const std::type_info & type ,
            unsigned n )
    : m_name( input_name ),
      m_meta_data( my_meta_data ),
      m_meta_data_ordinal( meta_data_ordinal ),
      m_type( type ), m_size( n ) {}

  virtual void add_property( unsigned ) = 0 ;

  virtual ~Property();

private:

  /** \brief  The \ref stk::mesh::MetaData "meta data manager"
   *          that owns this property.
   */
  MetaData & meta_data() const { return m_meta_data ; }

  const std::string      m_name ;       ///< Name of the property
  MetaData             & m_meta_data ;  ///< In which this property resides
  const unsigned         m_meta_data_ordinal ;
  const std::type_info & m_type ;       ///< Member data type
  const unsigned         m_size ;       ///< Number of items

  Property();
  Property( const Property & );
  Property & operator = ( const Property & );

  friend class MetaData ;

#endif /* DOXYGEN_COMPILE */
};

//----------------------------------------------------------------------

/** \ingroup stk_mesh_module
 *  \brief  Property with defined data type and multi-dimensions (if any)
 */

template< typename DataType >
class Property : public PropertyBase {
#ifndef DOXYGEN_COMPILE
private:
  friend class MetaData ;

  typedef std::map< unsigned , DataType > map_scalar ;

  map_scalar m_data_scalar ;

protected:

  Property( MetaData & my_meta_data, unsigned meta_data_ordinal ,
            const std::string & input_name, unsigned input_size = 1 )
    : PropertyBase( my_meta_data, meta_data_ordinal ,
                    input_name, typeid(DataType), input_size ) {}

  virtual void add_property( unsigned key ) { m_data_scalar[ key ]; }

  virtual ~Property() {}

#endif /* DOXYGEN_COMPILE */
public:

  /** \brief  Type of the properties' members */
  typedef DataType data_type ;

  /** \brief  Access the properties' members */
  virtual data_type * data( unsigned key )
  {
    const typename map_scalar::iterator i = m_data_scalar.find( key );
    return i != m_data_scalar.end() ? & (*i).second : static_cast<data_type*>(NULL) ;
  }

  /** \brief  Access the properties' members */
  virtual const data_type * data( unsigned key ) const
  {
    const typename map_scalar::const_iterator i = m_data_scalar.find( key );
    return i != m_data_scalar.end() ? & (*i).second : static_cast<data_type*>(NULL) ;
  }

};

#ifndef DOXYGEN_COMPILE

template< typename DataType >
class Property< std::vector< DataType > > : public Property<DataType> {
private:
  friend class MetaData ;

  typedef std::map< unsigned , std::vector< DataType > > map_array ;

  map_array m_data_array ;

  void add_property( unsigned key )
    { m_data_array[ key ].resize( Property<void>::size() ); }

  ~Property() {}
  Property();
  Property( const Property & );
  Property & operator = ( const Property & );

public:

  Property( MetaData & my_meta_data ,
            unsigned   meta_data_ordinal ,
            const std::string & name ,
            unsigned size )
    : Property<DataType>( my_meta_data, meta_data_ordinal, name, size ) {}

  typedef DataType data_type ;

  data_type * data( unsigned key )
  {
    const typename map_array::iterator i = m_data_array.find( key );
    return i != m_data_array.end() ? & (*i).second[0] : static_cast<data_type*>(NULL) ;
  }

  const data_type * data( unsigned key ) const
  {
    const typename map_array::const_iterator i = m_data_array.find( key );
    return i != m_data_array.end() ? & (*i).second[0] : static_cast<data_type*>(NULL) ;
  }
};

#endif /* DOXYGEN_COMPILE */


} // namespace mesh
} // namespace stk

#endif // stk_mesh_base_PropertyBase_hpp

