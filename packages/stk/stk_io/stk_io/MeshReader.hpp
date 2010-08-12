/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_io_MeshReader_hpp
#define stk_io_MeshReader_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>

class CellTopologyData;

namespace stk {
  namespace mesh {
    class MetaData;
    class BulkData;
    class Part;
  }
namespace io {

//----------------------------------------------------------------------

class MeshReader ;
class MeshReaderFilter ;

/** \brief  Allocate a mesh reader for the given file format
 *          and populate the meta data according to the user
 *          provided filter.
 */
MeshReader * new_reader( ParallelMachine ,
                         mesh::MetaData & ,
                         mesh::VectorField & coord_field ,
                         const MeshReaderFilter & filter ,
                         const std::string & file_format ,
                         const std::string & file_name );

//----------------------------------------------------------------------

class MeshReader {
public:

  /** \brief  Read the entities, relations, coordinates, attributes,
   *          and distribution factors from the mesh file.
   */
  virtual void read_model( mesh::BulkData & ) const = 0 ;

  virtual ~MeshReader();

  const std::vector< mesh::Part * > & parts() const { return m_parts ; }

  const std::vector< std::pair<stk::mesh::FieldBase*,std::string> >
    & attributes() const { return m_attributes ; }

  const std::vector< std::pair<stk::mesh::FieldBase*,std::string> >
    & transients() const { return m_transients ; }

protected:

  ParallelMachine     m_parallel ;
  mesh::MetaData    & m_meta_data ;
  mesh::VectorField & m_coord_field ;
  const std::string & m_file_format ;
  const std::string & m_file_name ;
  std::vector< mesh::Part * > m_parts ;
  std::vector< std::pair<stk::mesh::FieldBase*,std::string> > m_attributes ;
  std::vector< std::pair<stk::mesh::FieldBase*,std::string> > m_transients ;

  MeshReader( ParallelMachine     parallel ,
              mesh::MetaData    & meta_data ,
              mesh::VectorField & coord_field ,
              const std::string & file_format ,
              const std::string & file_name )
    : m_parallel( parallel ),
      m_meta_data( meta_data ),
      m_coord_field( coord_field ),
      m_file_format( file_format ),
      m_file_name( file_name )
    {}

private:
  MeshReader();
  MeshReader( const MeshReader & );
  MeshReader & operator = ( const MeshReader & );
};

//----------------------------------------------------------------------

/** \brief  Filter to apply when reading a mesh file.
 *
 *  Default behavior is to create parts and fields for every
 *  part and field found in the mesh file.
 *
 *  \todo REFACTOR  Differentiate between attribute and non-attribute
 *        fields.  Treat the distribution factor as an attribute field.
 */
class MeshReaderFilter {
public:

  virtual bool accept_part( mesh::MetaData & ,
                            const std::string & ,
                            mesh::EntityRank ,
                            const CellTopologyData * ) const ;

  /** \brief  Map input information to an attribute field.  */
  virtual mesh::FieldBase *
    map_attribute( mesh::Part & , const std::string & ,
                   const std::vector< const shards::ArrayDimTag * > & ) const ;

  /** \brief  Map input information to a transient field.  */
  virtual mesh::FieldBase *
    map_transient( mesh::Part & , const std::string & ,
                   const std::vector< const shards::ArrayDimTag * > & ) const ;

  virtual ~MeshReaderFilter();

  MeshReaderFilter() {}

private:

  MeshReaderFilter( const MeshReaderFilter & );
  MeshReaderFilter & operator = ( const MeshReaderFilter & );
};


} // namespace io
} // namespace stk

#endif /* stk_io_MeshReader_hpp */

