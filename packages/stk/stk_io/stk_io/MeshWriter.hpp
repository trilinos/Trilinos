/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_io_MeshWriter_hpp
#define stk_io_MeshWriter_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>

namespace stk {
namespace io {

//----------------------------------------------------------------------

class MeshWriter ;

/** \brief  Allocate a mesh writer for the given file format.
 *  \param part_map  Parts to output with a name for the output.
 *                   These parts must have a primary entity type.
 *  \param field_map Attribute fields to output with a name for the output.
 *  \param field_map Non-attribute fields to output with a name for the output.
 */
MeshWriter * new_writer(
  ParallelMachine ,
  mesh::MetaData    & meta_data ,
  mesh::VectorField & coord_field ,
  const std::vector< mesh::Part * >      & parts ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & attributes ,
  const std::vector< std::pair<mesh::FieldBase*,std::string> > & transients ,
  const std::string & file_format ,
  const std::string & file_name );

//----------------------------------------------------------------------

class MeshWriter {
public:

  /** \brief  Write the entities, relations, coordinates, attributes,
   *          and distribution factors from the mesh file.
   */
  virtual void write_model( mesh::BulkData & ) const = 0 ;

  /** \brief  Write transient fields tagged with the given time */
  virtual void write_transients( mesh::BulkData & , double time ) const = 0 ;

  virtual ~MeshWriter();

protected:

  ParallelMachine     m_parallel ;
  mesh::MetaData    & m_meta_data ;
  mesh::VectorField & m_coord_field ;
  const std::vector< mesh::Part * > m_parts ;
  const std::vector< std::pair<mesh::FieldBase*,std::string> > m_attributes ;
  const std::vector< std::pair<mesh::FieldBase*,std::string> > m_transients ;
  const std::string   m_file_format ;
  const std::string   m_file_name ;

  MeshWriter(
    ParallelMachine     parallel ,
    mesh::MetaData    & meta_data ,
    mesh::VectorField & coord_field ,
    const std::vector< mesh::Part * > & parts ,
    const std::vector< std::pair<mesh::FieldBase*,std::string> > & attributes ,
    const std::vector< std::pair<mesh::FieldBase*,std::string> > & transients ,
    const std::string & file_format ,
    const std::string & file_name )
    : m_parallel( parallel ),
      m_meta_data( meta_data ),
      m_coord_field( coord_field ),
      m_parts( parts ),
      m_attributes( attributes ),
      m_transients( transients ),
      m_file_format( file_format ),
      m_file_name( file_name )
    {}

private:
  MeshWriter();
  MeshWriter( const MeshWriter & );
  MeshWriter & operator = ( const MeshWriter & );
};

//----------------------------------------------------------------------

} // namespace io
} // namespace stk

#endif /* stk_io_MeshWriter_hpp */

