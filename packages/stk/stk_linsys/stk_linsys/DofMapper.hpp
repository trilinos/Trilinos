#ifndef stk_linsys_DofMapper_hpp
#define stk_linsys_DofMapper_hpp

#include <map>

#include <mpi.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_linsys/FieldIdMap.hpp>
#include <stk_linsys/FeiBaseIncludes.hpp>

namespace stk {
/** Linear-System Assembly
*/
namespace linsys {

/** \addtogroup stk_linsys_module
 *  \{
 */

/** Perform mappings between degrees-of-freedom and equation-indices.

A degree-of-freedom is specified by four things:
  (entity-type,entity-id,field,offset-into-field)

An equation-index is a member of a globally contiguous, zero-based index space.

A DOF-mapping allows the caller to provide a degree-of-freedom and obtain an equation-index.
A reverse DOF-mapping allows the caller to provide an equation-index and obtain a degree-of-freedom.

By default this DofMapper class provides DOF-mappings and reverse-DOF-mappings. Providing reverse
DOF-mappings consumes extra memory since it requires constructing an additional FEI object to do
the reverse lookups. If this is not desired, reverse-mappings can be disabled using DofMapper's
second constructor argument.

The FEI library is utilized for accumulating and storing the mappings. (fei::VectorSpace provides
DOF-mappings, and fei::ReverseMapper provides reverse-DOF-mappings.)

Since the FEI works entirely with primitive data types (e.g., int) and has no knowledge
of stk::mesh types, this DofMapper class essentially acts as a translation bridge
between stk::mesh and the FEI library.
*/
class DofMapper {
 public:
  /** Constructor that internally creates an fei::VectorSpace object.*/
  DofMapper(MPI_Comm comm, bool create_reverse_mappings=true);

  /** Constructor that accepts an existing fei::VectorSpace object.*/
//  DofMapper(fei::SharedPtr<fei::VectorSpace> vspace, bool create_reverse_mappings=true);

  /** Destructor */
  virtual ~DofMapper();

  /** Given a mesh, an entity-type and a field, store the resulting DOF mappings.
   *  This method iterates the buckets for the specified entity-type, and for each
   *  bucket that has the given field and is in the specified part-intersection,
   *  DOF-mappings are stored for each entity-id in the bucket.
   *
   *  Usage note: The specified part-intersection may be a single part such as
   *  mesh_bulk.mesh_meta_data().locally_used_part(), although some application
   *  scenarios will involve different parts...
   *
   *  This method may be called repeatedly, to add dof mappings for different parts,
   *  different entity-types, different fields, etc.
   */
  void add_dof_mappings( const stk::mesh::BulkData& mesh_bulk,
                         const stk::mesh::PartVector& part_intersection,
                         stk::mesh::EntityType ent_type,
                         const stk::mesh::FieldBase& field );

  /** This method internally calls fei::VectorSpace::initComplete(), which finalizes
   * and synchronizes the DOF-mappings (ensures that indices for shared-entities are
   * consistent, etc.). Also, if reverse-mappings are not disabled, this method
   * creates the reverse-mappings object. (The get_dof() method is not available until
   * after this has happened.)
   *
   * This is a collective method, must be called on all processors.
   */
  void finalize();

  /** Query whether reverse-DOF-mappings are enabled.
  * (See second constructor argument above.)
  */
  bool reverse_mappings_enabled() const { return m_reverse_mappings_enabled; }

  /** Return the integer id that the specified field is mapped to.
   * The integer id is the FEI's representation of the field.
   */
  int get_field_id(const stk::mesh::FieldBase& field) const;

  /** Return a global equation index for the specified entity type/id pair and field.
   *
   * Note: this method should be const, but it calls an fei query that is not const.
   * When the fei method is corrected, this method will be made const.
   *
   * Note2: this method may not work correctly until after 'finalize()' has been called.
   */
  int get_global_index(stk::mesh::EntityType ent_type,
                       stk::mesh::EntityId ent_id,
                       stk::mesh::FieldBase& field,
                       int offset_into_field=0);

  /** Given a global_index, return the specification for the DOF that it corresponds to.
   * Throw an exception if the global_index is not found, or if DofMapper::finalize() has
   * not been called.
   * Note: this method will be const after the corresponding fei query is corrected for constness.
   */
  void get_dof(int global_index,
               stk::mesh::EntityType& ent_type,
               stk::mesh::EntityId& ent_id,
               const stk::mesh::FieldBase*& field,
               int& offset_into_field) const;

  /** Return the underlying fei::VectorSpace object.
   */
  const fei::SharedPtr<fei::VectorSpace> get_fei_VectorSpace() const { return m_fei_vecspace; }

  /** Return the underlying fei::VectorSpace object.
   */
  fei::SharedPtr<fei::VectorSpace> get_fei_VectorSpace() { return m_fei_vecspace; }

  const FieldIdMap& get_FieldIdMap() const { return m_field_id_map; }

  FieldIdMap& get_FieldIdMap() { return m_field_id_map; }

 private:

  FieldIdMap m_field_id_map;

  //we store the fei::VectorSpace in a fei::SharedPtr because other fei APIs expect it
  //to be that way. e.g., the constructor for fei::MatrixGraph requires a VectorSpace in a
  //fei::SharedPtr...
  fei::SharedPtr<fei::VectorSpace> m_fei_vecspace;
  bool m_reverse_mappings_enabled;
  fei::ReverseMapper* m_fei_reversemap;

  DofMapper(const DofMapper &);
  void operator = (const DofMapper &);
};//class DofMapper

/** \} */

}//namespace linsys
}//namespace stk

#endif

