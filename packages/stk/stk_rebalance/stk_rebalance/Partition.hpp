
/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2002 Sandia Corporation.                        */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#ifndef stk_rebalance_Partition_hpp
#define stk_rebalance_Partition_hpp

/** @file Partition.h
 * @brief For partitioning of mesh objects over a processing grid.
 *
 * This file defines a single class, Partition.  This class describes
 * how a set of mesh objects is distributed over a processing grid.
 * The class contains a list of mesh objects and owning processors.
 * The class is initialized when each processor inserts it's list
 * of mesh objects.  During initialization it is assumed that the
 * current processor is the owning processor so the owning processor
 * list is initialized to the current processor number at that time.
 *
 * A different distribution is defined by changing the owning
 * processor of an object.
 *
 * The Partition class does not provide advanced communication
 * routines to transfer objects to the new owning processor.
 * The Partition class simply keeps track of the owning
 * processor information.
 *
 * If the mesh objects are redistributed by another class, then the
 * information contained in the Partition class will be out of
 * date and should be cleared and re-initialized.
 */


// STL components
#include <vector>
#include <utility>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>

namespace stk {
namespace rebalance {

/** Class for keeping track of a mesh object partition.
 * The class must be initialized with a list of mesh
 * objects unique to each processor.
 *
 * the Field references determine what will be used
 * for the physical coordinate reference and a weight.
 * These fields are made available for the computation
 * of a redistribution by a derived class.
 *
 * The destination processors ids are initialized to the
 * current processor.  Thus the initialization defines
 * a default partition that fits the current mesh object
 * distribution.  It is then assumed that a new partition
 * will be defined by updating the owing processors.
 * Subsequently the new partition will be realized by
 * a function that actually redistributes the mesh objects
 * to the new owning processors.  But these actions
 * are beyond the scope of this class.
 */

class Partition {

public:

  /** RegionInfo is a structure to organize the mesh object data.
   * A geometric decomposition can be constructed from one or more
   * regions; furthermore, for each region, the decomposition can
   * be based on any type of mesh object.  Each region has it's
   * own node coordinate field and weight field.
   *
   * Mesh objects are organized according to a vector of RegionInfo
   * structures.  The major index of a mesh object is the index
   * of the RegionInfo vector and the minor index is the index
   * within the mesh_objects vector.
   *
   * A unique global id is constructed by taking the major index
   * of the RegionInfo vector along with the global_id() of
   * the individual mesh object.  Together these two integers
   * should form a unique global identification across all
   * processors.
   */

  struct RegionInfo {
    std::vector<mesh::Entity *>   mesh_objects;
    const stk::mesh::Field<double>                 * nodal_coord_ref ;
    const stk::mesh::Field<double>                 * elem_weight_ref;
    std::vector<unsigned>             dest_proc_ids ;

    /** Default Constructor. */
    RegionInfo():
      nodal_coord_ref(NULL),
      elem_weight_ref(NULL) {}

    /** Destructor. */
    ~RegionInfo() {}

    /** Copy constructor. */
    RegionInfo( const RegionInfo & r ):
      mesh_objects      ( r.mesh_objects    ),
      nodal_coord_ref( r.nodal_coord_ref ),
      elem_weight_ref( r.elem_weight_ref ),
      dest_proc_ids  ( r.dest_proc_ids   ) {}
  };

  /** Default Constructor.
   */
  Partition():
    total_number_objects_(0),
    iter_initialized_(false) {}

  /** Default Copy Constructor.  */
  inline Partition(const Partition& P):
    total_number_objects_(P.total_number_objects_),
    region_obj_information_(P.region_obj_information_),
    object_iter_(P.object_iter_),
    object_iter_len_(P.object_iter_len_),
    iter_initialized_(P.iter_initialized_) {}

  /** Add another list of mesh objects to the partition.
   * The list of mesh objects is unique to the processor
   * and the owning processor will be assumed to be the
   * current processor calling Add_Mesh.
   *
   * The global id (returned by
   * the global_id() function on the mesh object)
   * must be unique across all of the processors
   * for each call to Add_mesh().
   *
   * Add_Mesh() must be called synchronously across
   * all processors, even if it means passing
   * a zero length array of mesh objects. This is
   * not strictly enforced, but failure to do so
   * could lead to mesh objects without unique
   * global ids.
   *
   * It is important that the length of nodal_coord_ref
   * is consistent across all mesh_objects.  It is
   * not possible to mix two dimensional and three
   * dimensional data.  Also the dimension must be
   * one, two, or three.
   *
   * If unspecified, default element weight is 1 for all objects.
   */
  //void add_mesh ( const std::vector<mesh::Entity *> &mesh_objects,
  //                const stk::mesh::Field<double>   * nodal_coord_ref,
  //                const stk::mesh::Field<double>   * elem_weight_ref=NULL);

  /** Replace old mesh with new one.
   * This is useful if there is only one mesh segment.  Function
   * calls reset_mesh_data followed by a call to add_mesh.
   */
  void replace_mesh ( const std::vector<mesh::Entity *> &mesh_objects,
                      const stk::mesh::Field<double>   * nodal_coord_ref,
                      const stk::mesh::Field<double>   * elem_weight_ref=NULL);

  /** Reset owning processor.
   *  Default destination for an object is the processor
   *  that already owns the object, which is this processor.
   *  The length of the dest_proc_ids vector is the same
   *  length as the mesh_objects vector.
   */
  void reset_dest_proc_data ();

  /** Destructor. */
  virtual ~Partition(){}

  /** Given mesh object, find owning processor.
   * This does a search of the internally stored mesh
   * objects and returns what the current owner
   * is set to.  This is internal to the Partition class
   * and may not reflect what the actual owner is as
   * defined by the stk library.
   */
  int proc_owner
  (const mesh::Entity        & mesh_obj,
   const int                 & index );

  /** Return the number of local ids per global ids (objects per region).*/
  unsigned num_moid() const;

  /** Various data access functions.
   * The indexes are the region id as defined by the
   * sequential number of the call to Add_Mesh.
   * The mesh object id is the sequential index of the
   * mesh object in the vector passed to Add_Mesh.
   *
   * MeshObjects_GlobalID() returns the result of global_id()
   * function of the indexed mesh object.
   */
  int globalID         (const unsigned moid) const;

  /** Return the owning processor.*/
  unsigned destination_proc(const unsigned moid) const;

  /** Return the owning processor.*/
  double object_weight   (const unsigned moid) const;

  /** Set the owning processor.*/
  void set_destination_proc(const unsigned moid,
                            const unsigned proc );

  /** Find the local ID of a given mesh object. */
  bool find_mesh_object(const mesh::Entity * obj, unsigned & moid) const;

  /** Return a mesh object pointer. */
  mesh::Entity *mesh_object(const unsigned moid ) const;

  /** Return the Field points to the object coordinates.*/
  const stk::mesh::Field<double> * object_coord_ref () const;

  /** Return the Field points to the object coordinates.*/
  const stk::mesh::Field<double> * object_weight_ref () const;

  /** Return the total number of mesh objects in all lists. */
  unsigned num_elems() const;

  /** Determine New Partition.
   * This is where all of the real work takes place.  This
   * virtual function should be specialized to determine
   * the new partition.  RebalancingNeeded is set if the new
   * partition is different than the old one.
   */
  virtual int determine_new_partition(bool &RebalancingNeeded)=0;

  /** Perform communication to create new partition.
   * Given a communication specification this
   * function will apply the new partition by
   * transfering the ownership of the registered
   * mesh objects according to the specification
   * determiened by the function Determine_New_Partition.
   * After move_mesh_objects is called, GeomDecomp
   * should be reinitialized with new vectors of
   * mesh objects before rebalancing is performed
   * again.
   */
  int get_new_partition(std::vector<mesh::EntityProc> &new_partition);

private:

  unsigned total_number_objects_;
  RegionInfo  region_obj_information_;



  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  /* Iterator functions.  These functions can be used
   * to iterate over the mesh objects instead of using
   * the access methods above.  There can be only one
   * itererator used at a time since this is not
   * a subclass.
   *
   * It could be made into an iterator
   * subclass by adding a pointer to the
   * Partition class in the private data if
   * there is ever a need for multiple iterators.
   */
private:
  unsigned object_iter_;
  unsigned object_iter_len_;
  bool iter_initialized_;
public:
  /** Reset iteration to first mesh object. */
  void iter_init();

  /** Check if at end of objects. */
  bool at_end() const;

  /** Iterator operator. */
  Partition & operator++();

  /** Return current mesh object. */
  mesh::Entity *iter_mesh_object() const;

  /** Return current owning processor. */
  unsigned iter_destination_proc() const;

  /** Return the current element weight. */
  double iter_object_weight() const;

  /** Set the current owning processor. */
  void iter_set_destination_proc (unsigned id);

  /** return current region index. */
  unsigned iter_current_key() const;

};


}
} // namespace stk

#endif
