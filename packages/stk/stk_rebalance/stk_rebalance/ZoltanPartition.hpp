/*----------------------------------------------------------------------*/
/*                                                                      */
/*       author: Jonathan Scott Rath                                    */
/*      author2: Michael W. Glass   (DEC/2000)                          */
/*     filename: ZoltanPartition.h                                      */
/*      purpose: header file for stk toolkit zoltan methods             */
/*                                                                      */
/*----------------------------------------------------------------------*/
/*    Copyright 2001,2010 Sandia Corporation.                           */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a         */
/*    non-exclusive license for use of this work by or on behalf        */
/*    of the U.S. Government.  Export of this program may require       */
/*    a license from the United States Government.                      */
/*----------------------------------------------------------------------*/

// Copyright 2001 Sandia Corporation, Albuquerque, NM.

#ifndef stk_rebalance_ZoltanPartition_hpp
#define stk_rebalance_ZoltanPartition_hpp

#include <utility>
#include <vector>
#include <string>

#include <Teuchos_ParameterList.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_rebalance/GeomDecomp.hpp>

//Forward declaration for pointer to a Zoltan structrue.
struct Zoltan_Struct;

namespace stk {
namespace rebalance {

typedef Teuchos::ParameterList Parameters;

class Zoltan : public GeomDecomp {

public:

  /** MeshInfo is a structure to organize the mesh entity data.
   * A geometric decomposition can be constructed from one or more
   * regions; furthermore, for each region, the decomposition can
   * be based on any type of mesh entity.  Each region has it's
   * own node coordinate field and weight field.
   *
   * Mesh entities are organized according to a vector of MeshInfo
   * structures.  The major index of a mesh entity is the index
   * of the MeshInfo vector and the minor index is the index
   * within the mesh_entities vector.
   *
   * A unique global id is constructed by taking the major index
   * of the MeshInfo vector along with the global_id() of
   * the individual mesh entity.  Together these two integers
   * should form a unique global identification across all
   * processors.
   */

  struct MeshInfo {
    std::vector<mesh::Entity *>      mesh_entities;
    const VectorField * nodal_coord_ref ;
    const ScalarField * elem_weight_ref;
    std::vector<unsigned>            dest_proc_ids ;

    /** Default Constructor. */
    MeshInfo():
      nodal_coord_ref(NULL),
      elem_weight_ref(NULL) {}

    /** Destructor. */
    ~MeshInfo() {}
  };

  virtual void set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                               const VectorField   * nodal_coord_ref,
                               const ScalarField   * elem_weight_ref=NULL);

  /** Reset owning processor.
   *  Default destination for an entity is the processor
   *  that already owns the entity, which is this processor.
   *  The length of the dest_proc_ids vector is the same
   *  length as the mesh_entities vector.
   */
  void reset_dest_proc_data ();

  /** Given mesh entity, find owning processor.
   * This does a search of the internally stored mesh
   * entities and returns what the current owner
   * is set to.  This is internal to the Partition class
   * and may not reflect what the actual owner is as
   * defined by the stk library.
   */
  int proc_owner
  (const mesh::Entity        & ,
   const int                 & index );

  static const std::string zoltan_parameters_name();
  static const std::string default_parameters_name();

  void init_default_parameters();


  /**
   * Constructor
   */

  explicit Zoltan(ParallelMachine pm, const unsigned ndim, Teuchos::ParameterList & rebal_region_parameters, std::string parameters_name=default_parameters_name());

  void init(const std::vector< std::pair<std::string, std::string> >
            &dynamicLoadRebalancingParameters);

  static double init_zoltan_library();
  /**
   * Destructor
   */

  virtual ~Zoltan();

  /** Return the total number of mesh entities in all lists. */
  unsigned num_elems() const
  { return total_number_entities_ ; }

  /** Return the number of local ids per global ids (entities per region).*/
  unsigned num_moid() const;

  /** Various data access functions.
   */
  int globalID         (const unsigned moid) const
  { return mesh_information_.mesh_entities[ moid ]->identifier(); }

  /** Return the owning processor.*/
  unsigned destination_proc(const unsigned moid) const;

  /** Return the owning processor.*/
  double entity_weight(const unsigned moid) const;

  /** Set the owning processor.*/
  void set_destination_proc(const unsigned moid,
                            const unsigned proc );

  /** Find the local ID of a given mesh entity. */
  bool find_mesh_entity(const mesh::Entity * obj, unsigned & moid) const;

  /** Return a mesh entity pointer. */
  mesh::Entity *mesh_entity(const unsigned moid ) const;

  /** Return the Field points to the entity coordinates.*/
  const VectorField * entity_coord_ref () const;

  /** Return the Field points to the entity coordinates.*/
  const ScalarField * entity_weight_ref () const;

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////


  /* Iterator functions.  These functions can be used
   * to iterate over the mesh entities instead of using
   * the access methods above.  There can be only one
   * itererator used at a time since this is not
   * a subclass.
   *
   * It could be made into an iterator
   * subclass by adding a pointer to the
   * Partition class in the private data if
   * there is ever a need for multiple iterators.
   */

  /** Reset iteration to first mesh entity. */
  void iter_init();

  /** Check if at end of entities. */
  bool at_end() const;

  /** Iterator operator. */
  Partition & operator++();

  /** Return current mesh entity. */
  mesh::Entity *iter_mesh_entity() const;

  /** Return current owning processor. */
  unsigned iter_destination_proc() const;

  /** Return the current element weight. */
  double iter_entity_weight() const;

  /** Set the current owning processor. */
  void iter_set_destination_proc (unsigned id);

  /** return current region index. */
  unsigned iter_current_key() const;

  /** Name Conversion Functions.
   * Long friendly string prarameters
   * need to be converted into short Zoltan names.  These
   * two functions do that.  Merge_Default_Values should be called
   * first because it merges with the long names.  Convert_Names_and_Values
   * should then be called second to convert long names to short names.
   */
  static void merge_default_values   (const Teuchos::ParameterList &from,
                                      Teuchos::ParameterList &to);
  static void convert_names_and_values (const Teuchos::ParameterList &from,
                                        Teuchos::ParameterList &to);

  /**
   * Register SIERRA Framework Zoltan call-back functions */
  int register_callbacks();

  virtual void determine_new_partition (bool & RebalancingNeeded);

  virtual int get_new_partition(stk::mesh::EntityProcVec &new_partition);

  bool partition_dependents_needed() const
  { return true; /* Zoltan partitions elements and leaves the rest to someone else */ }

  /**
   * Evaluate the performance/quality of dynamic load rebalancing
   */
  int evaluate ( int   print_stats,
                 int   *nobj,
                 double  *obj_wgt,
                 int   *ncuts,
                 double  *cut_wgt,
                 int   *nboundary,
                 int   *nadj         );

  /**
   * Decomposition Augmentation
   */
  virtual int point_assign( double    *coords,
                            int  *proc ) const;

  virtual int box_assign ( double min[],
                           double max[],
                           std::vector<int> &procs) const;

  /**
   * Inline functions to access private data
   */

  double zoltan_version()  const;
  const std::string & parameter_entry_name() const;

  Zoltan_Struct * zoltan() {
    return zoltan_id;
  }
  const Zoltan_Struct * zoltan() const {
    return zoltan_id;
  }
  unsigned spatial_dimension() const {
    return m_spatial_dimension_;
  }
private:
  /** Zoltan load balancing struct       */
  struct    Zoltan_Struct *zoltan_id;

  const     unsigned       m_spatial_dimension_;
  /** Name that was used to initialize this Zoltan_Struct
   * if the parameter constructor was used.
   */
  std::string  parameter_entry_Name;

  static const std::string zoltanparametersname;
  static const std::string defaultparametersname;
  Parameters               m_default_parameters_;
  MeshInfo  mesh_information_;

  unsigned total_number_entities_;
  unsigned entity_iter_;
  unsigned entity_iter_len_;
  bool iter_initialized_;
};

}
} // namespace stk

#endif
