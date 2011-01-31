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

/** \addtogroup stk_rebalance_module
 *  \{
 */

/** \file ZoltanPartition.hpp derived from \a Partition to implement Zoltan based rebalancing.
 *
 * \class Zoltan 
 *
 * \brief Class for implementing Zoltan based rebalancing
 *
 * Derived from the \a GeomDecomp class.
 *
 * The \a Partition class can be derived from to define all different ways of 
 * determining a new partition for rebalancing. This is one example where the
 * derived class uses Zoltan to determine the new partition.
 */
namespace stk {
namespace rebalance {

typedef Teuchos::ParameterList Parameters;

class Zoltan : public GeomDecomp {

public:

  /** \struct MeshInfo
   *
   *  \brief A structure to organize the mesh entity data.
   *
   * \param mesh_entities  Vector of mesh entities to balance
   *                       these can span multiple regons.
   *
   * \param nodal_coord_ref Coordinate reference defined on the
   *                        nodes of the entities in \a mesh_entities.
   *
   * \param elem_weight_ref Element weight used in determining the
   *                        new partition.
   *
   * \param dest_proc_ids   A vector of same length at \a mesh_entities
   *                        that will be filled with the new owner processor.
   *
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
    const VectorField              * nodal_coord_ref ;
    const ScalarField              * elem_weight_ref;
    std::vector<unsigned>            dest_proc_ids ;

    /** \brief Default Constructor. */
    MeshInfo():
      nodal_coord_ref(NULL),
      elem_weight_ref(NULL) {}

    /** \brief Destructor. */
    ~MeshInfo() {}
  };

  /** \brief Define mesh entities to balance.
   *
   * \param mesh_entities    Vector of mesh entities to rebalance
   *
   * \param nodal_coord_ref  Nodal coordinate field to determine new partition
   *                         if using geometric based partitioning.
   *
   * \param elem_weight_ref  Weighting of elements used in defining 
   *                         the new partition. If used, the total element
   *                         weight will be balanced across all of the 
   *                         processors. Can be NULL.
   */

  virtual void set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                               const VectorField   * nodal_coord_ref,
                               const ScalarField   * elem_weight_ref=NULL);

  /** \brief Reset owning processor.
   *
   *  Default destination for an entity is the processor
   *  that already owns the entity, which is this processor.
   *  The length of the dest_proc_ids vector is the same
   *  length as the mesh_entities vector.
   */
  void reset_dest_proc_data ();

  /** \brief Return name of Zoltan parameter block being used. */
  static const std::string zoltan_parameters_name();

  /** \brief Return name of default Zoltan parameter block being used. */
  static const std::string default_parameters_name();


  /** \brief Constructor
   *
   * \param pm     The parallel communicator.
   *
   * \param ndim   The spatial dimention of the mesh being balanced.
   *
   * \param rebal_region_parameters  This is a hierarchial map of strings to strings 
   *                                 that defines the parameters used to
   *                                 initialize \a Zoltan. See the 
   *                                 \a fill_default_values function 
   *                                 for a list of valid parameter names.
   *
   * \param parameters_name  The subset of parameter in \a rebal_region_parameters
   *                         to use.  The default name is \a 'DEFAULT".
   *
   * There are many parameters that effect the wrokings of Zoltan and more are
   * added with each release.  Examine the source code, ZoltanPartition.cpp for
   * a list of the latest releaseed and supported parameters.
   */

  explicit Zoltan(ParallelMachine pm, 
                  const unsigned ndim, 
                  Teuchos::ParameterList & rebal_region_parameters, 
                  std::string parameters_name=default_parameters_name());

  /**
   * \brief Destructor
   */

  virtual ~Zoltan();

  /** \brief Return the total number of mesh entities in all lists. */
  unsigned num_elems() const { return m_total_number_entities_ ; }

  /** \brief Return the owning processor.*/
  double entity_weight(const unsigned moid) const;

  /** \brief Set the owning processor.*/
  void set_destination_proc(const unsigned moid,
                            const unsigned proc );

  /** \brief Various data access functions.
   */
  int globalID         (const unsigned moid) const
  { return m_mesh_information_.mesh_entities[ moid ]->identifier(); }

  /** \brief Return the number of local ids per global ids (entities per region).*/
  unsigned num_moid() const;

  /** \brief Find the local ID of a given mesh entity. */
  bool find_mesh_entity(const mesh::Entity * entity, unsigned & moid) const;

  /** \brief Return a mesh entity pointer. */
  mesh::Entity *mesh_entity(const unsigned moid ) const;

  /** \brief Return the Field points to the entity coordinates.*/
  const VectorField * entity_coord_ref () const;

  /** \brief Return the Field points to the entity coordinates.*/
  const ScalarField * entity_weight_ref () const;

  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

  /** \brief Name Conversion Functions.
   * Long friendly string prarameters
   * need to be converted into short Zoltan names.  The
   * \a merge_default_values and \a convert_names_and_values
   * functions do that.  Merge_Default_Values should be called
   * first because it merges with the long names.  
   */
  static void merge_default_values   (const Teuchos::ParameterList &from,
                                      Teuchos::ParameterList &to);

  /** \brief Name Conversion Functions.
   * \a Convert_Names_and_Values
   * should then be called after \a merge_default_values
   *  to convert long names to short names.
   */
  static void convert_names_and_values (const Teuchos::ParameterList &from,
                                        Teuchos::ParameterList &to);


  /** \brief determine New Partition.
   * 
   * \param RebalancingNeeded  If true, then a new partition 
   *                           has been defined.  If false, the
   *                           new partition is the same as the old
   *                           one and no rebalancing is needed.
   *
   * This is where all of the real work takes place.  This
   * virtual function should be specialized to determine
   * the new partition.  \a RebalancingNeeded is set if the new
   * partition is different than the old one.
   */

  virtual void determine_new_partition (bool & RebalancingNeeded);

  /** \brief Perform communication to create new partition.
   *
   * \param  new_partition  New layout of mesh entities on the processing grid.
   *
   * Given a communication specification this
   * function will apply the new partition by
   * transferring the ownership of the registered
   * mesh entities according to the specification
   * determined by the function \a Determine_New_Partition.
   * After \a move_mesh_entities is called, GeomDecomp
   * should be reinitialized with new vectors of
   * mesh entities before rebalancing is performed
   * again.
   */
  virtual int get_new_partition(stk::mesh::EntityProcVec &new_partition);

  /** \brief Query whether element dependents need to be rebalanced outside this Partition. */
  bool partition_dependents_needed() const
  { return true; /* Zoltan partitions elements and leaves the rest to someone else */ }

  /**
   * \brief Evaluate the performance/quality of dynamic load rebalancing
   *
   * \param print_stats   Zoltan should print some info about the partition.
   *
   * \param nentity          Number of entities partitioned on this processor.
   *
   * \param entity_wgt       Total weight of entities on this processor.
   *
   * \param ncuts         Total number of cuts used in partition
   *
   * \param cut_wgt       Not sure about this one.
   * 
   * \param nboundary     Maybe the number of entities on the processor boundary.
   *
   * \param nadj          Not sure about this one.
   *
   * This is a non-virtual function specific to just the Zoltan interface.
   * It calls through to \a Zoltan_LB_Eval_Balance and \a Zoltan_LB_Eval_Graph
   * to return some information about the Zoltan defined partition.
   *
   */
  int evaluate ( int      print_stats,
                 int     *nentity,
                 double  *entity_wgt,
                 int     *ncuts,
                 double  *cut_wgt,
                 int     *nboundary,
                 int     *nadj  );

  /** \brief Return version of Zoltan linked into executable */
  double zoltan_version()  const;

  /** \brief Return the parameter list that is being used for this partition. */
  const std::string & parameter_entry_name() const;

  /** \brief Zoltan_Struct is an internal Zoltan handle. */
  Zoltan_Struct * zoltan() {
    return m_zoltan_id_;
  }

  /** \brief Zoltan_Struct is an internal Zoltan handle. */
  const Zoltan_Struct * zoltan() const {
    return m_zoltan_id_;
  }

  /** \brief Return the spatial dimention of the entities being rebalanced. */
  unsigned spatial_dimension() const {
    return m_spatial_dimension_;
  }

private:
  /** Zoltan load balancing struct       */
  struct    Zoltan_Struct *m_zoltan_id_;

  const     unsigned       m_spatial_dimension_;
  /** Name that was used to initialize this Zoltan_Struct
   * if the parameter constructor was used.
   */
  std::string              m_parameter_entry_name_;

  static const std::string m_zoltanparametersname_;
  static const std::string m_defaultparametersname_;
  Parameters               m_default_parameters_;
  MeshInfo                 m_mesh_information_;
  unsigned                 m_total_number_entities_;

  void init_default_parameters();
  void init(const std::vector< std::pair<std::string, std::string> >
            &dynamicLoadRebalancingParameters);
  static double init_zoltan_library();

  /** \brief Return the owning processor.*/
  unsigned destination_proc(const unsigned moid) const;

  /** \brief Register SIERRA Framework Zoltan call-back functions */
  int register_callbacks();


};

/** \} */

}
} // namespace stk

#endif
