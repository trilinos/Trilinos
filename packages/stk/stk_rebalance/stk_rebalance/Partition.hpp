
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
 * \brief For partitioning of mesh entities over a processing grid.
 *
 * This file defines a single class, Partition.  This class describes
 * how a set of mesh entities is distributed over a processing grid.
 * The class contains a list of mesh entities and owning processors.
 * The class is initialized when each processor inserts it's list
 * of mesh entities.  During initialization it is assumed that the
 * current processor is the owning processor so the owning processor
 * list is initialized to the current processor number at that time.
 *
 * A different distribution is defined by changing the owning
 * processor of an entity.
 *
 * The Partition class does not provide advanced communication
 * routines to transfer entities to the new owning processor.
 * The Partition class simply keeps track of the owning
 * processor information.
 *
 * If the mesh entities are redistributed by another class, then the
 * information contained in the Partition class will be out of
 * date and should be cleared and re-initialized.
 */

#include <stdexcept>

// STL components
#include <vector>
#include <utility>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk {
namespace rebalance {

typedef mesh::Field<double, mesh::Cartesian>  VectorField ;
typedef mesh::Field<double>                   ScalarField ;


/** @class Class for keeping track of a mesh entity partition.
 *
 * \brief Initialized with a list of mesh entities unique to each processor.
 *
 * the Field references determine what will be used
 * for the physical coordinate reference and a weight.
 * These fields are made available for the computation
 * of a redistribution by a derived class.
 *
 * The destination processors ids are initialized to the
 * current processor.  Thus the initialization defines
 * a default partition that fits the current mesh entity
 * distribution.  It is then assumed that a new partition
 * will be defined by updating the owing processors.
 * Subsequently the new partition will be realized by
 * a function that actually redistributes the mesh entities
 * to the new owning processors.  But these actions
 * are beyond the scope of this class.
 */

class Partition {

public:

  /** \brief Constructors.  */
  Partition(stk::ParallelMachine comm);
private:
  Partition(const Partition &p);
public:

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
                               const ScalarField   * elem_weight_ref=NULL)
  { throw std::runtime_error("Interface class Partition does not implement set_mesh_info."); }

  /** \brief Destructor. */
  virtual ~Partition();

  /** \brief Return the parallel communicator for this partition entity.*/
  ParallelMachine parallel() const
  { return comm_; }

  /** \brief Return the total number of mesh entities in all lists. */
  virtual unsigned num_elems() const = 0;

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
  virtual void determine_new_partition(bool &RebalancingNeeded) = 0;

  /** \brief Perform communication to create new partition.
   *
   * \param  new_partition  New layout of mesh objects on the processing grid.
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
  virtual int get_new_partition(stk::mesh::EntityProcVec &new_partition) = 0;

  /** \brief Query whether element dependents need to be rebalanced outside this Partition. */
  virtual bool partition_dependents_needed() const = 0;


protected:

  const stk::ParallelMachine comm_;
};

}
} // namespace stk

#endif
