
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
 * @brief For partitioning of mesh entities over a processing grid.
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


/** Class for keeping track of a mesh entity partition.
 * The class must be initialized with a list of mesh
 * entities unique to each processor.
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

  /** Constructors.
   */
  Partition(stk::ParallelMachine comm);
private:
  Partition(const Partition &p);
public:

  /** set MeshInfo
   * This is useful if partitioning is done over a subset of the mesh.
   */
  virtual void set_mesh_info ( const std::vector<mesh::Entity *> &mesh_entities,
                               const VectorField   * nodal_coord_ref,
                               const ScalarField   * elem_weight_ref=NULL)
  { throw std::runtime_error("Interface class Partition does not implement set_mesh_info."); }

  /** Destructor. */
  virtual ~Partition();

  /** Return the parallel communicator for this partition entity.*/
  ParallelMachine parallel() const
  { return comm_; }

  /** Return the total number of mesh entities in all lists. */
  virtual unsigned num_elems() const = 0;

  /** determine New Partition.
   * This is where all of the real work takes place.  This
   * virtual function should be specialized to determine
   * the new partition.  RebalancingNeeded is set if the new
   * partition is different than the old one.
   */
  virtual void determine_new_partition(bool &RebalancingNeeded) = 0;

  /** Perform communication to create new partition.
   * Given a communication specification this
   * function will apply the new partition by
   * transfering the ownership of the registered
   * mesh entities according to the specification
   * determiened by the function Determine_New_Partition.
   * After move_mesh_entities is called, GeomDecomp
   * should be reinitialized with new vectors of
   * mesh entities before rebalancing is performed
   * again.
   */
  virtual int get_new_partition(stk::mesh::EntityProcVec &new_partition) = 0;

  /** Query whether element dependents need to be rebalanced outside this Partition. */
  virtual bool partition_dependents_needed() const = 0;


protected:

  const stk::ParallelMachine comm_;
};

}
} // namespace stk

#endif
