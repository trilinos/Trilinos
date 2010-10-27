/*--------------------------------------------------------------------*/
/*    Copyright 2001, 2002, 2010 Sandia Corporation.                  */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Copyright 2001, 2002 Sandia Corporation, Albuquerque, NM.

#ifndef stk_rebalance_GeomDecomp_hpp
#define stk_rebalance_GeomDecomp_hpp

/** @file GeomDecomp.h
 * @brief Geometric support for partitioning of mesh entities.
 *
 * This file contains a single class GeomDecomp derived from
 * the Partition class.  GeomDecomp adds geometric computations
 * to be used in determining an optimal distribution of mesh
 * entities.  GeomDecomp defines some virtual functions that
 * are to be specialized by the class that actually determines
 * the optimal distribution.
 *
 * Presently, only the support class for Zoltan,
 * specializes the GeomDecomp class.  If another optimal
 * distribution method is to be implemented, it should
 * be based on the GeomDecomp class.  This will enable application
 * codes to include just the GeomDecomp base class without
 * having specialized knowledge of the derived classes.
 */


// STL components
#include <vector>
#include <string>
#include <utility>
#include <map>

#include <stk_rebalance/Partition.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace rebalance {

/** Class for determining the optimal partitioning of mesh entities.
 * The GeomDecomp class has no data associated with it, only
 * member functions.  All data is inherited from the Partition
 * class.
 *
 * GeomDecomp has two functions.  It adds functions to
 * compute geometry information for mesh entities, such as the
 * center point for a mesh entity.  And it defines virtual
 * functions to be specialized by other classes that interface
 * to partitioning packages such as Zoltan.
 */
class GeomDecomp: public Partition {
public:
  typedef mesh::Field<double,mesh::Cartesian> VectorField ;
public:

  GeomDecomp(ParallelMachine comm): Partition(comm) {}

  virtual ~GeomDecomp(){}

  /** Convert a single mesh entity to a single point.
   * This is used in the case where a mesh entity is
   * an element with many nodes.  Then something like
   * the element centroid can be used to define
   * a single coordinate point for it.
   * coor is the output cooreindate and is assumed to
   * be large enough to hold a single coordinate.
   * Note the maximum needed in all cases is length 3.
   *
   * entity_coordinates is used to return the nodal coordinates
   * that compute_obj_centroid averages.  The return value is
   * the mesh entities from which the coordinates were obtained.
   * compute_obj_centroid returns a vector of vectors containing
   * the coordinates of the nodes that were used to compute the
   * centroid.
   */
  static std::vector<const mesh::Entity *> entity_coordinates(const mesh::Entity           & obj,
                                                               const VectorField     & nodal_coor_ref,
                                                               std::vector<std::vector<double> >    & coordinates);
  static std::vector<std::vector<double> > compute_obj_centroid( const mesh::Entity     & obj,
                                                               const VectorField   & ref,
                                                               std::vector<double> & coor);
  static void obj_to_point( const mesh::Entity     & obj,
                            const VectorField & ref,
                            std::vector<double> &coor);

  /** Check existance of library entry name on domain library.
   * This is a simple one line convieience function
   */

  static bool confirm ( const std::string &param_set_name );

  /** Convert from a single mesh entity to a box.
   *  lo and hi are arrays of length one, two, or three
   *  depending on the dimension of the entity.  These
   * two points define the extent of the boy.
   * box_expansion_sum is an amount that will be used to
   * expand the box to make it larger than the smallest
   * containing box.
   * The min and max corners of the axis alligned box
   * are returned in lo and hi which are assumed
   * to be large enough to hold a coordinate point.
   * Note the maximum needed in all cases is length 3.
   */
  static void obj_to_box( const mesh::Entity & obj,
                          double                 lo[],
                          double                 hi[],
                          const VectorField & ref,
                          double                 box_expansion_sum );



  /** Determine which processor contains a point.
   * position is the input point in one, two, or three dimensions.
   * proc_id is the output which is which processor the point
   * should be assigned to.
   */
  virtual int point_assign (double *position, int *proc_id) const =0;

  /** Determine which processors contain a box.
   * The input is the extent of a box.  Output is a list
   * of processors which would own regions intersecting the
   * box.
   * The minimum and maximum axis alligned box points
   * are given in min and max.
   */
  virtual int box_assign   (double min[],
                            double max[],
                            std::vector<int> &procs) const =0;

  /**
   * Find all the processors which geometrically overlap a given
   * mesh entity (returned in the third argument).
   * NOTE:
   * The box_expansion_sum expands the surrounding "box" associated
   * with the mesh entity uniformly in each direction by the given sum.
   * Specifically, each component of (x,y,z)_min of the box is REDUCED
   * by the sum, and each component of (x,y,z)_max of the box is
   * INCREASED by the sum.
   */

  void ghost_procs (const VectorField & nodal_coord_ref ,
                    const mesh::Entity         & mesh_obj ,
                    std::vector<int>           & procs,
                    double                        box_expansion_sum= 0.0 ) const;


  int owning_proc (const VectorField & nodal_coord_ref ,
                   const mesh::Entity         & mesh_obj          ) const;

};
}
} // namespace stk

#endif
