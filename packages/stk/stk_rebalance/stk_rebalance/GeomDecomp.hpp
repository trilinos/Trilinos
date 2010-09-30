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
 * @brief Geometric support for partitioning of mesh objects.
 *
 * This file contains a single class GeomDecomp derived from
 * the Partition class.  GeomDecomp adds geometric computations
 * to be used in determining an optimal distribution of mesh
 * objects.  GeomDecomp defines some virtual functions that
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

#include <stk_rebalance/Partition.h>
#include <stk_mesh/Types.h>

namespace stk {
namespace rebalance {

/** Class for determining the optimal partitioning of mesh objects.
 * The GeomDecomp class has no data associated with it, only
 * member functions.  All data is inherited from the Partition
 * class.
 *
 * GeomDecomp has two functions.  It adds functions to
 * compute geometry information for mesh objects, such as the
 * center point for a mesh object.  And it defines virtual
 * functions to be specialized by other classes that interface
 * to partitioning packages such as Zoltan.
 */
class GeomDecomp: public Partition {
private:
  static conststd::string zoltanparametersname;
  static conststd::string defaultparametersname;

public:

  GeomDecomp(): Partition() {}

  virtual ~GeomDecomp(){}

  static conststd::string zoltan_parameters_name();
  static conststd::string default_parameters_name();

  static void init_default_parameters();

  /** Convert a single mesh object to a single point.
   * This is used in the case where a mesh object is
   * an element with many nodes.  Then something like
   * the element centroid can be used to define
   * a single coordinate point for it.
   * coor is the output cooreindate and is assumed to
   * be large enough to hold a single coordinate.
   * Note the maximum needed in all cases is length 3.
   *
   * meshobj_coordinates is used to return the nodal coordinates
   * that compute_obj_centroid averages.  The return value is 
   * the mesh objects from which the coordinates were obtained. 
   * compute_obj_centroid returns a vector of vectors containing
   * the coordinates of the nodes that were used to compute the 
   * centroid.
   */
  static std::vector<const mesh::Entity *> meshobj_coordinates(const mesh::Entity              & obj,
                                                         const stk::mesh::Field<double>Ref             & nodal_coor_ref,
                                                         std::vector<std::vector<Real> >  & coordinates);
  static std::vector<std::vector<Real> > compute_obj_centroid( const mesh::Entity     & obj,
                                                               const stk::mesh::Field<double>Ref    & ref,
                                                               std::vector<Real> & coor);
  static void obj_to_point( const mesh::Entity     & obj,
                            const stk::mesh::Field<double>       & ref,
                            std::vector<Real> &coor);

  /** Check existance of library entry name on domain library.
   * This is a simple one line convieience function
   */

  static bool confirm ( conststd::string &param_set_name );

  /** Convert from a single mesh object to a box.
   *  lo and hi are arrays of length one, two, or three
   *  depending on the dimension of the object.  These
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
                          Real                 lo[],
                          Real                 hi[],
                          const stk::mesh::Field<double>        & ref,
                          Real                 box_expansion_sum );



  /** Determine which processor contains a point.
   * position is the input point in one, two, or three dimensions.
   * proc_id is the output which is which processor the point
   * should be assigned to.
   */
  virtual Int point_assign (Real *position, Int *proc_id) const =0;

  /** Determine which processors contain a box.
   * The input is the extent of a box.  Output is a list
   * of processors which would own regions intersecting the
   * box.
   * The minimum and maximum axis alligned box points
   * are given in min and max.
   */
  virtual Int box_assign   (Real min[],
                            Real max[],
                            std::vector<Int> &procs) const =0;

  int point_assign_all_objs (std::vector<EntityProc>);


  /**
   * Find all the processors which geometrically overlap a given
   * mesh object (returned in the third argument).
   * NOTE:
   * The box_expansion_sum expands the surrounding "box" associated
   * with the mesh object uniformly in each direction by the given sum.
   * Specifically, each component of (x,y,z)_min of the box is REDUCED
   * by the sum, and each component of (x,y,z)_max of the box is
   * INCREASED by the sum.
   */

  void ghost_procs (const stk::mesh::Field<double>  & nodal_coord_ref ,
                    const mesh::Entity         & mesh_obj ,
                    std::vector<Int>           & procs,
                    Real                        box_expansion_sum= 0.0 ) const;


  Int owning_proc (const stk::mesh::Field<double>  & nodal_coord_ref ,
                   const mesh::Entity         & mesh_obj          ) const;

  virtual Diag::Writer &verbose_print(Diag::Writer &dout) const;

};

inline Diag::Writer &operator<<(Diag::Writer &dout, const GeomDecomp &geom_decomp) {
  return geom_decomp.verbose_print(dout);
}

}
} // namespace stk

#endif
