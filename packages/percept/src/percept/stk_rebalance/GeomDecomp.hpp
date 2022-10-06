// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.




#ifndef stk_rebalance_GeomDecomp_hpp
#define stk_rebalance_GeomDecomp_hpp

// STL components
#include <vector>
#include <string>
#include <utility>
#include <map>

#include <percept/PerceptMesh.hpp>

#include <percept/stk_rebalance/Partition.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
/** Rebalance API
*/
namespace rebalance {

/** \addtogroup stk_rebalance_module
 *  \{
 */

/** \file GeomDecomp.hpp
 *
 * \brief Geometric support for partitioning of mesh entities.
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


/** \class GeomDecomp 
 *
 * \brief Class for determining the optimal partitioning of mesh entities.
 *
 * Derived from the \a Partition class.
 *
 * The \a GeomDecomp class has no data associated with it, only
 * member functions.  All data is inherited from the \a Partition
 * class.
 *
 * \a GeomDecomp has two functions.  It adds functions to
 * compute geometry information for mesh entities, such as the
 * center point for a mesh entity.  And it defines virtual
 * functions to be specialized by other classes that interface
 * to partitioning packages such as \a Zoltan.
 */

class GeomDecomp: public Partition {
public:
  typedef mesh::Field<double> VectorField;
public:

  GeomDecomp( ParallelMachine comm): Partition(comm) {}

  virtual ~GeomDecomp(){}

  /** \brief Convert a single mesh entity to a single point.
   *
   * \param entity  Entity to take centroid of, element, side, face or node.
   *
   * \param ref  Coordinate field to average, usually defined on the nodes.
   *
   * \param coor Is the output coordinate and is assumed to
   *             be large enough to hold a single coordinate.
   *             Note the maximum needed in all cases is length 3.
   *
   * The \a entity_to_point function is used in the case where a mesh entity is
   * an element with many nodes.  Then something like
   * the element centroid can be used to define
   * a single coordinate point for it.
   *
   */
  static void entity_to_point(stk::mesh::BulkData& bulk_data, const mesh::Entity  & entity,
                            const stk::mesh::FieldBase   & ref,
                            std::vector<double> & coor);

  /** \brief Used to return the nodal entities that \a compute_entity_centroid averages.
   *
   * \param entity  Entity to take coordinates of, element, side, face or node.
   *
   * \param ref  Coordinate field to average, usually defined on the nodes.
   *
   * \param coor Is the output coordinates that entity_to_point would average to 
   *             determine a centroid.
   *
   * The return value is the mesh entities from which the coordinates were obtained.
   */

  static std::vector< mesh::Entity > entity_coordinates(stk::mesh::BulkData& bulk_data, const mesh::Entity     & entity,
                                                               const stk::mesh::FieldBase     & ref,
                                                               std::vector<std::vector<double> >    & coordinates);

  /** \brief  Returns a vector of vectors containing the coordinates of the nodes that were used to compute the centroid.
   *
   * \param entity  Entity to take coordinates of, element, side, face or node.
   *
   * \param ref  Coordinate field to average, usually defined on the nodes.
   *
   * \param coor Is the output coordinated and is assumed to
   *             be large enough to hold a single coordinate.
   *             Note the maximum needed in all cases is length 3.
   *
   * return value is the output coordinates that entity_to_point would average to 
   *             determine a centroid.
   *
   */
  static std::vector<std::vector<double> > compute_entity_centroid(stk::mesh::BulkData& bulk_data, const mesh::Entity     & entity,
                                                               const stk::mesh::FieldBase   & ref,
                                                               std::vector<double> & coor);
  /** \brief Check existence of library entry name on domain library.
   * This is a simple one line convenience function
   */

  static bool confirm ( const std::string &param_set_name );

};

/** \} */

}
} // namespace stk

#endif
