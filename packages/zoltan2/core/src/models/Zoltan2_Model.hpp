// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_Model.hpp
    \brief Defines the Model interface.
*/

#ifndef _ZOLTAN2_MODEL_HPP_
#define _ZOLTAN2_MODEL_HPP_

#include <Zoltan2_Standards.hpp>
#include <bitset>

namespace Zoltan2 {

/*! \brief An identifier for the general type of model.
 */
enum ModelType
{
  HypergraphModelType,
  GraphModelType,
  CoordinateModelType,
  IdentifierModelType,
  MAX_NUM_MODEL_TYPES
};



/*! \brief Flags are set by a Problem to tell a Model what transformations
 *          it may need to do on the user's input.
 */
enum ModelFlags{
  // General flags
  GENERATE_CONSECUTIVE_IDS, /*!< \brief algorithm requires consecutive ids */

  // Graph model flags
  BUILD_LOCAL_GRAPH, /*!< \brief model represents graph within only one rank*/
  SYMMETRIZE_INPUT_TRANSPOSE, /*!< \brief model must symmetrize input */
  SYMMETRIZE_INPUT_BIPARTITE, /*!< \brief model must symmetrize input */
  VERTICES_ARE_MATRIX_ROWS,   /*!< \brief use matrix rows as graph vertices */
  VERTICES_ARE_MATRIX_COLUMNS,/*!< \brief use columns as graph vertices */
  VERTICES_ARE_MATRIX_NONZEROS, /*!< \brief use nonzeros as graph vertices */
  VERTICES_ARE_MESH_NODES,    /*!< \brief use mesh nodes as vertices */
  VERTICES_ARE_MESH_ELEMENTS, /*!< \brief use mesh elements as vertices */
  REMOVE_SELF_EDGES,     /*!< \brief algorithm requires no self edges */
  BUILD_SUBSET_GRAPH,    /*!< \brief ignore invalid neighbors */

  NUM_MODEL_FLAGS
};

typedef std::bitset<NUM_MODEL_FLAGS> modelFlag_t;

/*! \brief The base class for all model classes.

  The Model is the computational model created by a Problem based on
  the user's input data and parameters.  Graphs, hypergraph, and
  collections of geometric coordinates are examples of computational
  models.

  The Problem passes the Model to an algorithm.
  The algorithm queries the Model for input to its calculation.

  \todo Add HypergraphModel, CoordinateModel

*/

template <typename Adapter>
class Model
{
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename Adapter::lno_t       lno_t;
  typedef typename Adapter::gno_t       gno_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
#endif

  /*! Destructor
   */
  virtual ~Model() {};

  /*! Constructor
   */
  Model() {}

  /*!  \brief Return the local number of objects.
   *
   * Return the local number of objects, which may be
   *  vertices, matrix rows, identifiers, coordinates,
   *  or mesh nodes or elements.
   */
  virtual size_t getLocalNumObjects() const = 0;

  /*! \brief Return the global number of objects.
   *
   *  Return the global number of objects, which may be
   *  vertices, matrix rows, identifiers, coordinates,
   *  or mesh nodes or elements.
   */
  virtual global_size_t getGlobalNumObjects() const = 0;

protected:

private:

};

}   //  namespace Zoltan2

#endif
