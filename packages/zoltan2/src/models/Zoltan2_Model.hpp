// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_Model.hpp
    \brief Defines the Model interface.
*/

#ifndef _ZOLTAN2_MODEL_HPP_
#define _ZOLTAN2_MODEL_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierMap.hpp>
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
  IDS_MUST_BE_GLOBALLY_CONSECUTIVE, /*!< \brief algorithm requires consecutive ids */

  // Graph model flags
  SYMMETRIZE_INPUT_TRANSPOSE, /*!< \brief model must symmetrize input */ 
  SYMMETRIZE_INPUT_BIPARTITE, /*!< \brief model must symmetrize input */
  VERTICES_ARE_MATRIX_ROWS,   /*!< \brief use matrix rows as graph vertices */
  VERTICES_ARE_MATRIX_COLUMNS,/*!< \brief use columns as graph vertices */
  VERTICES_ARE_MATRIX_NONZEROS, /*!< \brief use nonzeros as graph vertices */
  VERTICES_ARE_MESH_NODES,    /*!< \brief use mesh nodes as vertices */
  VERTICES_ARE_MESH_ELEMENTS, /*!< \brief use mesh elements as vertices */
  SELF_EDGES_MUST_BE_REMOVED, /*!< \brief algorithm requires no self edges */
  GRAPH_IS_A_SUBSET_GRAPH,    /*!< \brief ignore invalid neighbors */

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
  typedef typename Adapter::zgid_t       zgid_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::user_t      user_t;
  typedef typename Adapter::userCoord_t userCoord_t;
  typedef IdentifierMap<user_t> idmap_t;
#endif

  /*! Destructor
   */
  virtual ~Model() {};

  /*! Constructor
   */
  Model() : idMap_() {}

   /*! \brief Return the map from user global identifiers to internal
   *                Zoltan2 global numbers.
   *
   *  Every model must have an IdentifierMap, whether it needs for mapping 
   *  or not. The Map can simply indicate that Zoltan2 global numbers are 
   *  identical to the application's global IDs.
   */
  RCP<const idmap_t > getIdentifierMap() const { return idMap_; }

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

  /*! \brief Set the IdentifierMap used by the model.
   *
   *  The Model should set the identifier map with this call
   *  during the constructor.
   */
  void setIdentifierMap(RCP<const idmap_t> &map) { idMap_ = map; }

private:

  RCP<const idmap_t> idMap_;

};

}   //  namespace Zoltan2

#endif
