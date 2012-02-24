// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Model.hpp

    \brief The abstract interface for a computational model.
*/

#ifndef _ZOLTAN2_MODEL_HPP_
#define _ZOLTAN2_MODEL_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierMap.hpp>
#include <Zoltan2_StridedInput.hpp>

namespace Zoltan2 {

/*! \brief An identifier of the general type of model.
 */
enum ModelType {
  InvalidModel = 0,
  HypergraphModelType,
  GraphModelType,
  GeometryModelType,
  IdentifierModelType
};

/*! \brief Flags that may be used in building a graph model.
 */
enum GraphModelFlags{
  SYMMETRIZE_INPUT_TRANSPOSE,     /*!< model must symmetrize input */
  SYMMETRIZE_INPUT_BIPARTITE,     /*!< model must symmetrize input */
  VERTICES_ARE_MATRIX_ROWS,       /*!< use matrix rows as graph vertices */
  VERTICES_ARE_MATRIX_COLUMNS,    /*!< use matrix columns as vertices */
  VERTICES_ARE_MATRIX_NONZEROS,   /*!< use matrix nonzeros as vertices */
  VERTICES_ARE_MESH_NODES,        /*!< use mesh nodes as vertices */
  VERTICES_ARE_MESH_ELEMENTS,     /*!< use mesh elements as vertices */
  IDS_MUST_BE_GLOBALLY_CONSECUTIVE, /*!< algorithm requires consecutive ids */
  SELF_EDGES_MUST_BE_REMOVED,       /*!< algorithm requires no self edges */
  GRAPH_IS_A_SUBSET_GRAPH, /*!< ignore neighbors that are not valid vertices */
  NUM_GRAPH_MODEL_FLAGS
};

/*! \brief Flags that may be used in building an identifier model.
 */
enum IdentifierModelFlags{
  NEED_CONSECUTIVE_IDS,     /*!< algorithm requires consecutive ids */
  NUM_IDENTIFIER_MODEL_FLAGS
};

/*! \brief The base class for all model classes.
*/

template <typename Adapter>
class Model
{
public:

   /*! Zoltan2 identifier types, user identifier types
    */
  typedef typename Adapter::lno_t    lno_t;
  typedef typename Adapter::gno_t    gno_t;
  typedef typename Adapter::gid_t    gid_t;
  typedef typename Adapter::scalar_t    scalar_t;
  typedef typename Adapter::user_t    user_t;
  typedef StridedInput<lno_t, scalar_t> input_t;
  typedef IdentifierMap<user_t> idmap_t;

  /*! Pure virtual destructor
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

  /*! \brief Get a list of the global Ids for the local objects.
   *
   *  Set a view to the list of object global numbers, which may be
   *  vertex IDs, matrix row IDs, identifiers, coordinate IDs,
   *  or mesh node or element IDs.
   */
  virtual void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const = 0;

  /*! \brief Return the number of weights supplied for each object.
   */
  virtual int getNumWeights() const = 0;

protected:

  /*! Set the IdentifierMap used by the model.
   */
  void setIdentifierMap(RCP<const idmap_t> &map) { idMap_ = map; }

private:

  RCP<const idmap_t> idMap_;
};

}   //  namespace Zoltan2

#endif
