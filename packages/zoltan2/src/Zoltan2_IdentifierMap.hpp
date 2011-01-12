// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_IDENTIFIERMAP_HPP_
#define _ZOLTAN2_IDENTIFIERMAP_HPP_

/*! \file Zoltan2_IdentifierMap.hpp
    \brief The IdentifierMap class manages global IDs.
*/

#include <Tpetra_MultiVector_def.hpp>

/*! Z2
    \brief A namespace for objects that are common to the 3 proposed interfaces
*/

namespace Z2
{

/*! Z2::IdentifierMap
    \brief An IdentifierMap represents a global space of unique object identifiers.

    The IDs may identify geometric coordinates, mesh elements, mesh
    vertices, matrix rows, graph vertices, graph edges, matrix columns,
    or matrix non-zeroes (when 2-dimensional matrix operations are implemented).
    The IDs are supplied by the application which calls Zoltan2.

    The application global ID is some tuple of integral types.  The global ID
    used by Zoltan2 is a singleton integer.  To differentiate between the two
    values we use "global ID" or "GID" for the application ID and "global number"
    or "GNO" for the Zoltan2 ID.

    IdentifierMap provides mapping between GID, GNO and process local ID (LID).
    It identifies one process with each global ID/global number.  It also identifies one 
    value with each ID.  That value would usually be a part in a partitioning, 
    a color in a coloring, or an order in an ordering.

    IdentifierMap template parameters:
      LNO - the local ID type used by Zoltan2
      GNO - the global ID type used by Zoltan2
      AppGID - the datatype for the application global ID tuple
*/

template<typename LNO, typename GNO, typename AppGID>
  class IdentifierMap{

private:

  Tpetra::MultiVector<AppGID, LNO, GNO> gid;    /*!< application global IDs */
  Tpetra::Vector<int, LNO, GNO> value;          /*!< part, color or order */

  Tpetra::Map<LNO, GNO> map;                    /*!< distribution of vectors */

public:

  /*! Constructor */
  IdentifierMap(){}

  /*! Destructor */
  ~IdentifierMap(){}

  /*! Copy Constructor */
  IdentifierMap(const IdentifierMap &id){
  }

  /*! Assignment operator */
  IdentifierMap &operator=(const IdentifierMap &id){
  }

  /*! Methods to duplicate functionality of Zoltan_DD_*
   */

  void update();
  void find();
  void gid2gno();
  void gno2lno();
  void etc();

  /*! Assign global numbers randomly across processes.
   */
  bool randomize_ownership();
};

}   // end namespace Z2

#endif /* _ZOLTAN2_IDENTIFIERMAP_HPP_ */
