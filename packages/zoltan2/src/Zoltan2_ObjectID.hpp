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

#ifndef _ZOLTAN2_OBJECTID_HPP_
#define _ZOLTAN2_OBJECTID_HPP_

/*! \file Zoltan2_ObjectID.hpp
*/

#include <Tpetra_MultiVector_def.hpp>

namespace Zoltan2
{

/*! Zoltan2::ObjectID
    \brief An ObjectID represents a global space of unique object identifiers.

    The IDs may identify geometric coordinates, mesh elements, mesh
    vertices, matrix rows, graph vertices, graph edges, matrix columns,
    or matrix non-zeroes (when 2-dimensional matrix operations are implemented).
    The IDs are supplied by the application which calls Zoltan2.

    The application global ID is some tuple of integral types.  The global ID
    used by Zoltan2 is a singleton integer.  To differentiate between the two
    values we use "global ID" or "GID" for the application ID and "global number"
    or "GNO" for the Zoltan2 ID.

    ObjectID provides mapping between GID, GNO and process local ID (LID).
    It identifies one process with each global ID/global number.  It also identifies one 
    value with each ID.  That value would usually be a part in a partitioning, 
    a color in a coloring, or an order in an ordering.  But in general it is a tuple
    of a caller defined type.

    TODO: We are using the default Kokkos node for now.  It's a Tpetra::MultiVector
       template parameter.

    ObjectID template parameters:
      LNO - the local ID type used by Zoltan2
      GNO - the global ID type used by Zoltan2
      AppGID - the datatype for the application global ID tuple
*/

template<typename LNO, typename GNO, typename AppGID>
  class ObjectID{

private:

  Tpetra::MultiVector<AppGID, LNO, GNO> gid;    /*!< application global IDs */
  Tpetra::Vector<int, LNO, GNO> value;          /*!< part, color or order */

  Tpetra::Map<LNO, GNO> map;                    /*!< distribution of vectors */

public:

  /*! Constructor */
  ObjectID(){}

  /*! Destructor */
  ~ObjectID(){}

  /*! Copy Constructor */
  ObjectID(const ObjectID &id){
  }

  /*! Assignment operator */
  ObjectID &operator=(const ObjectID &id){
  }

  /*! Methods to duplicate functionality of Zoltan_DD_*
   */
  

};

}

#endif /* _ZOLTAN2_OBJECTID_HPP_ */
