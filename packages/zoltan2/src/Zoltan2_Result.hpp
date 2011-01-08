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

#ifndef _ZOLTAN2_RESULT_HPP_
#define _ZOLTAN2_RESULT_HPP_

/*! \file Zoltan2_Result.hpp
*/

namespace Zoltan2
{

/*! Zoltan2::Result
    \brief The Result object contains a solution computed by Zoltan2.

    The Result contains the result of partitioning, ordering or coloring
    the caller's object. 

    It can be queried by the caller.  It can be passed with an Objective
    to a method that will compute a quality metric.  It can be passed with
    an objective to a Zoltan2 method for further computation and improvement.

    TODO: Enhance methods so they can use the result even if refinement of
      a mesh has occured since the result was computed.

    TODO: some of the instance variables are actually smart pointers to objects
           created elsewhere

    TODO: It's not clear to me yet what the design of this class should be.
      Solutions differ depending on the problem input (ObjectSource) and the
      method applied.  So if an Epetra_RowMatrix is partitioned the result
      could be new Epetra_Maps.  Or it could be a mapping from global row
      ID to process ID.  If a graph read from a matrix market file is colored,
      the result could be a map from vertex global IDs to color numbers. If
      matrix rows are ordered, do we provide the order for the local rows or
      all rows.  Is this all in one class, or do we have subclasses for
      different combinations of input type and zoltan method.
*/

template<typename GNO, typename AppGID>
  class Result {

private:

  bool set;                 /*!< true iff a solution has been computed */
  std::map<AppGID, GNO> gid2gno;   /*!< mapping caller's gid to ours   */
  std::map<GNO, AppGID> gno2gid;   /*!< mapping our gid to the caller's*/
  
public:

  /*! Constructor */
  Result(){}

  /*! Destructor */
  ~Result(){}

  /*! Copy Constructor */
  Result(const Result &r){
  }

  /*! Assignment operator */
  Result &operator=(const Result &r){
  }

  /*    set/get for part information TODO */

  int get_parts_for_gids(vector<AppGID> &gids, vector<int> &parts);
  int get_parts_for_gnos(vector<GNO> &gnos, vector<int> &parts);
  int get_part_for_gid(AppGID &gid);
  int get_part_for_gno(GNO gno);

  /*    set/get for order information TODO */

  /*    set/get for color information TODO */

};

}

#endif /* _ZOLTAN2_RESULT_HPP_ */
