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

/*! \class Zoltan2::Result
    \brief The Result object contains a solution computed by Zoltan2.

    The Result contains the result of partitioning, ordering or coloring
    the caller's object. 

    It can be queried by the caller.  It can be passed with an Objective
    to a method that will compute a quality metric.  It can be passed to a
    Zoltan2 method for further computation and improvement.

    TODO: Enhance methods so they can use the result even if refinement of
      a mesh has occured since the result was computed.
*/

template<class GNO, class AppGID>
  class Result {

private:

  int num_parts;
  int num_colors;
  
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
#endif /* _ZOLTAN2_RESULT_HPP_ */
