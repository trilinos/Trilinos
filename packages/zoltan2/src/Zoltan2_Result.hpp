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
    \brief The Result class.
*/

/*! Z2
    \brief A namespace for the objects used in interface #1
*/

namespace Z2_Interface1
{

/*! Z2_Interface1::Result
    \brief The Result object contains a solution computed by Zoltan2.

    The Result contains the result of partitioning, ordering or coloring
    the caller's object. 

    It can be queried by the caller.  It can be passed with an Objective
    to a method that will compute a quality metric.  It can be passed with
    an objective to a Zoltan2 method for further computation and improvement.
*/

template<typename GNO>
  class Result {

private:

  bool set;                 /*!< true iff a solution has been computed */
  
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

  void complete() { set = true; }
};

/*! Z2_Interface1::PartitioningResult
    \brief A solution to a partitioning problem.
*/

template<typename GNO>
  class PartitioningResult<GNO> : public Result <GNO>{

private:
  std::map<GNO, int> part_list;

public:
  void get_part_list(std::vector<GNO> &gnos, std::vector<int> &parts);
  void set_part_list(std::vector<GNO> &gnos, std::vector<int> &parts);
};

/*! Z2_Interface1::OrderingResult
    \brief A solution to an ordering problem.
*/

template<typename GNO>
  class OrderingResult<GNO> : public Result <GNO>{

private:
  std::map<GNO, GNO> order;

public:
  void get_order(std::vector<GNO> &gnos, std::vector<GNO> &order);
  void set_order(std::vector<GNO> &gnos, std::vector<GNO> &order);
};

/*! Z2_Interface1::ColoringResult
    \brief A solution to a coloring problem.
*/

template<typename GNO>
  class ColoringResult<GNO> : public Result <GNO>{

private:
  std::map<GNO, int> color;
  int num_colors;
  int base;

public:
  void get_color(std::vector<GNO> &gnos, std::vector<int> &colors);
  void set_color(std::vector<GNO> &gnos, std::vector<int> &colors);
};

}  // end namespace Z2_Interface1

#endif /* _ZOLTAN2_RESULT_HPP_ */
