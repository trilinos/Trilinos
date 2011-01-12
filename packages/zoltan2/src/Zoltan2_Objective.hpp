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

#ifndef _ZOLTAN2_OBJECTIVE_HPP_
#define _ZOLTAN2_OBJECTIVE_HPP_

/*!
    \file Zoltan2_Objective.hpp
    \brief The Objective base class and derived classes.
*/

namespace Z2
{

/*! Z2::Objective
    \brief Objective is the base class for the problem descriptions.

  The Objective is the problem that Zoltan2 has to solve.  It has an
  InputAdapter which is a graph, matrix, or geometry.  It may have
  vertex and/or edge weights.  It may have constraints such as
  fixed vertices.  (TODO: Should we create a Constraints class to
  represent constraints like this?)  It could also have a Result to
  refine rather than starting from scratch.
*/

template<typename Scalar, typename LNO , typename GNO, typename AppGID>
  class Objective {

private:

  InputAdapter object;         /*!< The objects to be partitioned, colored, etc. */
  Result in_result;            /*!< An optional starting solution. */
  vector<Scalar> objWeights;   /*!< Optional weighting of objects. */
  vector<Scalar> edgeWeights;  /*!< Optional weighting of edges. */
  map<GNO, int> fixedVertices; /*!< Maybe this should be in a constraints object */

public:

  /*! Constructors */
  Objective(){}

  /*! Destructor */
  ~Objective(){}

  /*! Copy Constructor 
  
    \param o is the initializer of this object.
  */
  Objective(const Objective &o){
  }

  /*! Assignment operator 

    \param o is the right hand side of the copy operator.
  */
  Objective &operator=(const Objective &o){
  }

};

/*! Zoltan2::PartitioningObjective
    \brief A PartitioningObjective contains a partitioning problem.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class PartitioningObjective : public Objective<Scalar, LNO, GNO, AppGID> {
};

/*! Zoltan2::ColoringObjective
    \brief A ColoringObjective contains a coloring problem.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class ColoringObjective : public Objective<Scalar, LNO, GNO, AppGID> {
};

/*! Zoltan2::OrderingObjective
    \brief A OrderingObjective contains a matrix ordering problem.
*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class OrderingObjective : public Objective<Scalar, LNO, GNO, AppGID> {
};

} // namespace Zoltan2

#endif /* _ZOLTAN2_OBJECTIVE_HPP_ */
