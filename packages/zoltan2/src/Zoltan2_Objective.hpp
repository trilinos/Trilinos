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

/*! \file Zoltan2_Objective.hpp
    \brief The Objective base class and derived classes.

  The Objective is the problem that Zoltan2 has to solve.  It has an
  ObjectSource which is a graph, matrix, or geometry.  It may have
  vertex and/or edge weights.  It may have constraints such as
  fixed vertices.  (TODO: Should we create a Constraints class to
  represent constraints like this?)  It could also have a Result to
  refine rather than starting from scratch.

  A PartitioningObjective specifies number of parts and part sizes.
  A ColoringObjective specifies (TODO).
  An OrderingObjective specifies (TODO).
  A FineGrainMatrixObjective specifies (TODO).
*/

namespace Zoltan2
{

/*! \class Zoltan2::Objective
    \brief Objective is the base class for the problem descriptions.

*/

template<class Scalar, class LNO , class GNO, class AppGID>
  class Objective {

private:

  ObjectSource object;         // TODO some of these are RCPs to objects created
  Result result;               // by the user 
  vector<Scalar> objWeights;
  vector<Scalar> edgeWeights;
  map<GNO, int> fixedVertices;    // vertex GNO -> part

  // TODO: would the problem have parameters?  Not sure.
  //    Certainly the method being called to solve the problem would have parameters 
  //    Zoltan2::RCB(Objective &obj, Teuchos::ParameterList &params, Results &answer)

public:

  /*! Constructors */
  Objective(){}

  /*! Destructor */
  ~Objective(){}

  /*! Copy Constructor */
  Objective(const Objective &os){
  }

  /*! Assignment operator */
  Objective &operator=(const Objective &os){
  }

  // TODO virtual set/get problem components
};

// TODO - concrete objective classes

} // namespace Zoltan2

#endif /* _ZOLTAN2_OBJECTIVE_HPP_ */
