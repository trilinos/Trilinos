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

#ifndef _ZOLTAN2_FUNCTIONS_HPP_
#define _ZOLTAN2_FUNCTIONS_HPP_

/*! \file Zoltan2_Functions.hpp
    \brief The Functions class.
*/

#include <Zoltan2_Objective.hpp>
#include <Zoltan2_Result.hpp>
#include <Teuchos_ParameterList.hpp>

/*! Z2_Interface2
    \brief A namespace for the objects used in interface #2

  Interface #2: The Zoltan methods are functions in the Zoltan2 namespace.
*/

namespace Z2_Interface2
{

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int partition(Z2::PartitioningObjective<Scalar,LNO,GNO,AppGID> &input, 
                Teuchos::ParameterList &params,
                Z2::Interface2 &PartitioningResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int partition(Z2::PartitioningObjective<Scalar,LNO,GNO,AppGID> &input, 
                Z2::Interface2 &PartitioningResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int color(Z2::ColoringObjective<Scalar,LNO,GNO,AppGID> &input, 
                Teuchos::ParameterList &params,
                Z2::Interface2 &ColoringResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int color(Z2::ColoringObjective<Scalar,LNO,GNO,AppGID> &input, 
                Z2::Interface2 &ColoringResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int order(Z2::OrderingObjective<Scalar,LNO,GNO,AppGID> &input, 
                Teuchos::ParameterList &params,
                Z2::Interface2 &OrderingResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int order(Z2::OrderingObjective<Scalar,LNO,GNO,AppGID> &input, 
                Z2::Interface2 &OrderingResult<GNO>);

template <typename Scalar, typename LNO, typename GNO, typename AppGID>
  int evaluate(Z2::Objective<Scalar,LNO,GNO,AppGID> &input, 
                Z2::Interface2 &Result<GNO>);

};

}  // end namespace Z2_Interface2

#endif /* _ZOLTAN2_FUNCTIONS_HPP_ */
