#if 0
// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MatchingProblem.hpp
    \brief Defines the MatchingProblem class.
*/

#ifndef _ZOLTAN2_MATCHINGPROBLEM_HPP_
#define _ZOLTAN2_MATCHINGPROBLEM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_MatchingAlgorithms.hpp>
#include <Zoltan2_MatchingSolution.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>

#include <bitset>

using Teuchos::rcp_dynamic_cast;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief MatchingProblem sets up coloring problems for the user.
 *
 *  The MatchingProblem is the core of the Zoltan2 coloring API.
 *  Based on the the user's input and parameters, the MatchingProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the MatchingProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum, 
 *          we should have the word "View" in function names that return views.
 *
 *  \todo - Currently, only serial and shared-memory coloring are supported.
 */

template<typename Adapter>
class MatchingProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

#ifdef HAVE_ZOLTAN2_MPI
   typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
#endif

  /*! \brief Destructor
   */
  virtual ~MatchingProblem() {};


#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  MatchingProblem(Adapter *A, ParameterList *p, MPI_Comm comm) 
                      : Problem<Adapter>(A, p, comm) 
  {
    HELLO;
    createMatchingProblem();
  };
#endif

  /*! \brief Constructor that uses a default communicator
   */
  MatchingProblem(Adapter *A, ParameterList *p) : Problem<Adapter>(A, p) 
  {
    HELLO;
    createMatchingProblem();
  };

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute
  //  the most recent solution.
  
  void solve(bool updateInputData=true); 

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  MatchingSolution<Adapter> *getSolution() {
    // Get the raw ptr from the rcp
    return solution_.getRawPtr();
  };

private:
  void createMatchingProblem();

  RCP<MatchingSolution<Adapter> > solution_;

};


////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void MatchingProblem<Adapter>::solve(bool newData)
{
  HELLO;

  size_t nVtx = this->baseModel_->getLocalNumObjects();

  try
  {
      this->solution_ = rcp(new MatchingSolution<Adapter>(nVtx));
  }
  Z2_FORWARD_EXCEPTIONS;

  // Determine which algorithm to use based on defaults and parameters.
  // Need some exception handling here, too.

  std::string method = this->params_->template get<std::string>("color_method", "SerialGreedy");

  try
  {
  // TODO: Ignore case
  if (method.compare("SerialGreedy") == 0)
  {
      AlgSerialGreedy<Adapter> alg(this->graphModel_, this->params_,
                                   this->env_, this->comm_);
      alg.color(this->solution_);
  }
#if 0 // TODO later
  else if (method.compare("speculative") == 0) // Gebremedhin-Manne
  {
      AlgGM<base_adapter_t> alg(this->graphModel_, this->comm_);
      alg.color(this->solution_, this->params_);
  }
#endif
  }
  Z2_FORWARD_EXCEPTIONS;

}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void MatchingProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createMatchingProblem 
//  Method with common functionality for creating a MatchingProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void MatchingProblem<Adapter>::createMatchingProblem()
{
  HELLO;
  using Teuchos::ParameterList;

// std::cout << __func__zoltan2__ << " input adapter type " 
//       << this->inputAdapter_->inputAdapterType() << " " 
//       << this->inputAdapter_->inputAdapterName() << std::endl;

  // Create a copy of the user's communicator.

  // Only graph model supported.
  // TODO: Allow hypergraph later?

  ModelType modelType = GraphModelType; 

  // Select Model based on parameters and InputAdapter type

  std::bitset<NUM_MODEL_FLAGS> graphFlags;
  std::bitset<NUM_MODEL_FLAGS> idFlags;

  switch (modelType) {

  case GraphModelType:
    graphFlags.set(REMOVE_SELF_EDGES);
    graphFlags.set(BUILD_LOCAL_GRAPH);
    this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_, graphFlags));

    this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->graphModel_);

    break;


  case IdentifierModelType:
  case HypergraphModelType:
  case CoordinateModelType:
   std::cout << __func__zoltan2__ << " Model type " << modelType << " not yet supported." 
         << std::endl;
    break;

  default:
   std::cout << __func__zoltan2__ << " Invalid model" << modelType << std::endl;
    break;
  }
}
} //namespace Zoltan2

#endif
#endif
