// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Problem.cpp
    \brief Defines the Problem base class.
*/

#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>

using std::cout;
using std::endl;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.
template<typename Adapter>
class Problem {
public:
  
#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
   */
  Problem(Adapter *, ParameterList *params, MPI_Comm comm);
#endif

  /*! \brief Constructor where communicator is Teuchos default.
   */
  Problem(Adapter *, ParameterList *params);

  /*! \brief Destructor
   */
  virtual ~Problem() {};

  /*! \brief Reset the list of parameters
   */
  void resetParameters(ParameterList *params);

  /*! \brief Method that creates a solution.
   */
  virtual void solve(bool updateInputData) = 0;

protected:
  typedef typename Adapter::base_adapter_t base_adapter_t;

  Adapter* inputAdapter_;
  base_adapter_t *baseInputAdapter_;

  RCP<GraphModel<base_adapter_t> > graphModel_;  
  RCP<IdentifierModel<base_adapter_t> > identifierModel_;  

  RCP<const Model<base_adapter_t> > generalModel_;  

  RCP<ParameterList> params_;
  RCP<const Comm<int> > comm_;

  // The Problem has a non const Environment object.  This is because
  //   the Problem creates the Environment and may update it before
  //   finally calling the algorithm.

  RCP<Environment> env_;

  // The Problem needs a const version of the Environment.  No other
  //    methods are permitted to change the Environment.

  RCP<const Environment> envConst_;

private:

};

#ifdef HAVE_ZOLTAN2_MPI

template <typename Adapter>
  Problem<Adapter>::Problem( Adapter *input, ParameterList *params,
    MPI_Comm comm) : inputAdapter_(input), baseInputAdapter_(),
      graphModel_(), identifierModel_(), generalModel_(),
      params_(RCP<ParameterList>(params,false)), comm_(), env_(), envConst_()
{
  using Teuchos::OpaqueWrapper;
  using Teuchos::opaqueWrapper;

  baseInputAdapter_ = dynamic_cast<base_adapter_t *>(input);

  HELLO;

  env_ = rcp(new Environment(*params, Teuchos::DefaultComm<int>::getComm()));
  envConst_ = rcp_const_cast<const Environment>(env_);

  // The problem communicator may differ from the default application 
  //  communicator in the Environment.

  RCP<OpaqueWrapper<MPI_Comm> > wrapper = opaqueWrapper(comm);
  comm_ = rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
}

#endif

template <typename Adapter>
  Problem<Adapter>::Problem( Adapter *input, ParameterList *params):
    inputAdapter_(input), 
    baseInputAdapter_(dynamic_cast<base_adapter_t *>(input)),
    graphModel_(), identifierModel_(), generalModel_(),
    params_(RCP<ParameterList>(params,false)), comm_(), env_(), envConst_()
{
  HELLO;
  env_ = rcp(new Environment(*params, Teuchos::DefaultComm<int>::getComm()));
  envConst_ = rcp_const_cast<const Environment>(env_);
  comm_ = DefaultComm<int>::getComm();
}

template <typename Adapter>
  void Problem<Adapter>::resetParameters(ParameterList *params)
{
  env_ = rcp(new Environment(*params, Teuchos::DefaultComm<int>::getComm()));
  envConst_ = rcp_const_cast<const Environment>(env_);
  params_ = rcp<ParameterList>(params,false);
}

} // namespace Zoltan2

#endif
