// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Problem.hpp
    \brief Defines the Problem base class.
*/

#ifndef _ZOLTAN2_PROBLEM_HPP_
#define _ZOLTAN2_PROBLEM_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_CoordinateModel.hpp>

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

  /*! \brief If timer data was collected, print out global data.
   *
   *  If the parameter "timer_output_stream" or "timer_output_file"
   *  was set, then timing statistics are available and will be
   *  printed out to the requested output stream with this call.
   *
   *  All processes in the application must call this, even if
   *  they were not all in the problem communicator.
   *  All timers are reset back to zero after this call.
   *
   *  Timer starts, stops and displays are ignored if Zoltan2
   *  is compiled with Z2_OMIT_ALL_ERROR_CHECKING.
   */
#ifdef Z2_OMIT_ALL_ERROR_CHECKING
  void printTimers() const {return;}
#else
  void printTimers() const
  {
    timer_->printAndResetToZero();
  }
#endif


protected:

  // The Problem is templated on the input adapter.  We interact
  // with the input adapter through the base class interface.  
  // The Model objects are also templated on the input adapter and 
  // are explicitly instantiated for each base input type (vector, 
  // graph, matrix, mesh, identifier list, and coordinate list).

  typedef typename Adapter::base_adapter_t base_adapter_t;

  Adapter* inputAdapter_;
  base_adapter_t *baseInputAdapter_;

  RCP<GraphModel<base_adapter_t> > graphModel_;  
  RCP<IdentifierModel<base_adapter_t> > identifierModel_;  
  RCP<CoordinateModel<base_adapter_t> > coordinateModel_;  

  // Algorithms are passed a base model class, and query
  // the model through the base class interface (graph, hypergraph,
  // identifiers, or coordinates).

  RCP<const Model<base_adapter_t> > baseModel_;  

  RCP<ParameterList> params_;
  RCP<const Comm<int> > comm_;

  // The Problem has a non const Environment object.  This is because
  //   the Problem creates the Environment and may update it before
  //   finally calling the algorithm.

  RCP<Environment> env_;

  // The Problem needs a const version of the Environment.  No other
  //    methods are permitted to change the Environment.

  RCP<const Environment> envConst_;

  // If the user requested timing, this is the TimerManager.

  RCP<TimerManager> timer_;

private:
  void setupProblemEnvironment();

};

#ifdef HAVE_ZOLTAN2_MPI

template <typename Adapter>
  Problem<Adapter>::Problem( Adapter *input, ParameterList *params,
    MPI_Comm comm) : inputAdapter_(input), baseInputAdapter_(),
      graphModel_(), identifierModel_(), baseModel_(),
      params_(RCP<ParameterList>(params,false)), 
      comm_(), env_(), envConst_(), timer_()
{
  HELLO;
  RCP<Teuchos::OpaqueWrapper<MPI_Comm> > wrapper = 
    Teuchos::opaqueWrapper(comm);
  comm_ = rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
  setupProblemEnvironment();
}
#endif

template <typename Adapter>
  Problem<Adapter>::Problem( Adapter *input, ParameterList *params):
    inputAdapter_(input), 
    baseInputAdapter_(dynamic_cast<base_adapter_t *>(input)),
    graphModel_(), identifierModel_(), baseModel_(),
    params_(RCP<ParameterList>(params,false)), comm_(), env_(), envConst_(),
    timer_()
{
  HELLO;
  comm_ = DefaultComm<int>::getComm();
  setupProblemEnvironment();
}


template <typename Adapter>
  void Problem<Adapter>::setupProblemEnvironment()
{
  baseInputAdapter_ = dynamic_cast<base_adapter_t *>(inputAdapter_);

  env_ = rcp(new Environment(*params_, 
    Teuchos::DefaultComm<int>::getComm()));
  envConst_ = rcp_const_cast<const Environment>(env_);

  // Give a timer to the Environment of requested.
  bool isSet;
  string entryValue;

  getParameterValue<string>(*params_, "timing_output_stream", isSet,
    entryValue);

  if (isSet){
    if (entryValue == std::string("std::cout"))
      timer_ = rcp(new TimerManager(comm_, &cout));
    else if (entryValue == std::string("std::cerr"))
      timer_ = rcp(new TimerManager(comm_, &cerr));
  }
  else{
    getParameterValue<string>(*params_, "timing_output_file",
      isSet, entryValue);
    if (isSet){
      std::ofstream *dbgFile = new std::ofstream;
      if (comm_->getRank()==0){
        // Using Teuchos::TimeMonitor, node 0 prints global timing info.
        try{
          dbgFile->open(entryValue.c_str(), std::ios::out|std::ios::trunc);
        }
        catch(std::exception &e){
          throw std::runtime_error(e.what());
        }
      }
      timer_ = rcp(new TimerManager(comm_, dbgFile));
    }
  }

  env_->setTimer(timer_);
}

template <typename Adapter>
  void Problem<Adapter>::resetParameters(ParameterList *params)
{
  env_ = 
    rcp(new Environment(*params, Teuchos::DefaultComm<int>::getComm()));
  envConst_ = rcp_const_cast<const Environment>(env_);
  params_ = rcp<ParameterList>(params,false);

  // We assume the timing output parameters have not changed.

  if (timer_.getRawPtr() != NULL)
    env_->setTimer(timer_);
}

} // namespace Zoltan2

#endif
