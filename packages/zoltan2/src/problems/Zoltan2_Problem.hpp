// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
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
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_TimerManager.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Tuple.hpp>
#include <Zoltan2_IntegerRangeList.hpp>

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \brief ProblemRoot allows ptr storage and safe dynamic_cast of all
// problem types.

class ProblemRoot {
  public:
    virtual ~ProblemRoot() {} // required virtual declaration

    // could consider storing comm_ here...
    // this accessor means we can get comm without template upcast first
    virtual RCP<const Comm<int> > getComm() = 0;

   /*! \brief Method that creates a solution.
    */
    virtual void solve(bool updateInputData = true) = 0;
};

////////////////////////////////////////////////////////////////////////
//! \brief Problem base class from which other classes (PartitioningProblem, 
//!        ColoringProblem, OrderingProblem, MatchingProblem, etc.) derive.
     
template<typename Adapter>
class Problem : public ProblemRoot {
public:

  /*! \brief Constructor where Teuchos communicator is specified
   */
  Problem(const Adapter *input, ParameterList *params, 
          const RCP<const Comm<int> > &comm):
    inputAdapter_(rcp(input,false)),
    baseInputAdapter_(rcp(dynamic_cast<const base_adapter_t *>(input), false)),
    graphModel_(),
    identifierModel_(),
    baseModel_(),
    algorithm_(),
    params_(),
    comm_(),
    env_(rcp(new Environment(*params, comm))),
    envConst_(rcp_const_cast<const Environment>(env_)), 
    timer_()
  {
    comm_ = comm->duplicate();
    setupProblemEnvironment(params);
  }

  /*! \brief Destructor
   */
  virtual ~Problem() {};

  /*! \brief Return the communicator used by the problem
   */
  RCP<const Comm<int> > getComm() { return comm_; }

  /*! \brief Reset the list of parameters
   */
  void resetParameters(ParameterList *params);

  /*! \brief Return the communicator passed to the problem
   */

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
    if (!timer_.is_null())
      timer_->printAndResetToZero();
  }
#endif

  // Set up validators which are general to all probloems
  static void getValidParameters(ParameterList & pl)
  {
    // bool parameter
    pl.set("compute_metrics", false, "Compute metrics after computing solution",
      Environment::getBoolValidator());

    RCP<Teuchos::StringValidator> hypergraph_model_type_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "traditional", "ghosting" )));
    pl.set("hypergraph_model_type", "traditional", "construction type when "
      "creating a hypergraph model", hypergraph_model_type_Validator);

    // bool parameter
    pl.set("subset_graph", false, "If \"true\", the graph input is to be "
      "subsetted.  If a vertex neighbor is not a valid vertex, it will be "
      "omitted from the pList.  Otherwise, an invalid neighbor identifier "
      "is considered an error.", Environment::getBoolValidator());

    RCP<Teuchos::StringValidator> symmetrize_input_Validator = Teuchos::rcp(
      new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "no", "transpose", "bipartite" )));
    pl.set("symmetrize_input", "no", "Symmetrize input prior to pList.  "
      "If \"transpose\", symmetrize A by computing A plus ATranspose.  "
      "If \"bipartite\", A becomes [[0 A][ATranspose 0]].",
      symmetrize_input_Validator);

    // these sublists are used for parameters which do not get validated
    pl.sublist("zoltan_parameters");
    pl.sublist("parma_parameters");
  }

  /*! \brief Get the current Environment.
   *   Useful for testing.
   */
  const RCP<const Environment> & getEnvironment() const
  {
    return this->envConst_;
  }

protected:

  // The Problem is templated on the input adapter.  We interact
  // with the input adapter through the base class interface.  
  // The Model objects are also templated on the input adapter and 
  // are explicitly instantiated for each base input type (vector, 
  // graph, matrix, mesh, identifier list, and coordinate list).

  typedef typename Adapter::base_adapter_t base_adapter_t;

  RCP<const Adapter> inputAdapter_;
  RCP<const base_adapter_t> baseInputAdapter_;

  RCP<GraphModel<base_adapter_t> > graphModel_;  
  RCP<IdentifierModel<base_adapter_t> > identifierModel_;  
  RCP<CoordinateModel<base_adapter_t> > coordinateModel_;  

  // Algorithms are passed a base model class, and query
  // the model through the base class interface (graph, hypergraph,
  // identifiers, or coordinates).

  RCP<const Model<base_adapter_t> > baseModel_;  

  // Every problem needs an algorithm, right?
  RCP<Algorithm<Adapter> > algorithm_;

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
  void setupProblemEnvironment(ParameterList *pl);

};

template <typename Adapter>
  void Problem<Adapter>::setupProblemEnvironment(ParameterList * /* params */)
{
  ParameterList &processedParameters = env_->getParametersNonConst();
  params_ = rcp<ParameterList>(&processedParameters, false);

#ifndef Z2_OMIT_ALL_PROFILING
  ParameterList pl = *params_;

  // Give a timer to the Environment if requested.
  bool haveType=false, haveStream=false, haveFile=false;
  int choice = MACRO_TIMERS;   // default timer type

  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("timer_type");

  if (pe){
    choice = pe->getValue<int>(&choice);
    haveType = true;
  }

  TimerType tt = static_cast<TimerType>(choice);

  std::string fname;
  pe = pl.getEntryPtr("timer_output_file");
  if (pe){
    haveFile = true;
    fname = pe->getValue<std::string>(&fname);
    std::ofstream *dbgFile = new std::ofstream;
    if (comm_->getRank()==0){
      // Using Teuchos::TimeMonitor, node 0 prints global timing info.
      try{
        dbgFile->open(fname.c_str(), std::ios::out|std::ios::trunc);
      }
      catch(std::exception &e){
        throw std::runtime_error(e.what());
      }
    }
    timer_ = rcp(new TimerManager(comm_, dbgFile, tt));
  }
  else{
    choice = COUT_STREAM;  // default output stream
    pe = pl.getEntryPtr("timer_output_stream");
    if (pe){
      choice = pe->getValue<int>(&choice);
      haveStream = true;
    }

    OSType outputStream = static_cast<OSType>(choice);

    if (haveStream || haveType){
      if (outputStream == COUT_STREAM)
        timer_ = rcp(new TimerManager(comm_, &std::cout, tt));
      else if (outputStream == CERR_STREAM)
        timer_ = rcp(new TimerManager(comm_, &std::cerr, tt));
      else if (outputStream == NULL_STREAM){
        std::ofstream *of = NULL;
        timer_ = rcp(new TimerManager(comm_, of, tt));
      }
    }
  }

  if (haveType || haveStream || haveFile)
    env_->setTimer(timer_);
  
#endif

}

template <typename Adapter>
  void Problem<Adapter>::resetParameters(ParameterList *params)
{
  env_->resetParameters(*params);
  setupProblemEnvironment(params);

  // We assume the timing output parameters have not changed,
  // and carry on with the same timer.

  if (!timer_.is_null())
    env_->setTimer(timer_);
}

} // namespace Zoltan2

#endif
