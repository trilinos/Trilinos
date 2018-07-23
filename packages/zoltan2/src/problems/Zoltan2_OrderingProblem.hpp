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

/*! \file Zoltan2_OrderingProblem.hpp
    \brief Defines the OrderingProblem class.
*/

#ifndef _ZOLTAN2_ORDERINGPROBLEM_HPP_
#define _ZOLTAN2_ORDERINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_OrderingAlgorithms.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_EvaluateOrdering.hpp>

#include <Zoltan2_GraphModel.hpp>
#include <string>
#ifdef HAVE_ZOLTAN2_OVIS
#include <ovis.h>
#endif





#include <bitset>

using Teuchos::rcp_dynamic_cast;

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief OrderingProblem sets up ordering problems for the user.
 *
 *  The OrderingProblem is the core of the Zoltan2 ordering API.
 *  Based on the the user's input and parameters, the OrderingProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the OrderingProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo follow ordering with partitioning
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum, 
 *          we should have the word "View" in function names that return views.
 */

template<typename Adapter>
class OrderingProblem : public Problem<Adapter>
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
  virtual ~OrderingProblem() {}

  OrderingProblem(Adapter *A, ParameterList *p,
                  const RCP<const Teuchos::Comm<int> > &comm) :
    Problem<Adapter>(A, p, comm) 
  {
    HELLO;
    createOrderingProblem();
  }

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  OrderingProblem(Adapter *A, ParameterList *p, MPI_Comm mpicomm) :
  OrderingProblem(A, p, 
                  rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                            Teuchos::opaqueWrapper(mpicomm))))
  {}
#endif

  /*! \brief Constructor that uses a default communicator
   */
  OrderingProblem(Adapter *A, ParameterList *p) : 
  OrderingProblem(A, p, Tpetra::getDefaultComm())
  {}

  /*! \brief Set up validators specific to this Problem
  */
  static void getValidParameters(ParameterList & pl)
  {

#ifdef INCLUDE_ZOLTAN2_EXPERIMENTAL
    AlgND<Adapter>::getValidParameters(pl);
#endif

    RCP<Teuchos::StringValidator> order_method_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "rcm", "minimum_degree", "natural",
          "random", "sorted_degree", "scotch", "nd" )));
    pl.set("order_method", "rcm", "order algorithm",
      order_method_Validator);

    RCP<Teuchos::StringValidator> order_method_type_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "local", "global", "both" )));
    pl.set("order_method_type", "local", "local or global or both",
      order_method_type_Validator);

    RCP<Teuchos::StringValidator> order_package_Validator = Teuchos::rcp(
      new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "amd", "package2", "package3" )));
    pl.set("order_package", "amd", "package to use in ordering",
      order_package_Validator);
  }

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

  //!  \brief Get the local ordering solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  LocalOrderingSolution<lno_t> * getLocalOrderingSolution() {
    if(localOrderingSolution_ == Teuchos::null) {
      throw std::logic_error( "OrderingProblem was not created with local"
        " ordering. Set parameter order_method_type to local or both."
        " Or use getGlobalOrderingSolution()." );
    }
    return setupSolution(localOrderingSolution_);
  }

  //!  \brief Get the global ordering solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  GlobalOrderingSolution<gno_t> * getGlobalOrderingSolution() {
    if(globalOrderingSolution_ == Teuchos::null) {
      throw std::logic_error( "OrderingProblem was not created with global"
        " ordering. Set parameter order_method_type to global or both."
        " Or use getLocalOrderingSolution()." );
    }
    return setupSolution(globalOrderingSolution_);
  }

private:
  template<typename ordering_solution_t>
  ordering_solution_t *setupSolution(RCP<ordering_solution_t> solution) {
    // std::cout << "havePerm= " << solution->havePerm() <<  " haveInverse= "
    //   << solution->haveInverse() << std::endl;
    // Compute Perm or InvPerm, if one is missing.
    if (!(solution->havePerm()))
      solution->computePerm();
    if (!(solution->haveInverse()))
      solution->computeInverse();
    return solution.getRawPtr();
  }

  void createOrderingProblem();

  // local or global ordering is determined by which RCP is NULL
  RCP<LocalOrderingSolution<lno_t> > localOrderingSolution_;
  RCP<GlobalOrderingSolution<gno_t> > globalOrderingSolution_;

};

////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void OrderingProblem<Adapter>::solve(bool updateInputData)
{
  HELLO;

  size_t nVtx = this->baseModel_->getLocalNumObjects();

  // TODO: Assuming one MPI process now. nVtx = ngids = nlids
  try
  {
    std::string method_type = this->params_->template
      get<std::string>("order_method_type", "local");

    if(method_type == "local" || method_type == "both") {
      localOrderingSolution_ = rcp(new LocalOrderingSolution<lno_t>(nVtx));
    }
    if(method_type == "global" || method_type == "both") {
      globalOrderingSolution_ = rcp(new GlobalOrderingSolution<gno_t>(nVtx));
    }
  }
  Z2_FORWARD_EXCEPTIONS;

  // Determine which algorithm to use based on defaults and parameters.
  // TODO: Use rcm if graph model is defined, otherwise use natural.
  // Need some exception handling here, too.

  std::string method = this->params_->template
    get<std::string>("order_method", "rcm");
  
  // TODO: Ignore case
  try
  {

  // could be a template... seems maybe more awkward
  // added this to avoid duplicating local/global below
  // so many times.
  #define ZOLTAN2_COMPUTE_ORDERING                      \
      if(localOrderingSolution_ != Teuchos::null) {     \
        alg.localOrder(localOrderingSolution_);         \
      }                                                 \
      if(globalOrderingSolution_ != Teuchos::null) {    \
        alg.globalOrder(globalOrderingSolution_);       \
      }

  if (method.compare("rcm") == 0) {
    AlgRCM<base_adapter_t> alg(this->graphModel_, this->params_, this->comm_);
    ZOLTAN2_COMPUTE_ORDERING
  }
  else if (method.compare("natural") == 0) {
    AlgNatural<base_adapter_t> alg(this->identifierModel_, this->params_,
      this->comm_);
    ZOLTAN2_COMPUTE_ORDERING
  }
  else if (method.compare("random") == 0) {
    AlgRandom<base_adapter_t> alg(this->identifierModel_, this->params_,
      this->comm_);
    ZOLTAN2_COMPUTE_ORDERING
  }
  else if (method.compare("sorted_degree") == 0) {
    AlgSortedDegree<base_adapter_t> alg(this->graphModel_, this->params_,
      this->comm_);
    ZOLTAN2_COMPUTE_ORDERING
  }
  else if (method.compare("minimum_degree") == 0) {
    std::string pkg = this->params_->template get<std::string>(
      "order_package", "amd");
    if (pkg.compare("amd") == 0)
    {
      AlgAMD<base_adapter_t> alg(this->graphModel_,
        this->params_, this->comm_);
      ZOLTAN2_COMPUTE_ORDERING
    }
  }
  else if (method.compare("scotch") == 0) { // BDD Adding scotch ordering
    AlgPTScotch<Adapter> alg(this->envConst_, this->comm_,
      this->baseInputAdapter_);
    ZOLTAN2_COMPUTE_ORDERING
  }

#ifdef INCLUDE_ZOLTAN2_EXPERIMENTAL
  else if (method.compare("nd") == 0) {
    AlgND<Adapter> alg(this->envConst_, this->comm_, this->graphModel_,
      this->coordinateModel_,this->baseInputAdapter_);
    ZOLTAN2_COMPUTE_ORDERING
  }
#endif

  }
  Z2_FORWARD_EXCEPTIONS;
}

////////////////////////////////////////////////////////////////////////
//template <typename Adapter>
//void OrderingProblem<Adapter>::redistribute()
//{
//  HELLO;
//}

////////////////////////////////////////////////////////////////////////
//! createOrderingProblem 
//  Method with common functionality for creating a OrderingProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter>
void OrderingProblem<Adapter>::createOrderingProblem()
{
  HELLO;
  using Teuchos::ParameterList;

  // Determine which parameters are relevant here.
  // For now, assume parameters similar to Zoltan:
  //   MODEL = graph, hypergraph, geometric, ids
  //   ALGORITHM = rcm, random, amd

  ModelType modelType = IdentifierModelType; //default, change later
  std::string method = this->params_->template
    get<std::string>("order_method", "rcm");

  if ((method == std::string("rcm")) || 
      (method == std::string("sorted_degree")) || 
      (method == std::string("minimum_degree"))) {
    modelType = GraphModelType;
  }

#ifdef INCLUDE_ZOLTAN2_EXPERIMENTAL
  if ((method == std::string("nd")))
  {
    modelType = GraphModelType;
  }

#endif

  // Select Model based on parameters and InputAdapter type

  std::bitset<NUM_MODEL_FLAGS> graphFlags;
  std::bitset<NUM_MODEL_FLAGS> idFlags;


  //MMW: need to change this to allow multiple models
  //     as I did with partitioning, use modelAvail_

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
    this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
      this->baseInputAdapter_, this->envConst_, this->comm_, idFlags));

    this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
      this->identifierModel_);

    break;

  case HypergraphModelType:
  case CoordinateModelType:
    std::cout << __func__zoltan2__ 
              << " Model type " << modelType << " not yet supported." 
              << std::endl;
    break;

  default:
    std::cout << __func__zoltan2__ << " Invalid model" << modelType 
              << std::endl;
    break;
  }
}
} //namespace Zoltan2
#endif
