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

/*! \file Zoltan2_MatrixPartitioningProblem.hpp
    \brief Defines the MatrixPartitioningProblem class.
*/

#ifndef _ZOLTAN2_MATRIXPARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_MATRIXPARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
// #include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_MatrixPartitioningAlgs.hpp>
#include <Zoltan2_MatrixPartitioningSolution.hpp>
// #include <Zoltan2_EvaluatePartition.hpp>
// #include <Zoltan2_GraphModel.hpp>
// #include <Zoltan2_IdentifierModel.hpp>
// #include <Zoltan2_IntegerRangeList.hpp>
// #include <Zoltan2_MachineRepresentation.hpp>
// #include <Zoltan2_AlgSerialGreedy.hpp>


// TODO:
//
// Currently, we are implementing this matrix partitioning with several
// constraining assumptions.  We should try to relax several of these.
// In the meantime, we should report errors otherwise
//
// Assumptions:
//              1. 2D Cartesian partitioning (generalize to all 2D, 1D+2D)
//              2. Number of row processes = number of column processes (eventually relax this)
//              3. Number of processors is a square number -- follows from 2
//              4. Number of parts == number of processors 
//              5. Only supports matrix adapter (eventually maybe add graph, hypergraph)



namespace Zoltan2{

/*! \brief MatrixPartitioningProblem sets up partitioning problems for the user.
 *
 *  The MatrixPartitioningProblem is the core of the Zoltan2 partitioning API.
 *  Based on the the user's input and parameters, the MatrixPartitioningProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the MatrixPartitioningProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo repartition given an initial solution
 *  \todo follow partitioning with global or local ordering
 *  \todo allow unsetting of part sizes by passing in null pointers
 *  \todo add a parameter by which user tells us there are no self
 *        edges to be removed.
 *  \todo - Should Problems and Solution have interfaces for returning
 *          views and for returning RCPs?  Or just one?  At a minimum,
 *          we should have the word "View" in function names that return views.
 */

template<typename Adapter>
class MatrixPartitioningProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;


  //TODO FIND WAY TO LIMIT THIS TO ONLY WORK WITH MATRIX ADAPTER FOR NOW

  //! \brief Constructor where Teuchos communicator is specified
  MatrixPartitioningProblem(Adapter *A, ParameterList *p,
                      const RCP<const Teuchos::Comm<int> > &comm):
      Problem<Adapter>(A,p,comm), solution_(),
      inputType_(InvalidAdapterType),
      algName_()
  {
    for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
    initializeProblem();


  }

#ifdef HAVE_ZOLTAN2_MPI
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  /*! \brief Constructor where MPI communicator can be specified
   */
  MatrixPartitioningProblem(Adapter *A, ParameterList *p, MPI_Comm mpicomm):
  MatrixPartitioningProblem(A, p, 
                      rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                            Teuchos::opaqueWrapper(mpicomm))))
  {
  }


  //     Problem<Adapter>(A,p,comm), solution_(),
  //     inputType_(InvalidAdapterType),
  //     algName_()
  // {
  //   for(int i=0;i<MAX_NUM_MODEL_TYPES;i++) modelAvail_[i]=false;
  //   initializeProblem();
  // }
#endif

  //! \brief Constructor where communicator is the Teuchos default.
  MatrixPartitioningProblem(Adapter *A, ParameterList *p):
    MatrixPartitioningProblem(A, p, Tpetra::getDefaultComm())
  {

  }


  /*! \brief Destructor
   */
  ~MatrixPartitioningProblem() {};

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.
  //          By input data we mean coordinates, topology, or weights.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute
  //  the most recent solution.

  void solve(bool updateInputData=true );

  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  const PartitioningSolution<Adapter> &getSolution() 
  {
    return *(solution_.getRawPtr());
  };

  /*! \brief Set up validators specific to this Problem
  */
  static void getValidParameters(ParameterList & pl)
  {
    // Zoltan2_AlgMJ<Adapter>::getValidParameters(pl);
    // AlgPuLP<Adapter>::getValidParameters(pl);
    // AlgPTScotch<Adapter>::getValidParameters(pl);
    // AlgSerialGreedy<Adapter>::getValidParameters(pl);
    // AlgForTestingOnly<Adapter>::getValidParameters(pl);

    // This set up does not use tuple because we didn't have constructors
    // that took that many elements - Tuple will need to be modified and I
    // didn't want to have low level changes with this particular refactor
    // TO DO: Add more Tuple constructors and then redo this code to be
    //  Teuchos::tuple<std::string> algorithm_names( "rcb", "multijagged" ... );
    Array<std::string> algorithm_names(1);
    algorithm_names[0] = "2D Cartesian";
    RCP<Teuchos::StringValidator> algorithm_Validator = Teuchos::rcp(
      new Teuchos::StringValidator( algorithm_names ));
    pl.set("algorithm", "2D Cartesian", "partitioning algorithm",
      algorithm_Validator);


    // TODO: create set of objectives for matrix partitioning:  balance nonzeros, balance rows

    // RCP<Teuchos::StringValidator> partitioning_objective_Validator =
    //   Teuchos::rcp( new Teuchos::StringValidator(
    //    Teuchos::tuple<std::string>( "balance_object_count",
    //      "balance_object_weight", "multicriteria_minimize_total_weight",
    //      "multicriteria_minimize_maximum_weight",
    //      "multicriteria_balance_total_maximum", "minimize_cut_edge_count",
    //      "minimize_cut_edge_weight", "minimize_neighboring_parts",
    //      "minimize_boundary_vertices" )));
    // pl.set("partitioning_objective", "balance_object_weight",
    //   "objective of partitioning", partitioning_objective_Validator);

    pl.set("imbalance_tolerance", 1.1, "imbalance tolerance, ratio of "
      "maximum load over average load", Environment::getAnyDoubleValidator());

    // num_global_parts >= 1
    RCP<Teuchos::EnhancedNumberValidator<int>> num_global_parts_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(
        1, Teuchos::EnhancedNumberTraits<int>::max()) ); // no maximum
    pl.set("num_global_parts", 1, "global number of parts to compute "
      "(0 means use the number of processes)", num_global_parts_Validator);

    // num_local_parts >= 0
    RCP<Teuchos::EnhancedNumberValidator<int>> num_local_parts_Validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(
        0, Teuchos::EnhancedNumberTraits<int>::max()) ); // no maximum
    pl.set("num_local_parts", 0, "number of parts to compute for this "
      "process (num_global_parts == sum of all num_local_parts)", 
      num_local_parts_Validator);

    // RCP<Teuchos::StringValidator> partitioning_approach_Validator =
    //   Teuchos::rcp( new Teuchos::StringValidator(
    //     Teuchos::tuple<std::string>( "partition", "repartition",
    //       "maximize_overlap" )));
    // pl.set("partitioning_approach", "partition", "Partition from scratch, "
    //   "partition incrementally from current partition, of partition from "
    //   "scratch but maximize overlap  with the current partition",
    //   partitioning_approach_Validator);

    //TODO do I need new matrix model

    // RCP<Teuchos::StringValidator> model_Validator = Teuchos::rcp(
    //   new Teuchos::StringValidator(
    //     Teuchos::tuple<std::string>( "hypergraph", "graph",
    //       "geometry", "ids" )));
    // pl.set("model", "graph", "This is a low level parameter. Normally the "
    //   "library will choose a computational model based on the algorithm or "
    //   "objective specified by the user.", model_Validator);

  }

private:
  void initializeProblem();

  void createPartitioningProblem(bool newData);

  RCP<MatrixPartitioningSolution<Adapter> > solution_;

  BaseAdapterType inputType_;

  bool modelAvail_[MAX_NUM_MODEL_TYPES];

  std::string algName_;

};
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
  void MatrixPartitioningProblem<Adapter>::initializeProblem()
{
  HELLO;

  this->env_->debug(DETAILED_STATUS, "MatrixPartitioningProblem::initializeProblem");

  inputType_ = this->inputAdapter_->adapterType();

  if(inputType_ != MatrixAdapterType)
  {
    // TODO: Better error support
    std::cerr << "Error: only matrix adapter type supported" << std::endl;
    return;
  }

  this->env_->memory("After initializeProblem");
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void MatrixPartitioningProblem<Adapter>::solve(bool updateInputData)
{
  std::cout << "MatrixPartitioningProblem solve " << std::endl;


  this->env_->debug(DETAILED_STATUS, "Entering solve");

  // Create the computational model.

  this->env_->timerStart(MACRO_TIMERS, "create problem");

  createPartitioningProblem(updateInputData);

  this->env_->timerStop(MACRO_TIMERS, "create problem");

  // Create the solution. The algorithm will query the Solution
  // for part and weight information. The algorithm will
  // update the solution with part assignments and quality metrics.

  // Create the algorithm
  try 
  {


    //TODO NEED to add if statement back once parameters work
    // if (algName_ == std::string("2D Cartesian")) 
    // {

      // this->algorithm_ = rcp(new AlgMatrix<Adapter>(this->envConst_,
      // 								 this->comm_,
      // 								 this->baseInputAdapter_));


      this->algorithm_ = rcp(new AlgMatrix<Adapter>(this->envConst_,
								 this->comm_,
								 this->baseInputAdapter_));
    // }
    // else {
    //   throw std::logic_error("partitioning algorithm not supported");
    // }
  }
  Z2_FORWARD_EXCEPTIONS;

  // Create the solution
  this->env_->timerStart(MACRO_TIMERS, "create solution");
  MatrixPartitioningSolution<Adapter> *soln = NULL;

  try{
    soln = new MatrixPartitioningSolution<Adapter>(
      this->envConst_, this->comm_,
      this->algorithm_);
  }
  Z2_FORWARD_EXCEPTIONS;

  solution_ = rcp(soln);

  this->env_->timerStop(MACRO_TIMERS, "create solution");
  this->env_->memory("After creating Solution");

  // Call the algorithm

  try
  {
    this->algorithm_->partitionMatrix(solution_);
  }
  Z2_FORWARD_EXCEPTIONS;

  this->env_->debug(DETAILED_STATUS, "Exiting solve");
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void MatrixPartitioningProblem<Adapter>::createPartitioningProblem(bool newData)
{
  this->env_->debug(DETAILED_STATUS,
    "MatrixPartitioningProblem::createPartitioningProblem");

  using Teuchos::ParameterList;

  // A Problem object may be reused.  The input data may have changed and
  // new parameters or part sizes may have been set.
  //
  // Save these values in order to determine if we need to create a new model.

  // bool prevModelAvail[MAX_NUM_MODEL_TYPES];
  // for(int i=0;i<MAX_NUM_MODEL_TYPES;i++)
  // {
  //   prevModelAvail[i] = modelAvail_[i];
  // }


  // for(int i=0;i<MAX_NUM_MODEL_TYPES;i++)
  // {
  //   modelAvail_[i] = false;
  // }


  this->env_->debug(DETAILED_STATUS, "    parameters");
  Environment &env = *(this->env_);
  ParameterList &pl = env.getParametersNonConst();
  const Teuchos::ParameterEntry *pe;
  std::string defString("default");

  // Did the user specify a computational model?
  // Should I allow a model to be created?

  // std::string model(defString);
  // pe = pl.getEntryPtr("model");
  // if (pe)
  //   model = pe->getValue<std::string>(&model);


  // Did the user specify an algorithm?

  std::string algorithm(defString);
  pe = pl.getEntryPtr("algorithm");
  if (pe)
    algorithm = pe->getValue<std::string>(&algorithm);

  // Possible algorithm requirements that must be conveyed to the model:

  ///////////////////////////////////////////////////////////////////
  // Determine algorithm, model, and algorithm requirements.  This
  // is a first pass.  Feel free to change this and add to it.

  if (algorithm != defString)
  {

    // Figure out the model required by the algorithm
    if (algorithm == std::string("2D Cartesian"))
    {
      algName_ = algorithm;
    }
    else
    {
      // Parameter list should ensure this does not happen.
      throw std::logic_error("parameter list algorithm is invalid");
    }
  }

  // Object to be partitioned? (rows, columns, etc)

  // Perhaps should set this

  // std::string objectOfInterest(defString);
  // pe = pl.getEntryPtr("objects_to_partition");
  // if (pe)
  //   objectOfInterest = pe->getValue<std::string>(&objectOfInterest);

    this->env_->debug(DETAILED_STATUS, "createPartitioningProblem done");

}
////////////////////////////////////////////////////////////////////////////////


}  // namespace Zoltan2
#endif
