// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MappingProblem.hpp
    \brief Defines the MappingProblem class.
*/

#ifndef _ZOLTAN2_MAPPINGPROBLEM_HPP_
#define _ZOLTAN2_MAPPINGPROBLEM_HPP_

#include <Zoltan2_Standards.hpp>

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_MappingSolution.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_MachineRepresentation.hpp>

#include <Zoltan2_AlgBlockMapping.hpp>
#include <Zoltan2_TaskMapping.hpp>
#include <string>

namespace Zoltan2{

////////////////////////////////////////////////////////////////////////

/*! \brief MappingProblem enables mapping of a partition (either 
 *  computed or input) to MPI ranks.
 *
 *  The MappingProblem is the core of the Zoltan2 mappin API.
 *  Based on the the user's input and parameters, the MappingProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the MappingProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 */

template<typename Adapter, 
         typename MachineRep =   // Default MachineRep type
                  MachineRepresentation<typename Adapter::scalar_t,
                                        typename Adapter::part_t> >
class MappingProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

  typedef PartitioningSolution<Adapter> partsoln_t;
  typedef MappingSolution<Adapter> mapsoln_t;

  /*! \brief Destructor
   */
  virtual ~MappingProblem() {};

  /*! \brief Constructor that takes an Teuchos communicator
   */
  MappingProblem(Adapter *A_, Teuchos::ParameterList *p_,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &ucomm_,
                 partsoln_t *partition_ = NULL, 
                 MachineRep *machine_ = NULL) : 
    Problem<Adapter>(A_, p_, ucomm_) 
  {
    HELLO;
    createMappingProblem(partition_, machine_);
  };

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor that takes an MPI communicator
   */
  MappingProblem(Adapter *A_, Teuchos::ParameterList *p_, 
                 MPI_Comm mpicomm_,
                 partsoln_t *partition_ = NULL, 
                 MachineRep *machine_   = NULL) :
  MappingProblem(A_, p_,
                 rcp<const Comm<int> >(
                    new Teuchos::MpiComm<int>(
                        Teuchos::opaqueWrapper(mpicomm_))),
                 partition_, machine_)
  {}
#endif

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
  //  For the sake of performance, we ask the caller to set 
  //  \c updateInputData
  //  to false if he/she is computing a new solution using the same 
  //  input data, but different problem parameters, than that which was 
  //  used to compute the most recent solution.
  
  void solve(bool updateInputData=true); 

  /*! \brief Set up validators specific to this Problem
  */
  static void getValidParameters(ParameterList & pl)
  {
    MachineRepresentation <typename Adapter::scalar_t, 
        typename Adapter::part_t>::getValidParameters(pl);
    RCP<Teuchos::StringValidator> mapping_algorithm_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(
        Teuchos::tuple<std::string>( "geometric", "default", "block" )));
    pl.set("mapping_algorithm", "default", "mapping algorithm",
      mapping_algorithm_Validator);


    // bool parameter
    pl.set("distributed_input_adapter", true,
        "Whether the input adapter for mapping is distributed over processes or not",
        Environment::getBoolValidator());

    // bool parameter
    pl.set("divide_prime_first", false,
        "When partitioning into-non power of two, whether to partition for "
        "nonpowers of two at the beginning, or at the end",
        Environment::getBoolValidator());

    //TODO: This should be positive integer validator.
    pl.set("ranks_per_node", 1,
        "The number of MPI ranks per node",
        Environment::getAnyIntValidator());
    pl.set("reduce_best_mapping", true,
        "If true, nodes will calculate different mappings with rotations, and best "
        "one will be reduced. If not, the result will be the one with longest "
        "dimension partitioning.",
        Environment::getBoolValidator());
  }

  //!  \brief Get the solution to the problem.
  //
  //   \return  the solution to the most recent solve().

  mapsoln_t *getSolution() { return soln.getRawPtr(); };
  Teuchos::RCP<MachineRep> getMachine(){return machine; }
private:
  void createMappingProblem(partsoln_t *partition_, MachineRep *machine_);

  Teuchos::RCP<mapsoln_t> soln;

  Teuchos::RCP<partsoln_t> partition;
  Teuchos::RCP<MachineRep> machine;
};

////////////////////////////////////////////////////////////////////////
//  createMappingProblem 
//  Method with common functionality for creating a MappingProblem.
//  Individual constructors do appropriate conversions of input, etc.
//  This method does everything that all constructors must do.

template <typename Adapter, typename MachineRep>
void MappingProblem<Adapter, MachineRep>::createMappingProblem(
  partsoln_t *partition_,
  MachineRep *machine_)
{
  HELLO;

  // Save pointer to user's partitioning solution.  If not provided, create one.

  if (partition_) {
    // User provided a partitioning solution; use it.
    partition = Teuchos::rcp(partition_, false);
  }
  else {
    // User did not provide a partitioning solution;
    // Use input adapter to create a "fake" solution with the input partition.

    partition = rcp(new partsoln_t(this->env_, this->comm_,
                                   this->inputAdapter_->getNumWeightsPerID()));
    size_t nLocal = this->inputAdapter_->getLocalNumIDs();

    const part_t *inputPartsView = NULL;
    this->inputAdapter_->getPartsView(inputPartsView);
    if (nLocal && inputPartsView == NULL) {
      // User has not provided input parts in input adapter
      int me = this->comm_->getRank();
      ArrayRCP<part_t> inputParts = arcp(new part_t[nLocal], 0, nLocal, true);
      for (size_t i = 0; i < nLocal; i++) inputParts[i] = me;
      partition->setParts(inputParts);
    }
    else {
      // User has provided input parts; use those.
      ArrayRCP<part_t> inputParts = arcp(const_cast<part_t *>(inputPartsView),
                                         0, nLocal, false);
      partition->setParts(inputParts);
    }
  }

  // Save pointer to user's machine.  If not provided, create one.
  if (machine_) 
    machine = Teuchos::rcp(machine_, false);
  else {
    try {
      Teuchos::ParameterList pl = this->env_->getParameters();

      machine = Teuchos::rcp(new MachineRep(*(this->comm_), pl));
    }
    Z2_FORWARD_EXCEPTIONS;
  }
}

////////////////////////////////////////////////////////////////////////
template <typename Adapter, typename MachineRep>
void MappingProblem<Adapter, MachineRep>::solve(bool newData)
{
  HELLO;


  // Determine which algorithm to use based on defaults and parameters.
  std::string algName("block");  

  Teuchos::ParameterList pl = this->env_->getParametersNonConst();
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("mapping_algorithm");
  if (pe) algName = pe->getValue<std::string>(&algName);

  try {
    if (algName == "default") {
      throw(NotImplemented(__FILE__, __LINE__, __func__zoltan2__));
#ifdef KDDKDD_NOT_READH
      this->algorithm_ = rcp(new AlgDefaultMapping<Adapter,MachineRep>(
                                                   this->comm_, machine,
                                                   this->inputAdapter_,
                                                   partition, this->envConst_));
      this->soln = rcp(new mapsoln_t(this->env_, this->comm_, this->algorithm_));
      this->algorithm_->map(this->soln);
#endif
    }
    else if (algName == "block") {
      this->algorithm_ = rcp(new AlgBlockMapping<Adapter,MachineRep>(
                                                 this->comm_, machine,
                                                 this->inputAdapter_,
                                                 partition, this->envConst_));
      this->soln = rcp(new mapsoln_t(this->env_, this->comm_, this->algorithm_));
      this->algorithm_->map(this->soln);
    }
    else if (algName == "geometric") {

      bool is_input_distributed = true;
      const Teuchos::ParameterEntry *pe_input_adapter = 
        pl.getEntryPtr("distributed_input_adapter");
      if (pe_input_adapter) 
        is_input_distributed = pe_input_adapter->getValue<bool>(&is_input_distributed);


      int ranks_per_node = 1;
      pe_input_adapter = pl.getEntryPtr("ranks_per_node");
      if (pe_input_adapter) 
        ranks_per_node = pe_input_adapter->getValue<int>(&ranks_per_node);

      bool divide_prime_first = false;
      pe_input_adapter = pl.getEntryPtr("divide_prime_first");
      if (pe_input_adapter) 
        divide_prime_first = pe_input_adapter->getValue<bool>(&divide_prime_first);

      bool reduce_best_mapping = true;
      pe_input_adapter = pl.getEntryPtr("reduce_best_mapping");
      if (pe_input_adapter) 
        reduce_best_mapping = pe_input_adapter->getValue<bool>(&reduce_best_mapping);
      
      this->algorithm_ = 
            rcp(new CoordinateTaskMapper<Adapter,part_t>(this->comm_,
                                                         machine, 
                                                         this->inputAdapter_,
                                                         partition,
                                                         this->envConst_,
                                                         is_input_distributed, 
                                                         ranks_per_node, 
                                                         divide_prime_first, 
                                                         reduce_best_mapping));

      this->soln = rcp(new mapsoln_t(this->env_, this->comm_, this->algorithm_));

      this->algorithm_->map(this->soln);
    }
    else {
      // Add other mapping methods here
      throw std::logic_error("specified mapping_algorithm not supported");
    }
  }
  Z2_FORWARD_EXCEPTIONS;
}

} //namespace Zoltan2

#endif

#ifdef KDDKDD
Case 1
MappingProblem(
  InputAdapter
  partitioningSolution
  MachineRepresentation=NULL
// KDD Don't know how to properly template MachineRepresentation.  Proper types
// KDD probably depend on how it is to be used.  I imagine MJ needs 
// KDD pcoord_t to be scalar_t, right?  But how does user know that at the
// KDD time he calls this constructor?
)
{
  // Create MachineRepresentation if not provided
  // User would have called partitioning problem and provides a solution
  // Mapping vertices are the parts from the partitioning solution
  // Create MappingSolution that can return getRankForPart(part)
}


Case 2
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Compute mapping vertices based on InputAdapter's partition
  // Assuming input adapter's partition should be used.
  // KDD would use with Exodus/Nemesis input files or PamGen meshes

}


Case 3
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Call a partitioning algorithm to get mapping vertices that are the parts
  // Mapping vertices are computed from this internal partitioning solution
  // Maybe relevant models can be shared.
  // Similar to what's in PartitioningProblem now and to what LibTopoMap does

}


Case 4
MappingProblem(
  InputAdapter
  MachineRepresentation=NULL
)
{
  // Create MachineRepresentation if not provided
  // Call mapping with mapping vertices == IDs from the input adapter.
  // Similar in spirit to previous case, but runs more slowly since current
  // task mapping is done in serial.
  // Mehmet has done experiments with Scotch, comparing case 3 with 4.
  // Case 3 is faster; case 4 has higher quality.


}


// In general, the applyPartitioningSolution method should take an 
// optional MappingSolution.

// Should MappingSolution provide a re-numbered communicator reflecting the new mapping?

#endif
