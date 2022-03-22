// @HEADER
//
// ***********************************************************************
//
//                            Sphynx
//           Copyright 2020 National Technology & Engineering
//                  Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Seher Acer        (sacer@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//                    Karen Devine      (kddevin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_SPHYNXREFACTOREDPROBLEM_HPP_
#define _ZOLTAN2_SPHYNXREFACTOREDPROBLEM_HPP_


////////////////////////////////////////////////////////////////////////////////
// This file contains the implementation of SphynxRefactoredProblem.
//
// SphynxRefactoredProblem is a subset of PartitioningProblem in Zoltan2Core. This subset
// only consists of the functionality and data members needed by Sphynx.
//
// SphynxRefactoredProblem acts as an interface between user and the Sphynx algorithm.
// User creates the SphynxRefactoredProblem object on her adapter and calls solve() to
// get a partitioning solution.
//
////////////////////////////////////////////////////////////////////////////////

#include "Zoltan2_PartitioningSolution.hpp"
#include "Zoltan2_Sphynx.hpp"
#include <Zoltan2_PartitioningProblem.hpp>

namespace Zoltan2 {


/*! \brief Set up validators specific to this algorithm
 */
static void getSphynxValidParameters(ParameterList & pl)
{

  RCP<Teuchos::StringValidator> sphynx_preconditionner_type_method_Validator =
    Teuchos::rcp( new Teuchos::StringValidator(Teuchos::tuple<std::string>( "muelu", "jacobi", "polynomial")));

  pl.set("sphynx_preconditioner_type", "polynomial", "Sphynx preconditioner type", sphynx_preconditionner_type_method_Validator);


  RCP<Teuchos::StringValidator> sphynx_initial_guess_method_Validator =
    Teuchos::rcp( new Teuchos::StringValidator(Teuchos::tuple<std::string>( "random", "constants")));

  pl.set("sphynx_initial_guess", "random", "Sphynx initial guess", sphynx_initial_guess_method_Validator);

  RCP<Teuchos::StringValidator> sphynx_problem_type_method_Validator =
    Teuchos::rcp( new Teuchos::StringValidator(Teuchos::tuple<std::string>( "combinatorial", "normalized", "generalized")));

  pl.set("sphynx_problem_type", "combinatorial", "Sphynx problem type", sphynx_problem_type_method_Validator);

  RCP<Teuchos::EnhancedNumberValidator<int>> sphynx_verbosity_validator =
    Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1) );
  pl.set("sphynx_verbosity", 0, "Sphynx verbosity.", sphynx_verbosity_validator);

  // bool parameter
  pl.set("sphynx_skip_preprocessing", false, "Sphynx skip preprocessing.", Environment::getBoolValidator());
  pl.set("sphynx_use_full_ortho", true, "Sphynx use full ortho.", Environment::getBoolValidator());
}


  template <typename Adapter>
  class SphynxRefactoredProblem : public PartitioningProblem<Adapter>
  {

  public:

    using part_t = typename Adapter::part_t;
    using weight_t = typename Adapter::scalar_t;
    typedef typename Adapter::base_adapter_t base_adapter_t; // CHeck to Remove

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////// CONSTRUCTORS ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Constructor where Teuchos communicator is specified
    SphynxRefactoredProblem(Adapter *A,
                  Teuchos::ParameterList *p,
          const RCP<const Teuchos::Comm<int> > &comm):
      PartitioningProblem<Adapter>(A, p, comm)
    {
//      this->numberOfWeights_ = this->inputAdapter_->getNumWeightsPerID();
//      this->numberOfCriteria_ = (this->numberOfWeights_ > 1) ? this->numberOfWeights_ : 1;

//      Teuchos::ArrayRCP<part_t> *noIds =
//        new Teuchos::ArrayRCP<part_t> [this->numberOfCriteria_];
//      Teuchos::ArrayRCP<weight_t> *noSizes =
//        new Teuchos::ArrayRCP<weight_t> [this->numberOfCriteria_];

//      this->partIds_ = Teuchos::arcp(noIds, 0, this->numberOfCriteria_, true);
//      this->partSizes_ = Teuchos::arcp(noSizes, 0, this->numberOfCriteria_, true);

//        int nparts = -1;
//        const Teuchos::ParameterEntry *pe = params_->getEntryPtr("num_global_parts");
//        if(pe)
//          nparts = pe->getValue<int>(&nparts);

//        if(nparts == -1)
//          throw std::runtime_error("\nUser did not set num_global_parts"
//                                   "in the parameter list!n");


//        envParams_ = Teuchos::rcp(new Teuchos::ParameterList());
//        envParams_->set("num_global_parts", nparts);

//        env_ = Teuchos::rcp(new Environment(*envParams_, comm_));
//        envConst_ = Teuchos::rcp_const_cast<const Environment>(env_);

    }

#ifdef HAVE_ZOLTAN2_MPI
    // Constructor where MPI communicator can be specified
    SphynxRefactoredProblem(Adapter *A, ParameterList *p, MPI_Comm mpicomm):
      SphynxRefactoredProblem(A, p,
                              rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                                    Teuchos::opaqueWrapper(mpicomm))))
    {}
#endif

    // Constructor where communicator is the Teuchos default.
    SphynxRefactoredProblem(Adapter *A, ParameterList *p):
      SphynxRefactoredProblem(A, p, Tpetra::getDefaultComm())
    {}

    // Destructor
    ~SphynxRefactoredProblem() {};

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// FORWARD DECLARATIONS  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

//    void solve();
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

//    void solve(bool updateInputData=true) override;


    void createAlgorithm() override;
    void processAlgorithmName(const std::string& algorithm, const std::string& defString, const std::string& model,
                         Environment &env, bool& removeSelfEdges, bool& needConsecutiveGlobalIds) override;

    ///////////////////////////////////////////////////////////////////////////
    /////////////////////// MEMBER FUNCTIONS  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////


    const PartitioningSolution<Adapter> &getSolution() {
      return *(this->solution_.getRawPtr());
    };

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////// DATA MEMBERS ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

  private:
    Teuchos::RCP<Teuchos::ParameterList> envParams_;

  };

  ///////////////////////////////////////////////////////////////////////////
  /////////////////////// MORE MEMBER FUNCTIONS  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////


  template <typename Adapter>
  void SphynxRefactoredProblem<Adapter>::processAlgorithmName(const std::string& algorithm, const std::string& defString, const std::string& model,
                                                     Environment &env, bool& removeSelfEdges, bool& needConsecutiveGlobalIds) {
      this->algName_ = std::string("sphynx");
      std::cout << "SphynxRefactoredProblem::createAlgorithm " << std::endl;
}

  template <typename Adapter>
  void SphynxRefactoredProblem<Adapter>::createAlgorithm()
  {
      std::cout << "SphynxRefactoredProblem<Adapter>::createAlgorithm" << std::endl;
      // Create the algorithm
      if (this->algName_ == std::string("sphynx")) {
          this->algorithm_ = Teuchos::rcp(new Zoltan2::Sphynx<Adapter>(this->envConst_,
                                                                       this->params_,
                                                                       this->comm_,
                                                                       this->inputAdapter_));
      }
      else {
          throw std::logic_error("partitioning algorithm not supported");
      }
  }

} // namespace Zoltan2

#endif
