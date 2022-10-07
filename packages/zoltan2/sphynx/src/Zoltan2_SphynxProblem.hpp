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
#ifndef _ZOLTAN2_SPHYNXPROBLEM_HPP_
#define _ZOLTAN2_SPHYNXPROBLEM_HPP_


////////////////////////////////////////////////////////////////////////////////
// This file contains the implementation of SphynxProblem.
//
// SphynxProblem is a subset of PartitioningProblem in Zoltan2Core. This subset
// only consists of the functionality and data members needed by Sphynx.
//
// SphynxProblem acts as an interface between user and the Sphynx algorithm.
// User creates the SphynxProblem object on her adapter and calls solve() to
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

static void setSphynxValidatorsInList(
  const Teuchos::ParameterList &plSome,   // in: user's parameters
  const Teuchos::ParameterList &plAll,    // in: validators for all params
  Teuchos::ParameterList &plVal)          // out: validators for user's params
{
  ParameterList::ConstIterator next = plSome.begin();

  while (next != plSome.end()){

    const std::string &name = next->first;
    const ParameterEntry &entrySome = plSome.getEntry(name);
    const ParameterEntry &entryAll = plAll.getEntry(name);

    if (entrySome.isList()){
      plVal.sublist(name);     // create & get
      // Don't set validators for sublists; sublists are for TPL's parameters
    }
    else{
      plVal.setEntry(name, entryAll);
    }

    ++next;
  }
}

  template <typename Adapter>
  class SphynxProblem : public PartitioningProblem<Adapter>
  {

  public:

    using part_t = typename Adapter::part_t;
    using weight_t = typename Adapter::scalar_t;
    typedef typename Adapter::base_adapter_t base_adapter_t; // CHeck to Remove

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////// CONSTRUCTORS ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Constructor where Teuchos communicator is specified
    SphynxProblem(Adapter *A,
                  Teuchos::ParameterList *p,
                            RCP<Teuchos::ParameterList> sphynxParams,
          const RCP<const Teuchos::Comm<int> > &comm):
      PartitioningProblem<Adapter>(A, p, comm), sphynxParams_(sphynxParams)
    {
        // Validation of SphynxParameter
        ParameterList validParams;
        try{
            ParameterList allParameters;
            getSphynxValidParameters(allParameters);

            setSphynxValidatorsInList(*(sphynxParams_.get()), allParameters, validParams);
        }
        Z2_FORWARD_EXCEPTIONS

        sphynxParams_->validateParametersAndSetDefaults(validParams, 0);
        this->env_->convertStringToInt(*sphynxParams_.get());

        int nparts = -1;
        const Teuchos::ParameterEntry *pe = this->params_->getEntryPtr("num_global_parts");
        if(pe)
          nparts = pe->getValue<int>(&nparts);

        if(nparts == -1)
          throw std::runtime_error("\nUser did not set num_global_parts"
                                   "in the parameter list!n");
    }

#ifdef HAVE_ZOLTAN2_MPI
    // Constructor where MPI communicator can be specified
    SphynxProblem(Adapter *A, ParameterList *p, RCP<Teuchos::ParameterList> sphynxParams, MPI_Comm mpicomm):
      SphynxProblem(A, p,sphynxParams,
                              rcp<const Comm<int> >(new Teuchos::MpiComm<int>(
                                                    Teuchos::opaqueWrapper(mpicomm))))
    {}
#endif

    // Constructor where communicator is the Teuchos default.
    SphynxProblem(Adapter *A, ParameterList *p, RCP<Teuchos::ParameterList> sphynxParams):
      SphynxProblem(A, p,sphynxParams, Tpetra::getDefaultComm())
    {}

    // Destructor
    ~SphynxProblem() {};

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// FORWARD DECLARATIONS  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    void createAlgorithm() override;
    void processAlgorithmName(const std::string& algorithm, const std::string& defString, const std::string& model,
                         Environment &env, bool& removeSelfEdges, bool& isGraphType, bool& needConsecutiveGlobalIds) override;

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
    RCP<ParameterList> sphynxParams_;


  };

  ///////////////////////////////////////////////////////////////////////////
  /////////////////////// MORE MEMBER FUNCTIONS  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  template <typename Adapter>
  void SphynxProblem<Adapter>::processAlgorithmName(
      const std::string &algorithm, const std::string &defString,
      const std::string &model, Environment &env, bool &removeSelfEdges,
      bool &isGraphType, bool &needConsecutiveGlobalIds)
  {
    this->algName_ = std::string("sphynx");
  }

  template <typename Adapter>
  void SphynxProblem<Adapter>::createAlgorithm()
  {
      // Create the algorithm
      if (this->algName_ == std::string("sphynx")) {
          this->algorithm_ = Teuchos::rcp(new Zoltan2::Sphynx<Adapter>(this->envConst_,
                                                                       this->params_,
                                                                       this->sphynxParams_,
                                                                       this->comm_,
                                                                       this->inputAdapter_));
      }
      else {
          throw std::logic_error("partitioning algorithm not supported");
      }
  }

} // namespace Zoltan2

#endif
