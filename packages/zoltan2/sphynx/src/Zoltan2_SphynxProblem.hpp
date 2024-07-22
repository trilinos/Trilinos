// @HEADER
// *****************************************************************************
//                            Sphynx
//
// Copyright 2020 NTESS and the Sphynx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    RCP<Teuchos::StringValidator> sphynx_eigensolver_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(Teuchos::tuple<std::string>( "LOBPCG", "randomized", "GeneralizedDavidson", "BlockDavidson", "BlockKrylovSchur")));

    pl.set("sphynx_eigensolver", "LOBPCG", "Sphynx eigensolver", sphynx_eigensolver_Validator);

    RCP<Teuchos::StringValidator> sphynx_problem_type_method_Validator =
      Teuchos::rcp( new Teuchos::StringValidator(Teuchos::tuple<std::string>( "combinatorial", "normalized", "generalized")));

    pl.set("sphynx_problem_type", "combinatorial", "Sphynx problem type", sphynx_problem_type_method_Validator);

    RCP<Teuchos::EnhancedNumberValidator<int>> sphynx_verbosity_validator =
      Teuchos::rcp( new Teuchos::EnhancedNumberValidator<int>(0, 1) );
    pl.set("sphynx_verbosity", 0, "Sphynx verbosity.", sphynx_verbosity_validator);

    pl.set("sphynx_ortho_freq", 0, "Sphynx orthogonalization frequency");
    pl.set("sphynx_res_freq", 0, "Sphynx residual frequency");
    pl.set("sphynx_tolerance", 1e-1, "Sphynx tolerance");
    pl.set("sphynx_max_iterations", 1000, "Sphynx max iterations");
    pl.set("sphynx_block_size", 0, "Sphynx block size");
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
      using scalar_t = double; // Sphynx with scalar_t=double obtains better cutsize
      using lno_t = typename Adapter::lno_t;
      using gno_t = typename Adapter::gno_t;
      using node_t = typename Adapter::node_t;
      using mvector_t = typename Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>;  
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
      RCP<mvector_t> getSphynxEigenvectors();
      void setUserEigenvectors(const RCP<mvector_t> &userEvects);

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
      RCP<mvector_t> eigenVectors_;

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
        if( this->eigenVectors_!=Teuchos::null ){
          Teuchos::rcp_dynamic_cast<Zoltan2::Sphynx<Adapter>>(this->algorithm_)->setUserEigenvectors(eigenVectors_);
        }
      }
      else {
        throw std::logic_error("partitioning algorithm not supported");
      }
    }

  ///////////////////////////////////////////////////////////////////////////
  // Allows the user to manually set eigenvectors for the Sphynx partitioner
  // to use rather than solving for them with Anasazi. Mainly intended 
  // for debugging purposes.
  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    void SphynxProblem<Adapter>::setUserEigenvectors(const RCP<mvector_t> &userEvects)
    {
      eigenVectors_ = userEvects; 
    }

  ///////////////////////////////////////////////////////////////////////////
  // Returns an RCP containing a deep copy of the eigenvectors used by Sphynx.
  ///////////////////////////////////////////////////////////////////////////
  template <typename Adapter>
    Teuchos::RCP<Tpetra::MultiVector<double, typename Adapter::lno_t, typename Adapter::gno_t, typename Adapter::node_t> >
    SphynxProblem<Adapter>::getSphynxEigenvectors()
    {
      if(this->algorithm_!=Teuchos::null){
        return Teuchos::rcp_dynamic_cast<Zoltan2::Sphynx<Adapter>>(this->algorithm_)->getSphynxEigenvectors();
      }
      else{
        return Teuchos::null;
      }
    }
} // namespace Zoltan2

#endif
