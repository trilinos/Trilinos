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

namespace Zoltan2 {

  template <typename Adapter>
  class SphynxProblem
  {

  public:

    using part_t = typename Adapter::part_t;
    using weight_t = typename Adapter::scalar_t;

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////// CONSTRUCTORS ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    // Constructor where Teuchos communicator is specified
    SphynxProblem(Teuchos::RCP<Adapter> A, 
		  Teuchos::RCP<Teuchos::ParameterList> p,
		  const Teuchos::RCP<const Teuchos::Comm<int> > &comm):
      inputAdapter_(A),
      params_(p),
      comm_(comm),
      solution_(),
      numberOfWeights_(),
      numberOfCriteria_()
    {

      numberOfWeights_ = this->inputAdapter_->getNumWeightsPerID();
      numberOfCriteria_ = (numberOfWeights_ > 1) ? numberOfWeights_ : 1;

      Teuchos::ArrayRCP<part_t> *noIds = 
	new Teuchos::ArrayRCP<part_t> [numberOfCriteria_];
      Teuchos::ArrayRCP<weight_t> *noSizes = 
	new Teuchos::ArrayRCP<weight_t> [numberOfCriteria_];

      partIds_ = Teuchos::arcp(noIds, 0, numberOfCriteria_, true);
      partSizes_ = Teuchos::arcp(noSizes, 0, numberOfCriteria_, true);
    
      int nparts = -1;
      const Teuchos::ParameterEntry *pe = params_->getEntryPtr("num_global_parts");
      if(pe) 
	nparts = pe->getValue<int>(&nparts);
    
      if(nparts == -1)
	throw std::runtime_error("\nUser did not set num_global_parts"
				 "in the parameter list!n");
      

      envParams_ = Teuchos::rcp(new Teuchos::ParameterList());
      envParams_->set("num_global_parts", nparts);
    
      env_ = Teuchos::rcp(new Environment(*envParams_, comm_));
      envConst_ = Teuchos::rcp_const_cast<const Environment>(env_);

    }

#ifdef HAVE_ZOLTAN2_MPI
    // Constructor where MPI communicator can be specified
    SphynxProblem(Teuchos::RCP<Adapter> A, 
		  Teuchos::RCP<Teuchos::ParameterList> p, 
		  MPI_Comm mpicomm):
      SphynxProblem(A, p, Teuchos::rcp<const Teuchos::Comm<int>>
		    (new Teuchos::MpiComm<int>
		     (Teuchos::opaqueWrapper(mpicomm))))
    {}
#endif

    // Constructor where communicator is the Teuchos default.
    SphynxProblem(Teuchos::RCP<Adapter> A, 
		  Teuchos::RCP<Teuchos::ParameterList> p):
      SphynxProblem(A, p, Tpetra::getDefaultComm())
    {}

    // Destructor
    ~SphynxProblem() {};

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////// FORWARD DECLARATIONS  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    void solve();


    ///////////////////////////////////////////////////////////////////////////
    /////////////////////// MEMBER FUNCTIONS  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

  
    const PartitioningSolution<Adapter> &getSolution() {
      return *(solution_.getRawPtr());
    };  

    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////// DATA MEMBERS ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

  private:
  
    Teuchos::RCP<Adapter> inputAdapter_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Algorithm<Adapter> > algorithm_;

    Teuchos::RCP<Teuchos::ParameterList> envParams_;
    Teuchos::RCP<Environment> env_;
    Teuchos::RCP<const Environment> envConst_;

    Teuchos::RCP<PartitioningSolution<Adapter> > solution_;

    int numberOfWeights_;   // What user provides
    int numberOfCriteria_;  // What Sphynx uses

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<part_t> > partIds_;
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<weight_t> > partSizes_;
  
  };

  ///////////////////////////////////////////////////////////////////////////
  /////////////////////// MORE MEMBER FUNCTIONS  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////

  template <typename Adapter>
  void SphynxProblem<Adapter>::solve()
  {
    this->algorithm_ = Teuchos::rcp(new Zoltan2::Sphynx<Adapter>(this->envConst_,
								 this->params_,
								 this->comm_, 
								 this->inputAdapter_));
  

    PartitioningSolution<Adapter> *soln = NULL;

    try{
    
      soln = new PartitioningSolution<Adapter>(this->envConst_, this->comm_, numberOfWeights_,
					       partIds_.view(0, numberOfCriteria_),
					       partSizes_.view(0, numberOfCriteria_), this->algorithm_);
    }
    Z2_FORWARD_EXCEPTIONS;

    solution_ = Teuchos::rcp(soln);

    // Call the algorithm

    try {
      this->algorithm_->partition(solution_);
    }
    Z2_FORWARD_EXCEPTIONS;
  }


} // namespace Zoltan2

#endif
