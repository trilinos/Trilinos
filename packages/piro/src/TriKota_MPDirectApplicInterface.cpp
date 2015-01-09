// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "TriKota_MPDirectApplicInterface.hpp"
#include "DakotaModel.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace Dakota;
typedef EpetraExt::ModelEvaluator EEME;

TriKota::MPDirectApplicInterface::
MPDirectApplicInterface(
  ProblemDescDB& problem_db_,
  const Teuchos::RCP<Piro::Epetra::StokhosMPSolver>& model_,
  int p_index_, int g_index_) : 
  Dakota::DirectApplicInterface(problem_db_),
  model(model_),
  p_index(p_index_),
  g_index(g_index_),
  orientation(EEME::DERIV_MV_BY_COL)
{
  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (model != Teuchos::null) {

    // Make sure we support multi-point
    TEUCHOS_TEST_FOR_EXCEPTION(
      !model->createInArgs().supports(EEME::IN_ARG_p_mp, p_index), 
      std::logic_error,
      "Model does not support multi-point parameter " << g_index << "!");
    TEUCHOS_TEST_FOR_EXCEPTION(
      !model->createOutArgs().supports(EEME::OUT_ARG_g_mp, g_index), 
      std::logic_error,
      "Model does not support multi-point response " << g_index << "!");

    // Create product vectors 
    model_p = model->create_p_mp(p_index);
    model_g = model->create_g_mp(g_index);
    numParameters = (*model_p)[0].GlobalLength();
    numResponses  = (*model_g)[0].GlobalLength();
    block_size = model_p->map()->NumMyElements();

    *out << "TriKota:: ModeEval has " << numParameters <<
            " parameters and " << numResponses << " responses." << std::endl;
    
    // Create dg/dp
    EEME::DerivativeSupport supportDgDp = 
      model->createOutArgs().supports(EEME::OUT_ARG_DgDp_mp, g_index, p_index);
    supportsSensitivities =
      supportDgDp.supports(EEME::DERIV_TRANS_MV_BY_ROW) ||
      supportDgDp.supports(EEME::DERIV_MV_BY_COL);
    if (supportsSensitivities) {
      *out << "TriKota:: ModeEval supports gradients calculation." << std::endl;
      if (supportDgDp.supports(EEME::DERIV_TRANS_MV_BY_ROW)) {
        orientation = EEME::DERIV_TRANS_MV_BY_ROW;
        model_dgdp = model->create_p_mv_mp(p_index, numResponses);
      }
      else {
        orientation = EEME::DERIV_MV_BY_COL;
        model_dgdp = model->create_g_mv_mp(g_index, numParameters);
      }
    }

    *out << "TriKota:: Setting initial guess from Model Evaluator to Dakota " 
	 << std::endl;
    Model& first_model = *(problem_db_.model_list().begin());
    unsigned int num_dakota_vars =  first_model.acv();
    Dakota::RealVector drv(num_dakota_vars);
    TEUCHOS_TEST_FOR_EXCEPTION(
      num_dakota_vars > numParameters, std::logic_error,
      "TriKota Adapter Error: number of parameters in ModelEvaluator  "
      <<  numParameters 
      << "\n is less then the number of continuous variables\n"
      << " specified in the dakota.in input file " << num_dakota_vars << "\n" );
    for (unsigned int i=0; i<num_dakota_vars; i++) 
      drv[i] = (*model_p)[0][i];
    first_model.continuous_variables(drv);

  }
  else {
    *out << "Warning in TriKota::MPDirectmodellicInterface constructor\n" 
         << "\tModelEvaluator is null. This is OK iff Dakota has assigned"
         << " MPI_COMM_NULL to this Proc " << std::endl;
  }
}


int 
TriKota::MPDirectApplicInterface::
derived_map_ac(const Dakota::String& ac_name)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "derived_map_ac() called with multi-point interface!" << std::endl <<
    " Make sure the asynch keyword is in the dakota input spec!");
}

void
TriKota::MPDirectApplicInterface::
derived_map_asynch(const Dakota::ParamResponsePair& pair)
{
}

void 
TriKota::MPDirectApplicInterface::
derived_synch(Dakota::PRPQueue& prp_queue)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Stokhos::ProductEpetraVector;
  using Stokhos::ProductEpetraMultiVector;

  if (model != Teuchos::null) {
    unsigned int queue_size = prp_queue.size();
    unsigned int num_blocks = queue_size / block_size;
    unsigned int remainder = queue_size % block_size;

    *out << "TriKota:: Performing multi-point evaluation of " << queue_size
	 << " responses using a block size of " << block_size
	 << std::endl;
    
    Dakota::PRPQueueIter block_iter = prp_queue.begin();
    unsigned int block = 0;
    while (block_iter != prp_queue.end()) {

      // Load up block p
      Dakota::PRPQueueIter iter = block_iter;
      bool valFlag = true;
      bool gradFlag = true;
      bool hessFlag = true;
      unsigned int blk_sz = block_size;
      if (block == num_blocks && remainder > 0)
	blk_sz = remainder;
      for (unsigned int idx=0; idx<blk_sz; idx++) {
	const Dakota::Variables& vars = iter->prp_parameters();
	const Dakota::RealVector& xC  = vars.continuous_variables();
	unsigned int numVars = xC.length();
	for (unsigned int i=0; i<numVars; i++) 
	  (*model_p)[idx][i]=xC[i];

	const Dakota::ActiveSet& set  = iter->active_set();
	short asv = set.request_vector()[0];
	valFlag  = valFlag  && (asv & 1);
	gradFlag = gradFlag && (asv & 2);
	hessFlag = hessFlag && (asv & 4);

	Dakota::Response resp = iter->prp_response();
	Dakota::RealVector fnVals = resp.function_values_view();
	unsigned int numFns = fnVals.length();
	
	// test for consistency of problem definition between ModelEval and 
	// Dakota
	TEUCHOS_TEST_FOR_EXCEPTION(numVars > numParameters, std::logic_error,
				   "TriKota_Dakota Adapter Error: ");
	TEUCHOS_TEST_FOR_EXCEPTION(hessFlag, std::logic_error,
				   "TriKota_Dakota Adapter Error: ");
	TEUCHOS_TEST_FOR_EXCEPTION(numFns > numResponses, std::logic_error,
				   "TriKota_Dakota Adapter Error: ");
	TEUCHOS_TEST_FOR_EXCEPTION(gradFlag && !supportsSensitivities, 
				   std::logic_error,
				   "TriKota_Dakota Adapter Error: ");
	
	++iter;
      }

      // Put in copies of last point for remainder
      if (block == num_blocks && remainder > 0) {
	--iter;
	for (unsigned int idx=remainder; idx<block_size-remainder; idx++) {
	  const Dakota::Variables& vars = iter->prp_parameters();
	  const Dakota::RealVector& xC  = vars.continuous_variables();
	  unsigned int numVars = xC.length();
	  for (unsigned int i=0; i<numVars; i++) 
	    (*model_p)[idx][i]=xC[i];
	}
      }

      // Evaluate model
      EEME::InArgs inArgs = model->createInArgs();
      EEME::OutArgs outArgs = model->createOutArgs();
      inArgs.set_p_mp(p_index, model_p);
      if (valFlag)
	outArgs.set_g_mp(g_index, model_g);
      if (gradFlag) 
	outArgs.set_DgDp_mp(g_index, p_index,
			    EEME::MPDerivative(model_dgdp, orientation));
      *out << "TriKota:: evaluating block " << block;
      if (gradFlag)
	*out << " with gradients." << std::endl;
      else
	*out << " without gradients." << std::endl;
      model->evalModel(inArgs, outArgs);

      // Copy responses from block g
      iter = block_iter;
      for (unsigned int idx=0; idx<blk_sz; idx++) {
	const Dakota::Variables& vars = iter->prp_parameters();
	const Dakota::RealVector& xC  = vars.continuous_variables();
	unsigned int numVars = xC.length();
	Dakota::Response         resp = iter->prp_response(); // shared rep
	Dakota::RealVector fnVals     = resp.function_values_view();
	Dakota::RealMatrix fnGrads    = resp.function_gradients_view();
	unsigned int numFns = fnVals.length();
	
	if (valFlag)
	  for (unsigned int j=0; j<numFns; j++) 
	    fnVals[j]= (*model_g)[idx][j];
	if (gradFlag)  {
	  if (orientation == EEME::DERIV_MV_BY_COL) {
	    for (unsigned int i=0; i<numVars; i++)
	      for (unsigned int j=0; j<numFns; j++)
		fnGrads[j][i]= (*model_dgdp)[idx][i][j];
	  } 
	  else {
	    for (unsigned int i=0; i<numFns; i++)
	      for (unsigned int j=0; j<numVars; j++)
		fnGrads[i][j]= (*model_dgdp)[idx][i][j];
	  }
	}

	// indicate completion of job to ApplicationInterface schedulers
	int fn_eval_id = iter->eval_id();
	completionSet.insert(fn_eval_id);

	++iter;
	++block_iter;
      }
      ++block;
    }
    TEUCHOS_TEST_FOR_EXCEPT((remainder == 0 && block != num_blocks) || 
			    (remainder  > 0 && block != num_blocks+1))
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator() != MPI_COMM_NULL, 
      std::logic_error,
      "\nTriKota Parallelism Error: ModelEvaluator=null, " <<
      "but analysis_comm != MPI_COMMM_NULL");
  }
}
