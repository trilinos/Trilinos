// @HEADER
// ************************************************************************
// 
//        TriKota: A Trilinos Wrapper for the Dakota Framework
//                  Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <iostream>
#include "TriKota_DirectApplicInterface.hpp"
#include "DakotaModel.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace Dakota;
typedef EpetraExt::ModelEvaluator EEME;

// Define interface class
TriKota::DirectApplicInterface::DirectApplicInterface(
                                ProblemDescDB& problem_db_,
                                const Teuchos::RCP<EpetraExt::ModelEvaluator> App_,
				int p_index_, int g_index_)
  : Dakota::DirectApplicInterface(problem_db_),
    App(App_),
    p_index(p_index_),
    g_index(g_index_),
    orientation(EEME::DERIV_MV_BY_COL)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (App != Teuchos::null) {
    model_p = Teuchos::rcp(new Epetra_Vector(*(App->get_p_init(p_index))));
    //    model_g = Teuchos::rcp(new Epetra_Vector(*(App->get_g_map(g_index))));
    model_g = Teuchos::rcp(new Epetra_Vector((const Epetra_BlockMap&) *(App->get_g_map(g_index)), true));

    numParameters = model_p->GlobalLength();
    numResponses  = model_g->GlobalLength();

    *out << "TriKota:: ModeEval has " << numParameters <<
            " parameters and " << numResponses << " responses." << std::endl;

    EEME::DerivativeSupport supportDgDp = App->createOutArgs().supports(EEME::OUT_ARG_DgDp, g_index, p_index);
    supportsSensitivities = !(supportDgDp.none());

    // Create the MultiVector, then the Derivative object
    if (supportsSensitivities) {
      *out << "TriKota:: ModeEval supports gradients calculation." << std::endl;

      if (supportDgDp.supports(EEME::DERIV_TRANS_MV_BY_ROW)) {
        orientation = EEME::DERIV_TRANS_MV_BY_ROW;
        model_dgdp = Teuchos::rcp(new Epetra_MultiVector(model_p->Map(), numResponses ));
      }
      else if (supportDgDp.supports(EEME::DERIV_MV_BY_COL)) {
        orientation = EEME::DERIV_MV_BY_COL;
        model_dgdp = Teuchos::rcp(new Epetra_MultiVector(model_g->Map(), numParameters));
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(!supportDgDp.none(), std::logic_error,
              "TriKota Adapter Error: DgDp data type not implemented");
      }
    }

    *out << "TriKota:: Setting initial guess from Model Evaluator to Dakota " << std::endl;

    Model& first_model = *(problem_db_.model_list().begin());
    unsigned int num_dakota_vars =  first_model.acv();
    Dakota::RealVector drv(num_dakota_vars);

    TEUCHOS_TEST_FOR_EXCEPTION(num_dakota_vars > numParameters, std::logic_error,
                       "TriKota Adapter Error: number of parameters in ModelEvaluator  "
                       <<  numParameters << "\n is less then the number of continuous variables\n"
                       << " specified in the dakota.in input file " << num_dakota_vars << "\n" );

    for (unsigned int i=0; i<num_dakota_vars; i++) drv[i] = (*model_p)[i];
    first_model.continuous_variables(drv);

  }
  else {
    *out << "Warning in TriKota::DirectApplicInterface constructor\n" 
         << "\tModelEvaluator is null. This is OK iff Dakota has assigned"
         << " MPI_COMM_NULL to this Proc " << std::endl;
  }
}

int TriKota::DirectApplicInterface::derived_map_ac(const Dakota::String& ac_name)
{

  if (App != Teuchos::null) {

    // test for consistency of problem definition between ModelEval and Dakota
    TEUCHOS_TEST_FOR_EXCEPTION(numVars > numParameters, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEUCHOS_TEST_FOR_EXCEPTION(numFns > numResponses, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEUCHOS_TEST_FOR_EXCEPTION(hessFlag, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

    EEME::InArgs inArgs = App->createInArgs();
    EEME::OutArgs outArgs = App->createOutArgs();

    TEUCHOS_TEST_FOR_EXCEPTION(gradFlag && !supportsSensitivities, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

//    TEUCHOS_TEST_FOR_EXCEPTION(parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
//                       != App->MyMPIComm(), std::logic_error,
//     "TriKota Parallelism Error: derived_map_ac called with different MPI_comm");
 

    // Load parameters from Dakota to ModelEval data structure
    for (unsigned int i=0; i<numVars; i++) (*model_p)[i]=xC[i];

    // Evaluate model
    inArgs.set_p(p_index,model_p);
    outArgs.set_g(g_index,model_g);
    EEME::Derivative model_dgdp_deriv(model_dgdp, orientation); //may be null
    if (gradFlag) outArgs.set_DgDp(g_index,p_index,model_dgdp_deriv);
    App->evalModel(inArgs, outArgs);

    for (unsigned int j=0; j<numFns; j++) fnVals[j]= (*model_g)[j];

    if (gradFlag)  {
      if (orientation == EEME::DERIV_MV_BY_COL) {
        for (unsigned int i=0; i<numVars; i++)
          for (unsigned int j=0; j<numFns; j++)
            fnGrads[j][i]= (*model_dgdp)[i][j];
      } else {
        for (unsigned int i=0; i<numFns; i++)
          for (unsigned int j=0; j<numVars; j++)
            fnGrads[i][j]= (*model_dgdp)[i][j];
      }
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
                       != MPI_COMM_NULL, std::logic_error,
                      "\nTriKota Parallelism Error: ModelEvaluator=null, but analysis_comm != MPI_COMMM_NULL");
  }

  return 0;
}

int TriKota::DirectApplicInterface::derived_map_of(const Dakota::String& ac_name)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // WARNING:: For some reason, this method is not being called!
  *out << "Finished Dakota NLS Fitting!: " << std::setprecision(5) << std::endl;
  model_p->Print(*out << "\nParameters!\n");
  model_g->Print(*out << "\nResponses!\n");
  if (gradFlag)
    model_dgdp->Print(*out << "\nSensitivities!\n");

  return 0;
}
