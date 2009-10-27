#include <iostream>
#include "TriKota_DirectApplicInterface.hpp"
#include "Teuchos_VerboseObject.hpp"

using namespace Dakota;
typedef EpetraExt::ModelEvaluator EEME;

// Define interface class
TriKota::DirectApplicInterface::DirectApplicInterface(
                                ProblemDescDB& problem_db_,
                                const Teuchos::RCP<EpetraExt::ModelEvaluator> App_)
  : Dakota::DirectApplicInterface(problem_db_),
    App(App_),
    orientation(EEME::DERIV_MV_BY_COL)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (App != Teuchos::null) {
    model_p = Teuchos::rcp(new Epetra_Vector(*(App->get_p_init(0))));
    //    model_g = Teuchos::rcp(new Epetra_Vector(*(App->get_g_map(0))));
    model_g = Teuchos::rcp(new Epetra_Vector((const Epetra_BlockMap&) *(App->get_g_map(0)), true));

    numParameters = model_p->GlobalLength();
    numResponses  = model_g->GlobalLength();

    *out << "TriKota:: ModeEval has " << numParameters <<
            " parameters and " << numResponses << " responses." << std::endl;

    EEME::DerivativeSupport supportDgDp = App->createOutArgs().supports(EEME::OUT_ARG_DgDp, 0, 0);
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
        TEST_FOR_EXCEPTION(!supportDgDp.none(), std::logic_error,
              "TriKota Adapter Error: DgDp data type not implemented");
      }
    }

    model_dgdp = Teuchos::rcp(new Epetra_MultiVector(model_p->Map(), model_g->GlobalLength() ));

    *out << "TriKota:: Setting initial guess from Model Evaluator to Dakota " << std::endl;

    Model& first_model = *(problem_db_.model_list().begin());
    int num_dakota_vars =  first_model.acv();
    Dakota::RealVector drv(num_dakota_vars);

    TEST_FOR_EXCEPTION(num_dakota_vars > numParameters, std::logic_error,
                       "TriKota Adapter Error: number of parameters in ModelEvaluator  "
                       <<  numParameters << "\n is less then the number of continuous variables\n"
                       << " specified in the dakota.in input file " << num_dakota_vars << "\n" );

    for (int i=0; i<num_dakota_vars; i++) drv[i] = (*model_p)[i];
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
    TEST_FOR_EXCEPTION(numVars > numParameters, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEST_FOR_EXCEPTION(numADV != 0, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEST_FOR_EXCEPTION(numFns > numResponses, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEST_FOR_EXCEPTION(hessFlag, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

    EEME::InArgs inArgs = App->createInArgs();
    EEME::OutArgs outArgs = App->createOutArgs();

    TEST_FOR_EXCEPTION(gradFlag && !supportsSensitivities, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

//    TEST_FOR_EXCEPTION(parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
//                       != App->MyMPIComm(), std::logic_error,
//     "TriKota Parallelism Error: derived_map_ac called with different MPI_comm");
 

    // Load parameters from Dakota to ModelEval data structure
    for (int i=0; i<numVars; i++) (*model_p)[i]=xC[i];

    // Evaluate model
    inArgs.set_p(0,model_p);
    outArgs.set_g(0,model_g);
    EEME::Derivative model_dgdp_deriv(model_dgdp, orientation);
    if (gradFlag) outArgs.set_DgDp(0,0,model_dgdp_deriv);
    App->evalModel(inArgs, outArgs);

    for (int j=0; j<numFns; j++) fnVals[j]= (*model_g)[j];

    if (gradFlag)  {
      if (orientation == EEME::DERIV_MV_BY_COL) {
        for (int i=0; i<numVars; i++)
          for (int j=0; j<numFns; j++)
            fnGrads[j][i]= (*model_dgdp)[i][j];
      } else {
        for (int i=0; i<numFns; i++)
          for (int j=0; j<numVars; j++)
            fnGrads[i][j]= (*model_dgdp)[i][j];
      }
    }
  }
  else {
    TEST_FOR_EXCEPTION(parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
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
  model_dgdp->Print(*out << "\nSensitivities!\n");

  return 0;
}
