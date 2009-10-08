#include <iostream>
#include "TriKota_ThyraDirectApplicInterface.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"


using namespace Dakota;

// Define interface class
TriKota::ThyraDirectApplicInterface::ThyraDirectApplicInterface(
              ProblemDescDB& problem_db_,
              const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App_)
              //const Teuchos::RCP<Thyra::ModelEvaluator<double> > App_)
  : Dakota::DirectApplicInterface(problem_db_),
    App(App_)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (App != Teuchos::null) {
    model_p = Thyra::createMember<double>(App->get_p_space(0));
    model_g = Thyra::createMember<double>(App->get_g_space(0));

    Thyra::DetachedSpmdVectorView<double> my_p(model_p);
    Thyra::DetachedSpmdVectorView<double> my_g(model_g);
    
    numParameters = my_p.subDim();
    numResponses  = my_g.subDim();

    // Create the MultiVector
    model_dgdp = Thyra::createMembers<double>(App->get_p_space(0), numResponses);

    *out << "TriKota:: ModeEval has " << numParameters <<
            " parameters and " << numResponses << " responses." << std::endl;

    *out << "TriKota:: Setting initial guess from Model Evaluator to Dakota " << std::endl;
    const Thyra::ConstDetachedSpmdVectorView<double> my_pinit(App->getNominalValues().get_p(0));
    for (int i=0; i<numParameters; i++) my_p[i] = my_pinit[i];

    Model& first_model = *(problem_db_.model_list().begin());
    int num_dakota_vars =  first_model.acv();
    Dakota::RealVector drv(num_dakota_vars);

    TEST_FOR_EXCEPTION(num_dakota_vars > numParameters, std::logic_error,
                       "TriKota Adapter Error: number of parameters in ModelEvaluator  "
                       <<  numParameters << "\n is less then the number of continuous variables\n"
                       << " specified in the dakota.in input file " << num_dakota_vars << "\n" );

    for (int i=0; i<num_dakota_vars; i++) drv[i] = my_p[i];
    first_model.continuous_variables(drv);

  }
  else {
    *out << "Warning in TriKota::ThyraDirectApplicInterface constructor\n" 
         << "\tModelEvaluator is null. This is OK iff Dakota has assigned"
         << " MPI_COMM_NULL to this Proc " << std::endl;
  }
}

int TriKota::ThyraDirectApplicInterface::derived_map_ac(const Dakota::String& ac_name)
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

    Thyra::ModelEvaluatorBase::InArgs<double> inArgs = App->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<double> outArgs = App->createOutArgs();

    bool supportsSensitivities = 
      !outArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp, 0, 0).none();
    //   false;

    TEST_FOR_EXCEPTION(gradFlag && !supportsSensitivities, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

//    TEST_FOR_EXCEPTION(parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
//                       != App->MyMPIComm(), std::logic_error,
//     "TriKota Parallelism Error: derived_map_ac called with different MPI_comm");
 

    // Load parameters from Dakota to ModelEval data structure
    Thyra::DetachedSpmdVectorView<double> my_p(model_p);
    Thyra::DetachedSpmdVectorView<double> my_g(model_g);
    for (int i=0; i<numVars; i++) my_p[i]=xC[i];

    if (gradFlag) {
      for (int j=0; j<numFns; j++) {
        Thyra::DetachedSpmdVectorView<double> my_dgdp_j(model_dgdp->col(j));
        for (int i=0; i<numVars; i++) {
           my_dgdp_j[i] = 0.0;
        }
      }
    }

    // Evaluate model
    inArgs.set_p(0,model_p);
    outArgs.set_g(0,model_g);
    if (gradFlag) outArgs.set_DgDp(0,0,model_dgdp);
    App->evalModel(inArgs, outArgs);

    for (int j=0; j<numFns; j++) fnVals[j]= my_g[j];

    if (gradFlag) {
      for (int j=0; j<numFns; j++) {
        Thyra::DetachedSpmdVectorView<double> my_dgdp_j(model_dgdp->col(j));
        for (int i=0; i<numVars; i++) {
          fnGrads[j][i]= my_dgdp_j[i];
        }
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

int TriKota::ThyraDirectApplicInterface::derived_map_of(const Dakota::String& ac_name)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // WARNING:: For some reason, this method is not being called!
  *out << "Finished Dakota NLS Fitting!: " << std::setprecision(5) << std::endl;
/*
  model_p->Print(*out << "\nParameters!\n");
  model_g->Print(*out << "\nResponses!\n");
  model_dgdp->Print(*out << "\nSensitivities!\n");
*/

  return 0;
}
