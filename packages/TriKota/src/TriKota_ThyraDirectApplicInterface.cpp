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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include <iostream>
#include "TriKota_ThyraDirectApplicInterface.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"


using namespace Dakota;
typedef Thyra::ModelEvaluatorBase MEB;


// Define interface class
TriKota::ThyraDirectApplicInterface::ThyraDirectApplicInterface(
  ProblemDescDB& problem_db_,
  const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App_,
  int p_index_,
  int g_index_)
  : Dakota::DirectApplicInterface(problem_db_),
    App(App_),
    p_index(p_index_),
    g_index(g_index_),
    orientation(MEB::DERIV_MV_BY_COL)
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  if (App != Teuchos::null) {
    model_p = Thyra::createMember<double>(App->get_p_space(p_index));
    model_g = Thyra::createMember<double>(App->get_g_space(g_index));

    Thyra::DetachedVectorView<double> my_p(model_p);
    Thyra::DetachedVectorView<double> my_g(model_g);
    
    numParameters = my_p.subDim();
    numResponses  = my_g.subDim();

    *out << "TriKota:: ModeEval has " << numParameters <<
            " parameters and " << numResponses << " responses." << std::endl;

    MEB::DerivativeSupport supportDgDp = 
      App->createOutArgs().supports(MEB::OUT_ARG_DgDp, g_index, p_index);
    supportsSensitivities = !(supportDgDp.none());

    // Create the MultiVector, then the Derivative object
    if (supportsSensitivities) {
      *out << "TriKota:: ModeEval supports gradients calculation." << std::endl;

      if (supportDgDp.supports(MEB::DERIV_TRANS_MV_BY_ROW)) {
        orientation = MEB::DERIV_TRANS_MV_BY_ROW;
        model_dgdp = Thyra::createMembers<double>(App->get_p_space(p_index), 
						  numResponses);
      }
      else if (supportDgDp.supports(MEB::DERIV_MV_BY_COL)) {
        orientation = MEB::DERIV_MV_BY_COL;
        model_dgdp = Thyra::createMembers<double>(App->get_g_space(g_index), 
						  numParameters);
      }
      else {
        TEST_FOR_EXCEPTION(!supportDgDp.none(), std::logic_error,
              "TriKota Adapter Error: DgDp data type not implemented");
      }
    }

    *out << "TriKota:: Setting initial guess from Model Evaluator to Dakota " << std::endl;
    const Thyra::ConstDetachedVectorView<double> my_pinit(App->getNominalValues().get_p(p_index));
    for (unsigned int i=0; i<numParameters; i++) my_p[i] = my_pinit[i];

    Model& first_model = *(problem_db_.model_list().begin());
    unsigned int num_dakota_vars =  first_model.acv();
    Dakota::RealVector drv(num_dakota_vars);

    TEST_FOR_EXCEPTION(
      num_dakota_vars > numParameters, std::logic_error,
      "TriKota Adapter Error: number of parameters in ModelEvaluator  " <<  
      numParameters << 
      "\n is less then the number of continuous variables\n" << 
      " specified in the dakota.in input file " << num_dakota_vars << "\n" );

    for (unsigned int i=0; i<num_dakota_vars; i++) drv[i] = my_p[i];
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

    // Test for consistency of problem definition between ModelEval and Dakota
    TEST_FOR_EXCEPTION(numVars > numParameters, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEST_FOR_EXCEPTION(numFns > numResponses, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");
    TEST_FOR_EXCEPTION(hessFlag, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

    MEB::InArgs<double> inArgs = App->createInArgs();
    MEB::OutArgs<double> outArgs = App->createOutArgs();

    TEST_FOR_EXCEPTION(gradFlag && !supportsSensitivities, std::logic_error,
                       "TriKota_Dakota Adapter Error: ");

    // Load parameters from Dakota to ModelEval data structure
    {
      Thyra::DetachedVectorView<double> my_p(model_p);
      for (unsigned int i=0; i<numVars; i++) my_p[i]=xC[i];
    }

    // Evaluate model
    inArgs.set_p(p_index,model_p);
    outArgs.set_g(g_index,model_g);
    if (gradFlag) outArgs.set_DgDp(g_index,p_index,
      MEB::DerivativeMultiVector<double>(model_dgdp,orientation));
    App->evalModel(inArgs, outArgs);

    Thyra::DetachedVectorView<double> my_g(model_g);
    for (unsigned int j=0; j<numFns; j++) fnVals[j]= my_g[j];

    if (gradFlag) {
      if (orientation == MEB::DERIV_MV_BY_COL) {
        for (unsigned int j=0; j<numVars; j++) {
          Thyra::DetachedVectorView<double>
             my_dgdp_j(model_dgdp->col(j));
          for (unsigned int i=0; i<numFns; i++)  fnGrads[i][j]= my_dgdp_j[i];
        }
      }
      else {
        for (unsigned int j=0; j<numFns; j++) {
          Thyra::DetachedVectorView<double>
             my_dgdp_j(model_dgdp->col(j));
          for (unsigned int i=0; i<numVars; i++) fnGrads[j][i]= my_dgdp_j[i]; 
        }
      }
    }
  }
  else {
    TEST_FOR_EXCEPTION(
      parallelLib.parallel_configuration().ea_parallel_level().server_intra_communicator()
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
  return 0;
}
