// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_ADAPTERS_STK_MAIN_DRIVER_USER_APP_UTILITIES_HPP
#define PANZER_ADAPTERS_STK_MAIN_DRIVER_USER_APP_UTILITIES_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

#include "Panzer_NodeType.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_ModelEvaluatorFactory.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "NOX_Utils.H"
#include "NOX_Observer_Print.hpp"

#include "Tempus_IntegratorBasic.hpp"

#include "user_app_ClosureModel_Factory_TemplateBuilder.hpp"
#include "user_app_EquationSetFactory.hpp"
#include "user_app_BCStrategy_Factory.hpp"
#include "user_app_NOXObserverFactory.hpp"
#include "user_app_TempusObserverFactory.hpp"
#include "user_app_ResponseEvaluatorFactory_HOFlux.hpp"

#include <iostream>
#include <tuple>

namespace user_app {

  void addResponsesToModelEvaluatorFactory(const Teuchos::ParameterList& response_sublist,
                                           panzer_stk::ModelEvaluatorFactory<double>& me_factory);

  // void computeAndPrintResponses(const Teuchos::RCP<Thyra::ModelEvaluator<double>>& me,
  //                               const auto& x,
  //                               const auto& x_dot,
  //                               const auto& t,
  //                               const auto& global_data,
  //                               std::ostream& out);

  std::tuple< Teuchos::RCP<Thyra::ModelEvaluator<double>>,
              Teuchos::RCP<panzer::GlobalData>,
              Teuchos::RCP<panzer_stk::STK_Interface>,
              Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>, // normal response lib
              Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>>, // stk io response lib
              Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>>,
              Teuchos::RCP<panzer::GlobalIndexer>
              >
  buildModelEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& input_params,
                      const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

  Teuchos::RCP<Tempus::IntegratorBasic<double>>
  buildTimeIntegrator(const Teuchos::RCP<Teuchos::ParameterList>& input_params,
                      const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                      Teuchos::RCP<Thyra::ModelEvaluator<double>> me,
                      Teuchos::RCP<panzer_stk::STK_Interface> mesh,
                      Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> rLibrary,
                      Teuchos::RCP<panzer::ResponseLibrary<panzer::Traits>> stkIOResponseLibrary,
                      Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits>> linObjFactory,
                      Teuchos::RCP<panzer::GlobalIndexer> globalIndexer,
                      const bool overrideNoxOutput);

  std::tuple<int,int> findParameterIndex(const std::string& p_name,const Thyra::ModelEvaluator<double>& me);

  std::tuple<int,int> findResponseIndex(const std::string& g_name,const Thyra::ModelEvaluator<double>& me);
}

/*
void user_app::computeAndPrintResponses(const Teuchos::RCP<Thyra::ModelEvaluator<double>>& physics,
                                        const auto& x,
                                        const auto& x_dot,
                                        const auto& t,
                                        const auto& global_data,
                                        std::ostream& os)
{
  if(physics->Ng()>0) {

    os << "Ng = " << physics->Ng() << std::endl;
    os << "Np = " << physics->Np() << std::endl;
    Thyra::ModelEvaluatorBase::InArgs<double> respInArgs = physics->createInArgs();
    Thyra::ModelEvaluatorBase::OutArgs<double> respOutArgs = physics->createOutArgs();

    TEUCHOS_ASSERT(physics->Ng()==respOutArgs.Ng());

    respInArgs.set_x(x);
    respInArgs.set_x_dot(x_dot);

    for(int g=0;g<respOutArgs.Ng();++g) {
      Teuchos::RCP<Thyra::VectorBase<double> > response = Thyra::createMember(*physics->get_g_space(g));
      respOutArgs.set_g(g,response);

      for (int p = 0; p < physics->Np(); ++p) {
        bool derivative_supported = respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,g,p).supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM);
        os << "DgDp(" << g << "," << p << ") supports = " << derivative_supported << std::endl;
        if (derivative_supported) {
          auto dgdp = physics->create_DgDp_op(g,p);
          respOutArgs.set_DgDp(g,p,Thyra::ModelEvaluatorBase::Derivative<double>(dgdp));
        }
      }
    }

    physics->evalModel(respInArgs, respOutArgs);

    {
      auto parameter_library = global_data->pl;
      TEUCHOS_ASSERT(nonnull(parameter_library));
      for (auto i=parameter_library->begin();i != parameter_library->end(); ++i) {
        os << "Sacado::ParameterLibrary: " << i->first << std::endl;
      }
    }

    {
      auto nominalValues = physics->getNominalValues();
      for (int i=0; i < respInArgs.Np();++i) {
        auto p = nominalValues.get_p(i);
        auto p_names = physics->get_p_names(i);
        os << "ModelEvaluator::NominalValues Parameter Value: \"" << (*p_names)[0] << "\" = " << Thyra::get_ele(*p,0) << std::endl;
      }
    }

    for (int i=0;i<respOutArgs.Ng();i++) {
      Teuchos::RCP<Thyra::VectorBase<double> > response = respOutArgs.get_g(i);
      TEUCHOS_ASSERT(response!=Teuchos::null);
      os << "Response Value: \"" << physics->get_g_names(i)[0] << "\" = " << Thyra::get_ele(*response,0) << std::endl;
      for (int j=0; j < respOutArgs.Np(); ++j) {
        // os << "  dg(" << i << ")/dp(" << j << ") supports(GRAD_FORM) = " << respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,j).supports(Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM) << std::endl;
        // os << "  dg(" << i << ")/dp(" << j << ") supports(JAC_FORM)  = " << respOutArgs.supports(Thyra::ModelEvaluatorBase::OUT_ARG_DgDp,i,j).supports(Thyra::ModelEvaluatorBase::DERIV_MV_JACOBIAN_FORM) << std::endl;
      }
    }
  }
}
*/

#endif
