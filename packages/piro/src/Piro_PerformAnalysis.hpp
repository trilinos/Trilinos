// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_PERFORMANALYSIS
#define PIRO_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Piro_ROL_ObserverBase.hpp"

namespace Piro {

  template<class CharT, class Traits=std::char_traits<CharT>>
  class RolOutputBuffer : public std::basic_streambuf<CharT,Traits> {
     public:
     const std::stringstream& getStringStream() const;

     protected:
     inline virtual int overflow(int c = Traits::eof());

     private:
     std::stringstream ss;
  };

  //! \name Top-level Thyra analysis driver
  //@{
  //! \brief Performs analysis of a solved model.
  //! \details This function can either call one of the package-specific drivers
  //!          or perform a \link Piro_Thyra_solve_driver_grp forward solve\endlink.
  //! \ingroup Piro_Thyra_analysis_driver_grp
  int PerformAnalysis(
     Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& analysisParams,
     Teuchos::RCP< Thyra::VectorBase<double> >& result,
     Teuchos::RCP< ROL_ObserverBase<double> > observer = Teuchos::null
     );
  //@}

  //! \brief Performs analysis of a steady state solved model using ROL.
  //! \details Requires that the ROL package is available.
  //! \ingroup Piro_Thyra_analysis_driver_grp
  int PerformROLSteadyAnalysis(
     Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& rolParams,
     Teuchos::RCP< Thyra::VectorBase<double> >& p,
     Teuchos::RCP< ROL_ObserverBase<double> > observer = Teuchos::null
     );

  //! \brief Performs analysis of a transient solved model using ROL.
  //! \details Requires that the ROL package is available.
  //! \ingroup Piro_Thyra_analysis_driver_grp
  int PerformROLTransientAnalysis(
     Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& rolParams,
     Teuchos::RCP< Thyra::VectorBase<double> >& p,
     Teuchos::RCP< ROL_ObserverBase<double> > observer = Teuchos::null
     );

  //! \brief Performs analysis of a solved model using ROL.
  //! \details Requires that the ROL package is available.
  //! \ingroup Piro_Thyra_analysis_driver_grp
  int PerformROLAnalysis(
     Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& rolParams,
     Teuchos::RCP< Thyra::VectorBase<double> >& p,
     Teuchos::RCP< ROL_ObserverBase<double> > observer = Teuchos::null
     );
  //@}

  //! \name Analysis parameter list validation
  //@{
  //! Valid parameters for the list sent to PerformAnalysis
  //! \ingroup Piro_analysis_driver_grp
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidPiroAnalysisParameters();

  //! Valid parameters for the list sent to PerformROLAnalysis
  //! \ingroup Piro_analysis_driver_grp
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidPiroAnalysisROLParameters(int num_parameters);
  //@}
}

#endif //PIRO_PERFORMANALYSIS
