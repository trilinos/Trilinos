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

#ifndef PIRO_PERFORMANALYSIS
#define PIRO_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Piro_ROL_ObserverBase.hpp"

namespace Piro {

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

  //! \brief Performs analysis of a solved model using Dakota via %TriKota.
  //! \details Requires that the %TriKota package is available.
  //! \ingroup Piro_Thyra_analysis_driver_grp
  int PerformDakotaAnalysis(
     Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& dakotaParams,
     Teuchos::RCP< Thyra::VectorBase<double> >& p
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

  //! Valid parameters for the list sent to PerformDakotaAnalysis
  //! \ingroup Piro_analysis_driver_grp
  Teuchos::RCP<const Teuchos::ParameterList>
    getValidPiroAnalysisDakotaParameters();
  //@}
}

#endif //PIRO_PERFORMANALYSIS
