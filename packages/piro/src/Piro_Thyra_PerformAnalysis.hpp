// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
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

#ifndef PIRO_THYRA_PERFORMANALYSIS
#define PIRO_THYRA_PERFORMANALYSIS

#include "Piro_ConfigDefs.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace Piro {
  namespace Thyra {


// Piro Thyra Piro::Thyra  Deprecated namespace. Too many
//    collisions between  Piro::Thyra and Thyra namespaces.

TEUCHOS_DEPRECATED
  int PerformAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& analysisParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

TEUCHOS_DEPRECATED
  int PerformMoochoAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& moochoParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

TEUCHOS_DEPRECATED
  int PerformDakotaAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& dakotaParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

TEUCHOS_DEPRECATED
  int PerformOptiPackAnalysis(
     ::Thyra::ModelEvaluatorDefaultBase<double>& piroModel,
     Teuchos::ParameterList& optipackParams,
     Teuchos::ParameterList& globipackParams,
     Teuchos::RCP< ::Thyra::VectorBase<double> >& p
     );

TEUCHOS_DEPRECATED
   Teuchos::RCP<const Teuchos::ParameterList>
     getValidPiroAnalysisParameters();

  }
}

#endif //PIRO_THYRA_PERFORMANALYSIS
