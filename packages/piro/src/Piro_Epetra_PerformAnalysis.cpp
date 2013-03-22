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

#include "Piro_Epetra_PerformAnalysis.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Piro_PerformAnalysis.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"

int Piro::Epetra::PerformAnalysis(
    EpetraExt::ModelEvaluator &piroModel,
    Teuchos::ParameterList &analysisParams,
    Teuchos::RCP<Epetra_Vector> &p)
{
  const Teuchos::RCP<const EpetraExt::ModelEvaluator> epetraModel = Teuchos::rcpFromRef(piroModel);
  const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = Teuchos::null;
  Thyra::EpetraModelEvaluator thyraModel(epetraModel, lowsFactory);

  Teuchos::RCP<Thyra::VectorBase<double> > p_thyra;
  const int status = ::Piro::PerformAnalysis(thyraModel, analysisParams, p_thyra);

  if (Teuchos::nonnull(p_thyra)) {
    const int l = 0;
    const Teuchos::RCP<const Epetra_Map> p_map = epetraModel->get_p_map(l);
    p = Thyra::get_Epetra_Vector(*p_map, p_thyra);
  } else {
    p = Teuchos::null;
  }

  return status;
}
