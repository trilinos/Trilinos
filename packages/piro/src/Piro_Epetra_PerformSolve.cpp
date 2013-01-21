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

#include "Piro_Epetra_PerformSolve.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Piro_PerformSolve.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include <iterator>

namespace Piro {

namespace Epetra {

void PerformSolve(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Epetra_Vector> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Epetra_MultiVector> > > &sensitivities)
{
  const Teuchos::RCP<const EpetraExt::ModelEvaluator> epetraSolver = Teuchos::rcpFromRef(piroSolver);
  const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = Teuchos::null;
  const Thyra::EpetraModelEvaluator thyraSolver(epetraSolver, lowsFactory);

  typedef Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<double> > > ThyraResponseArray;
  ThyraResponseArray thyraResponses;

  typedef Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<double> > > > ThyraSensitivityArray;
  ThyraSensitivityArray thyraSensitivities;

  PerformSolveImpl(thyraSolver, solveParams, thyraResponses, thyraSensitivities);

  responses.clear();
  responses.reserve(thyraResponses.size());
  for (ThyraResponseArray::const_iterator it_begin = thyraResponses.begin(),
      it_end = thyraResponses.end(),
      it = it_begin;
      it != it_end;
      ++it) {
    const int j = std::distance(it_begin, it);
    const Epetra_Map g_map = *piroSolver.get_g_map(j);
    const Teuchos::RCP<const Thyra::VectorBase<double> > g_thyra = *it;
    const Teuchos::RCP<const Epetra_Vector> g =
      Teuchos::nonnull(g_thyra) ? Thyra::get_Epetra_Vector(g_map, g_thyra) : Teuchos::null;
    responses.push_back(g);
  }

  sensitivities.clear();
  sensitivities.resize(thyraSensitivities.size());
  for (ThyraSensitivityArray::const_iterator it_begin = thyraSensitivities.begin(),
      it_end = thyraSensitivities.end(),
      it = it_begin;
      it != it_end;
      ++it) {
    const int j = std::distance(it_begin, it);
    const Epetra_Map g_map = *piroSolver.get_g_map(j);
    for (ThyraSensitivityArray::value_type::const_iterator jt = it->begin(), jt_end = it->end(); jt != jt_end; ++jt) {
      const Teuchos::RCP<const Thyra::MultiVectorBase<double> > dgdp_thyra = *jt;
      const Teuchos::RCP<const Epetra_MultiVector> dgdp =
        Teuchos::nonnull(dgdp_thyra) ? Thyra::get_Epetra_MultiVector(g_map, dgdp_thyra) : Teuchos::null;
      sensitivities[j].push_back(dgdp);
    }
  }
}

} // namespace Epetra

} // namespace Piro
