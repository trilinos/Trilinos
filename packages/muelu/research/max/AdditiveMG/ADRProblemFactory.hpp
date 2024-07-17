// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Direct translation of parts of Galeri matrix generator.
*/
#ifndef ADRPROBLEMFACTORY_HPP
#define ADRPROBLEMXFACTORY_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Assert.hpp>

#include "ADR_Problem.hpp"
#include "CreateADRMatrix.hpp"

#include <iostream>

namespace ADR {

namespace Xpetra {

using Teuchos::RCP;

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
RCP<Problem<Map, Matrix, MultiVector> > BuildProblem(const std::string& MatrixType, const RCP<const Map>& map, Teuchos::ParameterList& list) {
  RCP<Problem<Map, Matrix, MultiVector> > P;

  if (MatrixType == "ADR1D")
    P.reset(new ADR1DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>(list, map));
  else if (MatrixType == "ADR2D")
    P.reset(new ADR2DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>(list, map));
  else if (MatrixType == "ADR3D")
    P.reset(new ADR3DProblem<Scalar, LocalOrdinal, GlobalOrdinal, Map, Matrix, MultiVector>(list, map));

  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                               "`MatrixType' has incorrect value (" << MatrixType << ") in input to function CreateCrsMatrix()."
                                                                    << "Check the documentation for a list of valid choices");

  P->setObjectLabel(MatrixType);

  return P;
}

}  // namespace Xpetra

}  // namespace ADR

#endif  // ifndef ADRPROBLEMFACTORY_HPP
