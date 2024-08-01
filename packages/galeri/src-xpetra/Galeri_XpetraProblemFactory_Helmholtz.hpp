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
#ifndef GALERI_XPETRAMATRIXFACTORY_HELMHOLTZ_HPP
#define GALERI_XPETRAMATRIXFACTORY_HELMHOLTZ_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Galeri_StencilProblems_Helmholtz.hpp"
#include "Galeri_HelmholtzFEM2DProblem.hpp"
#include "Galeri_HelmholtzFEM3DProblem.hpp"

#include <iostream>

namespace Galeri {

  namespace Xpetra {

    using Teuchos::RCP;

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<Problem_Helmholtz<Map,Matrix,MultiVector> > BuildProblem_Helmholtz(const std::string &MatrixType, const RCP<const Map>& map, Teuchos::ParameterList& list) {
      RCP<Problem_Helmholtz<Map,Matrix,MultiVector> > P;

      if (MatrixType == "Helmholtz1D")
        P.reset(new Helmholtz1DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));

      else if (MatrixType == "Helmholtz2D")
        P.reset(new Helmholtz2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));

      else if (MatrixType == "Helmholtz3D")
        P.reset(new Helmholtz3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));

      else if (MatrixType == "HelmholtzFEM2D")
        P.reset(new HelmholtzFEM2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));

      else if (MatrixType == "HelmholtzFEM3D")
        P.reset(new HelmholtzFEM3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));

      else
        TEUCHOS_TEST_FOR_EXCEPTION(true,
                                   std::logic_error,
                                   "`MatrixType' has incorrect value (" << MatrixType << ") in input to function CreateCrsMatrix()."
                                   << "Check the documentation for a list of valid choices");

      P->setObjectLabel(MatrixType);
      return P;

    } // BuildProblem()

  } // namespace Xpetra
} // namespace Galeri

#endif //ifndef GALERI_XPETRAMATRIXFACTORY_HELMHOLTZ_HPP
