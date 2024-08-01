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
#ifndef GALERI_XPETRAMATRIXFACTORY_HPP
#define GALERI_XPETRAMATRIXFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Galeri_StencilProblems.hpp"
#include "Galeri_Elasticity2DProblem.hpp"
#include "Galeri_Elasticity3DProblem.hpp"

#include <iostream>

namespace Galeri {

  namespace Xpetra {

    using Teuchos::RCP;

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix, typename MultiVector>
    RCP<Problem<Map,Matrix,MultiVector> > BuildProblem(const std::string &MatrixType, const RCP<const Map>& map, Teuchos::ParameterList& list) {
      RCP<Problem<Map,Matrix,MultiVector> > P;

      if      (MatrixType == "Laplace1D")    P.reset(new Laplace1DProblem   <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Laplace2D")    P.reset(new Laplace2DProblem   <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Star2D")       P.reset(new Star2DProblem      <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "BigStar2D")    P.reset(new BigStar2DProblem   <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "AnisotropicDiffusion")    P.reset(new AnisotropicDiffusion2DProblem   <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Laplace3D")    P.reset(new Laplace3DProblem   <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Brick3D")      P.reset(new Brick3DProblem     <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Elasticity2D") P.reset(new Elasticity2DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Elasticity3D") P.reset(new Elasticity3DProblem<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else if (MatrixType == "Identity")     P.reset(new IdentityProblem    <Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix,MultiVector>(list, map));
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                   "`MatrixType' has incorrect value (" << MatrixType << ") in input to function CreateCrsMatrix()."
                                   << "Check the documentation for a list of valid choices");

      P->setObjectLabel(MatrixType);

      return P;
    }

  } // namespace Xpetra

} // namespace Galeri

#endif //ifndef GALERI_XPETRAMATRIXFACTORY_HPP
