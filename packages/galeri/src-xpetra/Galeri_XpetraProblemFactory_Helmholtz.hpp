// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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
    RCP<Problem<Map,Matrix,MultiVector> > BuildProblem(const std::string &MatrixType, const RCP<const Map>& map, Teuchos::ParameterList& list) {
      RCP<Problem<Map,Matrix,MultiVector> > P;

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
