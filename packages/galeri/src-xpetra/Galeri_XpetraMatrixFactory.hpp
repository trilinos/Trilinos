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
#ifndef GALERI_XPETRAMATRIXFACTORY_HPP
#define GALERI_XPETRAMATRIXFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Galeri_XpetraMatrixTypes.hpp"

#include <iostream>

namespace Galeri {
  
  namespace Xpetra {

    using Teuchos::RCP;

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    RCP<Matrix>
    CreateCrsMatrix(const std::string &MatrixType, const RCP<const Map> & map, Teuchos::ParameterList& list) //TODO: rename CreateCrsMatrix to CreateMatrix or CreateMatrix ?
    {
      RCP<Matrix> returnMatrix;
      if (MatrixType == "Laplace1D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        if (nx == -1)
          {
            nx = map->getGlobalNumElements();
          }

        returnMatrix = TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, 2.0, -1.0, -1.0);

      } else if (MatrixType == "Laplace2D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        GlobalOrdinal ny = list.get("ny", (GlobalOrdinal) -1);
        if (nx == -1 || ny == -1)
        {
          GlobalOrdinal n = map->getGlobalNumElements();
          nx = (GlobalOrdinal)sqrt((double)n);
          ny = nx;
          TEUCHOS_TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
        }
        bool keepBCs = list.get("keepBCs",false);

        returnMatrix = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, 4.0, -1.0, -1.0, -1.0,
        -1.0,keepBCs);

      } else if (MatrixType == "Star2D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        GlobalOrdinal ny = list.get("ny", (GlobalOrdinal) -1);

        Scalar a = list.get("a", 8.0);
        Scalar b = list.get("b", -1.0);
        Scalar c = list.get("c", -1.0);
        Scalar d = list.get("d", -1.0);
        Scalar e = list.get("e", -1.0);
        Scalar z1 = list.get("z1", -1.0);
        Scalar z2 = list.get("z2", -1.0);
        Scalar z3 = list.get("z3", -1.0);
        Scalar z4 = list.get("z4", -1.0);

        returnMatrix = Star2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, a, b, c, d, e, z1, z2, z3, z4);

      } else if (MatrixType == "BigStar2D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        GlobalOrdinal ny = list.get("ny", (GlobalOrdinal) -1);

        Scalar a = list.get("a", 20.0);
        Scalar b = list.get("b", -8.0);
        Scalar c = list.get("c", -8.0);
        Scalar d = list.get("d", -8.0);
        Scalar e = list.get("e", -8.0);
        Scalar z1 = list.get("z1", 2.0);
        Scalar z2 = list.get("z2", 2.0);
        Scalar z3 = list.get("z3", 2.0);
        Scalar z4 = list.get("z4", 2.0);
        Scalar bb = list.get("bb", 1.0);
        Scalar cc = list.get("cc", 1.0);
        Scalar dd = list.get("dd", 1.0);
        Scalar ee = list.get("ee", 1.0);

        returnMatrix = BigStar2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, a, b, c, d, e, z1, z2, z3, z4, bb, cc, dd, ee);

      } else if (MatrixType == "Laplace3D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        GlobalOrdinal ny = list.get("ny", (GlobalOrdinal) -1);
        GlobalOrdinal nz = list.get("nz", (GlobalOrdinal) -1);
        if (nx == -1 || ny == -1 || nz == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
            ny = nx; nz = nx;
            TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
          } 

        returnMatrix = Cross3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, nz, 6.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

      } else if (MatrixType == "Brick3D") {

        GlobalOrdinal nx = list.get("nx", (GlobalOrdinal) -1);
        GlobalOrdinal ny = list.get("ny", (GlobalOrdinal) -1);
        GlobalOrdinal nz = list.get("nz", (GlobalOrdinal) -1);
        if (nx == -1 || ny == -1 || nz == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
            ny = nx; nz = nx;
            TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
          } 

        returnMatrix = Brick3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

      } else if (MatrixType == "Identity") {

        Scalar a = list.get("a", 1.0);
        returnMatrix = Identity<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, a);

      } else {

        TEUCHOS_TEST_FOR_EXCEPTION(true,
                           std::logic_error,
                           "`MatrixType' has incorrect value (" << MatrixType
                           << ") in input to function CreateCrsMatrix()."
                           << "Check the documentation for a list of valid choices");
      } //if-else

      returnMatrix->setObjectLabel(MatrixType);
      return returnMatrix;

    } // CreateCrsMatrix()

  } // namespace Xpetra
} // namespace Galeri

#endif //ifndef GALERI_XPETRAMATRIXFACTORY_HPP
