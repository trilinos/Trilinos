/*
  Direct translation of parts of Galeri matrix generator.
*/
#ifndef MUELU_MATRIXFACTORY_HPP
#define MUELU_MATRIXFACTORY_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

#include "MueLu_MatrixTypes.hpp"

#include <iostream>

namespace MueLu {
  
  namespace Gallery {

    using Teuchos::RCP;

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename Matrix>
    RCP<Matrix>
    CreateCrsMatrix(const std::string &MatrixType, const RCP<const Map> & map, Teuchos::ParameterList& list) //TODO: rename CreateCrsMatrix to CreateMatrix or CreateOperator ?
    {
      RCP<Matrix> returnMatrix;
      if (MatrixType == "Laplace1D") {

        GlobalOrdinal nx = list.get("nx", -1);
        if (nx == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal)sqrt((Scalar)n);
            TEST_FOR_EXCEPTION(nx*nx != n, std::logic_error, "You need to specify nx.");
          }

        returnMatrix = TriDiag<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, 2.0, -1.0, -1.0);

      } else if (MatrixType == "Laplace2D") {

        GlobalOrdinal nx = list.get("nx", -1);
        GlobalOrdinal ny = list.get("ny", -1);
        if (nx == -1 || ny == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal)sqrt((Scalar)n);
            ny = nx;
            TEST_FOR_EXCEPTION(nx*ny != n, std::logic_error, "You need to specify nx and ny.");
          }

        returnMatrix = Cross2D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, 4.0, -1.0, -1.0, -1.0, -1.0);

      } else if (MatrixType == "Star2D") {

        GlobalOrdinal nx = list.get("nx", -1);
        GlobalOrdinal ny = list.get("ny", -1);

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

        GlobalOrdinal nx = list.get("nx", -1);
        GlobalOrdinal ny = list.get("ny", -1);

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

        GlobalOrdinal nx = list.get("nx", -1);
        GlobalOrdinal ny = list.get("ny", -1);
        GlobalOrdinal nz = list.get("nz", -1);
        if (nx == -1 || ny == -1 || nz == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal) Teuchos::ScalarTraits<Scalar>::pow(n, 0.33334);
            ny = nx; nz = nx;
            TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
          } 

        returnMatrix = Cross3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, nz, 6.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

      } else if (MatrixType == "Brick3D") {

        GlobalOrdinal nx = list.get("nx", -1);
        GlobalOrdinal ny = list.get("ny", -1);
        GlobalOrdinal nz = list.get("nz", -1);
        if (nx == -1 || ny == -1 || nz == -1)
          {
            GlobalOrdinal n = map->getGlobalNumElements();
            nx = (GlobalOrdinal) Teuchos::ScalarTraits<Scalar>::pow(n, 0.33334);
            ny = nx; nz = nx;
            TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
          } 

        returnMatrix = Brick3D<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, nx, ny, nz, 26.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0);

      } else if (MatrixType == "Identity") {

        Scalar a = list.get("a", 1.0);
        returnMatrix = Identity<Scalar,LocalOrdinal,GlobalOrdinal,Map,Matrix>(map, a);

      } else {

        TEST_FOR_EXCEPTION(true,
                           std::logic_error,
                           "`MatrixType' has incorrect value (" << MatrixType
                           << ") in input to function CreateCrsMatrix()."
                           << "Check the documentation for a list of valid choices");
      } //if-else

      returnMatrix->setObjectLabel(MatrixType);
      return returnMatrix;

    } // CreateCrsMatrix()

  } // namespace Gallery
} // namespace MueLu

#define MUELU_MATRIXFACTORY_SHORT

#endif //ifndef MUELU_MATRIXFACTORY_HPP
