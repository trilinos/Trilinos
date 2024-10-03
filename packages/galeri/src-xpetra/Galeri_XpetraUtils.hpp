// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
/*
  Direct translation of Galeri coordinate generator.
*/
#ifndef GALERI_XPETRAUTILS_HPP
#define GALERI_XPETRAUTILS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"

#include "Galeri_VectorTraits.hpp"
#include "Galeri_Exception.h"

#include <iostream>
namespace Galeri {
  namespace Xpetra {

    template<typename T>
    std::string toString(T &a) {
      std::ostringstream ss;
      ss << a;
      return ss.str();
    }
  }
}

namespace Galeri {
  namespace Xpetra {
  class Utils {

    public:

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename RealValuedMultiVector>
    static Teuchos::RCP<RealValuedMultiVector>
    CreateCartesianCoordinates(std::string const &coordType, RCP<const Map> const & map, Teuchos::ParameterList& list) {
      using Galeri::Xpetra::VectorTraits;
      typedef typename RealValuedMultiVector::scalar_type real_type;

      Teuchos::RCP<RealValuedMultiVector> coordinates;

      GlobalOrdinal ix, iy, iz;
      real_type delta_x, delta_y, delta_z;

      real_type lx = Teuchos::as<real_type>(list.get<double>("lx", 1) * list.get<double>("stretchx", 1));
      real_type ly = Teuchos::as<real_type>(list.get<double>("ly", 1) * list.get<double>("stretchy", 1));
      real_type lz = Teuchos::as<real_type>(list.get<double>("lz", 1) * list.get<double>("stretchz", 1));

      GlobalOrdinal nx = -1;
      GlobalOrdinal ny = -1;
      GlobalOrdinal nz = -1;

      if (list.isParameter("nx")) {
        if (list.isType<int>("nx"))
          nx = Teuchos::as<GlobalOrdinal>(list.get<int>("nx"));
        else
          nx = list.get<GlobalOrdinal>("nx");
      }
      if (list.isParameter("ny")) {
        if (list.isType<int>("ny"))
          ny = Teuchos::as<GlobalOrdinal>(list.get<int>("ny"));
        else
          ny = list.get<GlobalOrdinal>("ny");
      }
      if (list.isParameter("nz")) {
        if (list.isType<int>("nz"))
          nz = Teuchos::as<GlobalOrdinal>(list.get<int>("nz"));
        else
          nz = list.get<GlobalOrdinal>("nz");
      }

      size_t NumMyElements = map->getLocalNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getLocalElementList();

      if (coordType == "1D") {
        coordinates = VectorTraits<Map,RealValuedMultiVector>::Build(map, 1, false);
        Teuchos::ArrayRCP<Teuchos::ArrayRCP<real_type> > Coord(1);
        Coord[0] = coordinates->getDataNonConst(0);

        delta_x = lx / Teuchos::as<real_type>(nx - 1);

        for (size_t i = 0; i < NumMyElements; ++i) {
          ix = MyGlobalElements[i];
          Coord[0][i] = delta_x * Teuchos::as<real_type>(ix);
        }

      } else if (coordType == "2D") {
        coordinates = VectorTraits<Map,RealValuedMultiVector>::Build(map, 2, false);
        Teuchos::ArrayRCP<Teuchos::ArrayRCP<real_type> > Coord(2);
        Coord[0] = coordinates->getDataNonConst(0);
        Coord[1] = coordinates->getDataNonConst(1);

        delta_x = lx / Teuchos::as<real_type>(nx - 1);
        delta_y = ly / Teuchos::as<real_type>(ny - 1);

        for (size_t i = 0; i < NumMyElements; ++i) {
          ix = MyGlobalElements[i] % nx;
          iy = (MyGlobalElements[i] - ix) / nx;

          Coord[0][i] = delta_x * Teuchos::as<real_type>(ix);
          Coord[1][i] = delta_y * Teuchos::as<real_type>(iy);
        }

      } else if (coordType == "3D") {
        coordinates = VectorTraits<Map,RealValuedMultiVector>::Build(map, 3, false);
        Teuchos::ArrayRCP<Teuchos::ArrayRCP<real_type> > Coord(3);
        Coord[0] = coordinates->getDataNonConst(0);
        Coord[1] = coordinates->getDataNonConst(1);
        Coord[2] = coordinates->getDataNonConst(2);

        delta_x = lx / Teuchos::as<real_type>(nx - 1);
        delta_y = ly / Teuchos::as<real_type>(ny - 1);
        delta_z = lz / Teuchos::as<real_type>(nz - 1);

        for (size_t i = 0; i < NumMyElements; i++) {
          GlobalOrdinal ixy = MyGlobalElements[i] % (nx * ny);
          iz = (MyGlobalElements[i] - ixy) / (nx * ny);

          ix = ixy % nx;
          iy = (ixy - ix) / nx;

          Coord[0][i] = delta_x * Teuchos::as<real_type>(ix);
          Coord[1][i] = delta_y * Teuchos::as<real_type>(iy);
          Coord[2][i] = delta_z * Teuchos::as<real_type>(iz);
        }

      } else {
        throw(std::runtime_error("in Galeri::Xpetra::Utils : `coordType' has incorrect value (" + coordType + ")"));
      } //if (coordType == ...

      return coordinates;

    } // CreateCartesianCoordinates()

    template <typename GlobalOrdinal>
    static void getSubdomainData(GlobalOrdinal n, GlobalOrdinal m, GlobalOrdinal i, GlobalOrdinal& start, GlobalOrdinal& end) {
      TEUCHOS_TEST_FOR_EXCEPTION(i>=m, std::runtime_error, "Incorrect input parameter to getSubdomainData: m = " + toString(m) + ", i = " + toString(i));

      // If the number of points is not multiple of the number of subdomains, we assign
      // extra points to the first few subdomains
      GlobalOrdinal minWidth = Teuchos::as<GlobalOrdinal>(floor(Teuchos::as<double>(n)/m));
      GlobalOrdinal numWides = n - m * minWidth;
      if (i < numWides) {
        start = i * (minWidth + 1);
        end = start + minWidth + 1;
      } else {
        start = numWides * (minWidth + 1) + (i - numWides) * minWidth;
        end = start + minWidth;
      }
    }

  }; // class Utils
  } // namespace Xpetra
} // namespace Galeri

#endif //ifndef GALERI_XPETRAUTILS_HPP
