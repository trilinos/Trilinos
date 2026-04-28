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
#include "Tpetra_Access.hpp"

#include <iostream>
namespace Galeri {
namespace Xpetra {

template <typename T>
std::string toString(T& a) {
  std::ostringstream ss;
  ss << a;
  return ss.str();
}
}  // namespace Xpetra
}  // namespace Galeri

namespace Galeri {
namespace Xpetra {
class Utils {
 public:
  template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map, typename RealValuedMultiVector>
  static Teuchos::RCP<RealValuedMultiVector>
  CreateCartesianCoordinates(std::string const& coordType, RCP<const Map> const& map, Teuchos::ParameterList& list) {
    using Galeri::Xpetra::VectorTraits;
    typedef typename RealValuedMultiVector::scalar_type real_type;

    Teuchos::RCP<RealValuedMultiVector> coordinates;

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

    using Node       = typename Map::node_type;
    using exec_space = typename Node::execution_space;
    using range_type = Kokkos::RangePolicy<LocalOrdinal, exec_space>;
    {
      int dim = 0;
      if (coordType == "1D")
        dim = 1;
      else if (coordType == "2D")
        dim = 2;
      else if (coordType == "3D")
        dim = 3;
      coordinates = VectorTraits<Map, RealValuedMultiVector>::Build(map, dim, false);
      typename RealValuedMultiVector::dual_view_type::t_dev_um lclCoords;
      lclCoords   = coordinates->getLocalViewDevice(::Tpetra::Access::OverwriteAll);
      auto lclMap = map->getLocalMap();
      if (coordType == "1D") {
        delta_x = lx / Teuchos::as<real_type>(nx - 1);
        Kokkos::parallel_for(
            "coords1d",
            range_type(0, NumMyElements),
            KOKKOS_LAMBDA(const LocalOrdinal i) {
              auto ix         = lclMap.getGlobalElement(i);
              lclCoords(i, 0) = delta_x * (real_type)ix;
            });
      } else if (coordType == "2D") {
        delta_x = lx / Teuchos::as<real_type>(nx - 1);
        delta_y = ly / Teuchos::as<real_type>(ny - 1);
        Kokkos::parallel_for(
            "coords2d",
            range_type(0, NumMyElements),
            KOKKOS_LAMBDA(const LocalOrdinal i) {
              auto ix = lclMap.getGlobalElement(i) % nx;
              auto iy = (lclMap.getGlobalElement(i) - ix) / nx;

              lclCoords(i, 0) = delta_x * (real_type)ix;
              lclCoords(i, 1) = delta_y * (real_type)iy;
            });
      } else if (coordType == "3D") {
        delta_x = lx / Teuchos::as<real_type>(nx - 1);
        delta_y = ly / Teuchos::as<real_type>(ny - 1);
        delta_z = lz / Teuchos::as<real_type>(nz - 1);
        Kokkos::parallel_for(
            "coords2d",
            range_type(0, NumMyElements),
            KOKKOS_LAMBDA(const LocalOrdinal i) {
              auto ixy = lclMap.getGlobalElement(i) % (nx * ny);
              auto iz  = (lclMap.getGlobalElement(i) - ixy) / (nx * ny);
              auto ix  = ixy % nx;
              auto iy  = (ixy - ix) / nx;

              lclCoords(i, 0) = delta_x * (real_type)ix;
              lclCoords(i, 1) = delta_y * (real_type)iy;
              lclCoords(i, 2) = delta_z * (real_type)iz;
            });
      }
    }
    return coordinates;

  }  // CreateCartesianCoordinates()

  template <typename GlobalOrdinal>
  static void getSubdomainData(GlobalOrdinal n, GlobalOrdinal m, GlobalOrdinal i, GlobalOrdinal& start, GlobalOrdinal& end) {
    TEUCHOS_TEST_FOR_EXCEPTION(i >= m, std::runtime_error, "Incorrect input parameter to getSubdomainData: m = " + toString(m) + ", i = " + toString(i));

    // If the number of points is not multiple of the number of subdomains, we assign
    // extra points to the first few subdomains
    GlobalOrdinal minWidth = Teuchos::as<GlobalOrdinal>(floor(Teuchos::as<double>(n) / m));
    GlobalOrdinal numWides = n - m * minWidth;
    if (i < numWides) {
      start = i * (minWidth + 1);
      end   = start + minWidth + 1;
    } else {
      start = numWides * (minWidth + 1) + (i - numWides) * minWidth;
      end   = start + minWidth;
    }
  }

};  // class Utils
}  // namespace Xpetra
}  // namespace Galeri

#endif  // ifndef GALERI_XPETRAUTILS_HPP
