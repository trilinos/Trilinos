/*
  Direct translation of Galeri coordinate generator.
*/
#ifndef MUELU_GALLERYUTILS_HPP
#define MUELU_GALLERYUTILS_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef XPETRA_ENABLED
// needed for the specialized traits:
#include "Xpetra_Map.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#endif

#include <iostream>

namespace MueLu {
  
  class GalleryUtils {

#   include "Xpetra_UseShortNames.hpp"

    public:

    template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Map>
    static Teuchos::RCP<MultiVector>
    CreateCartesianCoordinates(std::string const &coordType, RCP<const Map> const & map, Teuchos::ParameterList& list)
    {
      Teuchos::RCP<MultiVector> coordinates;

      Scalar delta_x, delta_y, delta_z;

      Scalar lx = list.get("lx", 1.0);
      Scalar ly = list.get("ly", 1.0);
      Scalar lz = list.get("lz", 1.0);

      GlobalOrdinal nx = list.get("nx", -1);
      GlobalOrdinal ny = list.get("ny", -1);
      GlobalOrdinal nz = list.get("nz", -1);

      GlobalOrdinal ix, iy, iz;

      LocalOrdinal NumMyElements = map->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

      if (coordType == "1D") {
        coordinates = MultiVectorFactory::Build(map,1);
        Teuchos::ArrayRCP<ArrayRCP<Scalar> > Coord(1);
        Coord[0] = coordinates->getDataNonConst(0);

        delta_x = lx / (nx - 1);

        for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i) {
          ix = MyGlobalElements[i];
          Coord[0][i] = delta_x * ix;
        }

      } else if (coordType == "2D") {

        coordinates = MultiVectorFactory::Build(map,2);
        Teuchos::ArrayRCP<ArrayRCP<Scalar> > Coord(2);
        Coord[0] = coordinates->getDataNonConst(0);
        Coord[1] = coordinates->getDataNonConst(1);

        delta_x = lx / (nx - 1);
        delta_y = ly / (ny - 1);

        for (LocalOrdinal i = 0 ; i < NumMyElements ; ++i)
        {
          ix = MyGlobalElements[i] % nx;
          iy = (MyGlobalElements[i] - ix) / nx;

          Coord[0][i] = delta_x * ix;
          Coord[1][i] = delta_y * iy;
        }

      } else if (coordType == "3D") {

        coordinates = MultiVectorFactory::Build(map,3);
        Teuchos::ArrayRCP<ArrayRCP<Scalar> > Coord(3);
        Coord[0] = coordinates->getDataNonConst(0);
        Coord[1] = coordinates->getDataNonConst(1);
        Coord[2] = coordinates->getDataNonConst(2);

        delta_x = lx / (nx - 1);
        delta_y = ly / (ny - 1);
        delta_z = lz / (nz - 1);

        for (LocalOrdinal i = 0 ; i < NumMyElements ; i++)
        {
          GlobalOrdinal ixy = MyGlobalElements[i] % (nx * ny);
          iz = (MyGlobalElements[i] - ixy) / (nx * ny);

          ix = ixy % nx;
          iy = (ixy - ix) / ny;

          Coord[0][i] = delta_x * ix;
          Coord[1][i] = delta_y * iy;
          Coord[2][i] = delta_z * iz;
        }

      } else {

        throw(std::runtime_error("in MueLu::GalleryUtils : `coordType' has incorrect value (" + coordType + ")"));

      } //if (coordType == ...

      return coordinates;

    } // CreateCartesianCoordinates()

  }; // class GalleryUtils
} // namespace MueLu

#define MUELU_GALLERYUTILS_SHORT

#endif //ifndef MUELU_GALLERYUTILS_HPP
