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
#ifndef GALERI_XPETRACARTESIAN2D_HPP
#define GALERI_XPETRACARTESIAN2D_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Galeri_Exception.h"
#include "Galeri_MapTraits.hpp"
#include "Galeri_XpetraUtils.hpp"

#ifdef HAVE_GALERI_TPETRA //TODO: this macro is not defined
#include <Tpetra_Map.hpp>
#endif

namespace Galeri {

  namespace Xpetra {

    namespace Maps {

      typedef size_t global_size_t;

      //TODO: avoid using GlobalOrdinal everywhere?

      template <class LocalOrdinal, class GlobalOrdinal, class Map>
      Teuchos::RCP<Map> Cartesian1D(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const GlobalOrdinal nx,
                                    const GlobalOrdinal mx,
                                    Teuchos::ParameterList & list) {
        if (nx <= 0 || mx <= 0 || (mx > nx))
          throw Exception(__FILE__, __LINE__,
                          "Incorrect input parameter to Maps::Cartesian1D()",
                          "nx = " + toString(nx) +
                          ", mx = " + toString(mx));

        typedef GlobalOrdinal GO;
        typedef LocalOrdinal  LO;

        int myPID = comm->getRank();

        GO startx, endx;
        Utils::getSubdomainData<GO>(nx, mx, myPID, startx, endx);

        list.set("lnx", Teuchos::as<LO>(endx - startx));
        list.set("lny", -1);
        list.set("lnz", -1);

        size_t numMyElements = endx - startx;
        std::vector<GO> myGlobalElements(numMyElements);

        size_t count = 0;
        for (GO i = startx; i < endx; i++)
          myGlobalElements[count++] = i;

        const Teuchos::ArrayView<const GO> elementList(myGlobalElements);

        global_size_t numGlobalElements = nx;
        return MapTraits<GO,Map>::Build(numGlobalElements, elementList, 0/*indexBase*/, comm /*TODO:node*/);
      }

      template <class LocalOrdinal, class GlobalOrdinal, class Map>
      Teuchos::RCP<Map> Cartesian2D(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const GlobalOrdinal nx, const GlobalOrdinal ny,
                                    const GlobalOrdinal mx, const GlobalOrdinal my,
                                    Teuchos::ParameterList & list) {
        if (nx <= 0 || ny <= 0 || mx <= 0 || my <= 0 || (mx > nx) || (my > ny))
          throw(Exception(__FILE__, __LINE__,
                          "Incorrect input parameter to Maps::Cartesian2D()",
                          "nx = " + toString(nx) +
                          ", ny = " + toString(ny) +
                          ", mx = " + toString(mx) +
                          ", my = " + toString(my)));

        typedef GlobalOrdinal GO;
        typedef LocalOrdinal  LO;

        int myPID = comm->getRank();

        GO startx, starty, endx, endy;
        Utils::getSubdomainData(nx, mx, myPID % mx, startx, endx);
        Utils::getSubdomainData(ny, my, myPID / mx, starty, endy);

        list.set("lnx", Teuchos::as<LO>(endx - startx));
        list.set("lny", Teuchos::as<LO>(endy - starty));
        list.set("lnz", -1);

        size_t numMyElements = (endx - startx) * (endy - starty);
        std::vector<GO> myGlobalElements(numMyElements);

        size_t count = 0;
        for (GO i = startx; i < endx; i++)
          for (GO j = starty; j < endy; j++)
            myGlobalElements[count++] = j*nx + i;

        const Teuchos::ArrayView<const GO> elementList(myGlobalElements);

        global_size_t numGlobalElements = nx * ny;
        return MapTraits<GO,Map>::Build(numGlobalElements, elementList, 0/*indexBase*/, comm /*TODO:node*/);
      }

      template <class LocalOrdinal, class GlobalOrdinal, class Map>
      Teuchos::RCP<Map> Cartesian3D(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const GlobalOrdinal nx, const GlobalOrdinal ny, const GlobalOrdinal nz,
                                    const GlobalOrdinal mx, const GlobalOrdinal my, const GlobalOrdinal mz,
                                    Teuchos::ParameterList & list) {
        if (nx <= 0 || ny <= 0 || nz <= 0 ||
            mx <= 0 || my <= 0 || mz <= 0 ||
            (mx > nx) || (my > ny) || (mz > nz))
          throw Exception(__FILE__, __LINE__,
                          "Incorrect input parameter to Maps::Cartesian3D()",
                          "nx = " + toString(nx) +
                          ", ny = " + toString(ny) +
                          ", nz = " + toString(nz) +
                          ", mx = " + toString(mx) +
                          ", my = " + toString(my) +
                          ", mz = " + toString(mz));

        typedef GlobalOrdinal GO;
        typedef LocalOrdinal  LO;

        GO mxy = mx * my;

        int myPID = comm->getRank();

        GO startx, starty, startz, endx, endy, endz;
        Utils::getSubdomainData(nx, mx, (myPID % mxy) % mx, startx, endx);
        Utils::getSubdomainData(ny, my, (myPID % mxy) / mx, starty, endy);
        Utils::getSubdomainData(nz, mz,  myPID / mxy      , startz, endz);

        list.set("lnx", Teuchos::as<LO>(endx - startx));
        list.set("lny", Teuchos::as<LO>(endy - starty));
        list.set("lnz", Teuchos::as<LO>(endz - startz));

        size_t numMyElements = (endx - startx) * (endy - starty) * (endz - startz);
        std::vector<GO> myGlobalElements(numMyElements);

        size_t count = 0;
        for (GO i = startx; i < endx; i++)
          for (GO j = starty; j < endy; j++)
            for (GO k = startz; k < endz; k++)
              myGlobalElements[count++] = k*(nx*ny) + j*nx + i;

        const Teuchos::ArrayView<const GO> elementList(myGlobalElements);

        global_size_t numGlobalElements = nx * ny * nz;
        return MapTraits<GO,Map>::Build(numGlobalElements, elementList, 0/*indexBase*/, comm /*TODO:node*/);
      }

    } // namespace Maps
  } // namespace Xpetra
} // namespace Galeri

#endif
