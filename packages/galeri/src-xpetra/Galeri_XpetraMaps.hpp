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

  Differences with Galeri1:
   - This function only supports mapType=Cartesian2D and Cartesian3D
   - Parameters that are not set by user but computed inside of this function are saved on the parameter list. This allows users to retrieve these parameters after the creation of the map.
     (In Galeri1, such parameters was set to -1 instead)
*/

//TODO: Check is some GlobalOrdinal are avoidable.

#ifndef GALERI_XPETRAMAPS_HPP
#define GALERI_XPETRAMAPS_HPP

#include <string.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_ConfigDefs.h"

#include "Galeri_XpetraCartesian1D.hpp"
#include "Galeri_XpetraCartesian2D.hpp"
#include "Galeri_XpetraCartesian3D.hpp"

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#include <Xpetra_Map.hpp> // for enum UnderlyingLib
#ifdef HAVE_XPETRA_TPETRA
#  include <Xpetra_TpetraMap.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#  include <Xpetra_EpetraMap.hpp>
#endif
#endif // HAVE_GALERI_XPETRA

//
// DECL
//

namespace Galeri {
  namespace Xpetra {

    using Teuchos::RCP;

    //! Map creation function (for Tpetra, Epetra, Xpetra::TpetraMap and Xpetra::EpetraMap)
    template <class LocalOrdinal, class GlobalOrdinal, class Map>
    RCP<Map> CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList & list);

#ifdef HAVE_GALERI_XPETRA
    //! Map creation function (for Xpetra::Map with an UnderlyingLib parameter)
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList & list);
#endif

  } // namespace Xpetra
} // namespace Galeri


//
// DEF
//

namespace Galeri {
  namespace Xpetra {

    using Teuchos::RCP;

#ifdef HAVE_GALERI_XPETRA
    //! Map creation function (for Xpetra::Map with UnderlyingLib parameter)
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
#ifdef HAVE_XPETRA_TPETRA
      if (lib == ::Xpetra::UseTpetra)
        return CreateMap< ::Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(mapType, comm, list);
#endif
      XPETRA_FACTORY_ERROR_IF_EPETRA(lib);
      XPETRA_FACTORY_END;
    }

    //! Map creation function (for Xpetra::Map with UnderlyingLib parameter)
    template <>
    RCP< ::Xpetra::Map<int, int, Kokkos::DefaultNode::DefaultNodeType> > CreateMap<int, int, Kokkos::DefaultNode::DefaultNodeType>(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {

      typedef int LocalOrdinal;
      typedef int GlobalOrdinal;
      typedef Kokkos::DefaultNode::DefaultNodeType Node;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == ::Xpetra::UseTpetra)
        return CreateMap<int, int, ::Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(mapType, comm, list);
#endif
#ifdef HAVE_XPETRA_EPETRA
      if (lib == ::Xpetra::UseEpetra)
        return CreateMap<int, int, ::Xpetra::EpetraMap>(mapType, comm, list);
#endif
      XPETRA_FACTORY_END;
    }
#endif // HAVE_GALERI_XPETRA


  template <class LocalOrdinal, class GlobalOrdinal, class Map>
  RCP<Map> CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
    GlobalOrdinal n = -1, nx = -1, ny = -1, nz = -1, mx = -1, my = -1, mz = -1;
    // Get matrix dimensions
    if (list.isParameter("n"))  n  = list.get<GlobalOrdinal>("n");
    if (list.isParameter("nx")) nx = list.get<GlobalOrdinal>("nx");
    if (list.isParameter("ny")) ny = list.get<GlobalOrdinal>("ny");
    if (list.isParameter("nz")) nz = list.get<GlobalOrdinal>("nz");
    if (list.isParameter("mx")) nz = list.get<GlobalOrdinal>("mx");
    if (list.isParameter("my")) nz = list.get<GlobalOrdinal>("my");
    if (list.isParameter("mz")) nz = list.get<GlobalOrdinal>("mz");

    if (mapType == "Cartesian1D") {
      if (nx == -1) {
        if (n <= 0)
          throw Exception(__FILE__, __LINE__, "If nx is not set, then n must be set");

        nx = n;
        list.set("nx", nx);
      }
      list.set("n", nx);

      if (mx == -1) {
        mx = Teuchos::as<GlobalOrdinal>(comm->getSize());
        list.set("mx", mx);
      }

      return Maps::Cartesian1D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, mx);

    } else if (mapType == "Cartesian2D") {
      if (nx == -1 || ny == -1) {
        if (n <= 0)
          throw Exception(__FILE__, __LINE__, "If nx or ny are not set, then n must be set");

        nx = ny = Teuchos::as<GlobalOrdinal>(sqrt(Teuchos::as<double>(n)));
        if (nx * ny != n)
          throw Exception(__FILE__, __LINE__,
                          "The number of global elements (n) must be",
                          "a perfect square, otherwise set nx and ny");

        list.set("nx", nx);
        list.set("ny", ny);
      }
      list.set("n", nx * ny);

      if (mx == -1 || my == -1) {
        // Simple method to try to find mx and my such that mx*my = NumProc
        mx = Teuchos::as<GlobalOrdinal>(sqrt(Teuchos::as<double>(comm->getSize())+.001));
        my = comm->getSize()/mx;

        while ((mx*my) != comm->getSize())
          my = comm->getSize()/(--mx);

        list.set("mx", mx);
        list.set("my", my);
      }

      return Maps::Cartesian2D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, ny, mx, my);

    } else if (mapType == "Cartesian3D") {
      if (nx == -1 || ny == -1 || nz == -1) {
        if (n <= 0)
          throw Exception(__FILE__, __LINE__, "If nx or ny or nz are not set, then n must be set");

        nx = ny = nz = Teuchos::as<GlobalOrdinal>(pow(Teuchos::as<double>(n), 0.333334));
        if (nx * ny * nz != n)
          throw Exception(__FILE__, __LINE__,
                          "The number of global elements n (" +
                          toString(n) + ") must be",
                          "a perfect cube, otherwise set nx, ny and nz") ;

        list.set("nx", nx);
        list.set("ny", ny);
        list.set("nz", nz);
      }
      list.set("n", nx * ny * nz);

      if (mx == -1 || my == -1 || mz == -1) {
        mx = my = mz = Teuchos::as<GlobalOrdinal>(pow(Teuchos::as<double>(comm->getSize()), 0.333334));

        if (mx * my * mz != comm->getSize())  {
          // Simple method to find a set of processor assignments
          mx = my = mz = 1;

          int ProcTemp = comm->getSize();
          GlobalOrdinal factors[50];
          for (GlobalOrdinal jj = 0; jj < 50; jj++) factors[jj] = 0;
          for (GlobalOrdinal jj = 2; jj < 50; jj++) {
            bool flag = true;
            while (flag) {
              GlobalOrdinal temp = ProcTemp/jj;
              if (temp*jj == ProcTemp) {
                factors[jj]++;
                ProcTemp = temp;

              } else {
                flag = false;
              }
            }
          }
          mx = ProcTemp;
          for (GlobalOrdinal jj = 50-1; jj > 0; jj--) {
            while (factors[jj] != 0) {
              if      ((mx <= my) && (mx <= mz)) mx = mx*jj;
              else if ((my <= mx) && (my <= mz)) my = my*jj;
              else                               mz = mz*jj;
              factors[jj]--;
            }
          }
        }
        list.set("mx", mx);
        list.set("my", my);
        list.set("mz", mz);

      } else {
        if (mx * my * mz != comm->getSize())
          throw Exception(__FILE__, __LINE__,
                          "mx * my * mz != number of processes!",
                          "mx = " + toString(mx) + ", my = " + toString(my)
                          + ", mz = " + toString(mz));
      }

      return Maps::Cartesian3D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, ny, nz, mx, my, mz);

    } else {
      throw Exception(__FILE__, __LINE__,
                      "`mapType' has incorrect value (" + mapType + ")",
                      "in input to function CreateMap()",
                      "Check the documentation for a list of valid choices");
    }
  } // CreateMap()

  } // namespace Xpetra
} // namespace Galeri

#endif
