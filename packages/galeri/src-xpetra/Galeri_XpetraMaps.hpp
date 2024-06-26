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

  Differences with Galeri1:
   - This function only supports mapType=Cartesian2D and Cartesian3D
   - Parameters that are not set by user but computed inside of this function are saved on the parameter list. This allows users to retrieve these parameters after the creation of the map.
     (In Galeri1, such parameters was set to -1 instead)
*/

//TODO: Check is some GlobalOrdinal are avoidable.

#ifndef GALERI_XPETRAMAPS_HPP
#define GALERI_XPETRAMAPS_HPP

#include <string>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_ConfigDefs.h"

#include "Galeri_XpetraCartesian.hpp"

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

#ifdef HAVE_XPETRA_EPETRA

    template <class LocalOrdinal, class GlobalOrdinal> struct privateCreateMapEpetra {
      static RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
         throw "Galeri::Xpetra::privateCreateMapEpetra: no default implementation";
      }
    };

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
    template <> struct privateCreateMapEpetra<int, int> {
      static RCP< ::Xpetra::Map<int, int, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
        return Galeri::Xpetra::CreateMap<int, int, ::Xpetra::EpetraMapT<int, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> >(mapType, comm, list);
      }
    };
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
    template <> struct privateCreateMapEpetra<int, long long> {
      static RCP< ::Xpetra::Map<int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
        return Galeri::Xpetra::CreateMap<int, long long, ::Xpetra::EpetraMapT<long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> >(mapType, comm, list);
      }
    };
#endif

#endif

#ifdef HAVE_GALERI_XPETRA
    //! Map creation function (for Xpetra::Map with UnderlyingLib parameter)
    template <class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< ::Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
#ifdef HAVE_XPETRA_TPETRA
      if (lib == ::Xpetra::UseTpetra)
        return CreateMap<LocalOrdinal, GlobalOrdinal, ::Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(mapType, comm, list);
#endif
#ifdef HAVE_XPETRA_EPETRA
      if (lib == ::Xpetra::UseEpetra) {
        return CreateMap<LocalOrdinal, GlobalOrdinal, ::Xpetra::EpetraMapT<GlobalOrdinal, Node> >(mapType, comm, list);
      }
#endif
      XPETRA_FACTORY_END;
    }

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
    //! Map creation function (for Xpetra::Map with UnderlyingLib parameter)
    RCP< ::Xpetra::Map<int, int, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > CreateMap(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {

      typedef int LocalOrdinal;
      typedef int GlobalOrdinal;
      typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == ::Xpetra::UseTpetra)
        return CreateMap<int, GlobalOrdinal, ::Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(mapType, comm, list);
#endif
#ifdef HAVE_XPETRA_EPETRA
      if (lib == ::Xpetra::UseEpetra)
        return CreateMap<int, GlobalOrdinal, ::Xpetra::EpetraMapT<GlobalOrdinal, Node> >(mapType, comm, list);
#endif
      XPETRA_FACTORY_END;
    }
#endif // XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
    //! Map creation function (for Xpetra::Map with UnderlyingLib parameter)
    RCP< ::Xpetra::Map<int, long long, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType> > CreateMap64(::Xpetra::UnderlyingLib lib, const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {

      typedef int LocalOrdinal;
      typedef long long GlobalOrdinal;
      typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifdef HAVE_XPETRA_TPETRA
      if (lib == ::Xpetra::UseTpetra)
        return CreateMap<int, GlobalOrdinal, ::Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> >(mapType, comm, list);
#endif
#ifdef HAVE_XPETRA_EPETRA
      if (lib == ::Xpetra::UseEpetra)
        return CreateMap<int, GlobalOrdinal, ::Xpetra::EpetraMapT<GlobalOrdinal, Node> >(mapType, comm, list);
#endif
      XPETRA_FACTORY_END;
    }
#endif // XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
#endif // HAVE_GALERI_XPETRA


  template <class LocalOrdinal, class GlobalOrdinal, class Map>
  RCP<Map> CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list) {
    GlobalOrdinal n = -1, nx = -1, ny = -1, nz = -1, mx = -1, my = -1, mz = -1;

    // Get matrix dimensions
    if (list.isParameter("n")) {
      if (list.isType<int>("n"))
        n = Teuchos::as<GlobalOrdinal>(list.get<int>("n"));
      else
        n = list.get<GlobalOrdinal>("n");
    }
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
    if (list.isParameter("mx")) {
      if (list.isType<int>("mx"))
        mx = Teuchos::as<GlobalOrdinal>(list.get<int>("mx"));
      else
        mx = list.get<GlobalOrdinal>("mx");
    }
    if (list.isParameter("my")) {
      if (list.isType<int>("my"))
        my = Teuchos::as<GlobalOrdinal>(list.get<int>("my"));
      else
        my = list.get<GlobalOrdinal>("my");
    }
    if (list.isParameter("mz")) {
      if (list.isType<int>("mz"))
        mz = Teuchos::as<GlobalOrdinal>(list.get<int>("mz"));
      else
        mz = list.get<GlobalOrdinal>("mz");
    }

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

      return Maps::Cartesian1D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, mx, list);

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

      return Maps::Cartesian2D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, ny, mx, my, list);

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

      return Maps::Cartesian3D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, ny, nz, mx, my, mz, list);

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
