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
   - This function only supports mapType=Cartesian2D
   - Parameters that are not set by user but computed inside of this function are saved on the parameter list. This allows users to retrieve these parameters after the creation of the map.
     (In Galeri1, such parameters was set to -1 instead)
*/
#ifndef GALERI_XPETRAMAPS_HPP
#define GALERI_XPETRAMAPS_HPP


//
// DECL
//

#include <string.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_ConfigDefs.h"

#ifdef HAVE_GALERI_XPETRA
#include "Xpetra_Map.hpp" // for enum UnderlyingLib
#endif

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

#include <string.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Galeri_ConfigDefs.h"
#include "Galeri_XpetraMaps.hpp"

//#include "Galeri_Utils.h"

//#include "Maps/Galeri_XpetraLinear.hpp"
//#include "Maps/Galeri_XpetraNodeCartesian2D.hpp"
#include "Maps/Galeri_XpetraCartesian2D.hpp"
//#include "Maps/Galeri_XpetraCartesian3D.hpp"
//#include "Maps/Galeri_XpetraRandom.hpp"
//#include "Maps/Galeri_XpetraInterlaced.hpp"

#ifdef HAVE_GALERI_XPETRA
#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Exceptions.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMap.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif
#endif // HAVE_GALERI_XPETRA

namespace Galeri {

  namespace Xpetra {

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
  RCP<Map> CreateMap(const std::string & mapType, const Teuchos::RCP<const Teuchos::Comm<int> > & comm, Teuchos::ParameterList & list)
  {
    // global parameters
    // GlobalOrdinal n = list.get("n", (GlobalOrdinal) -1);

    // Cycle on map type //
    /* TODO
    if (mapType == "Linear")
      {
        // get the parameters
        return(Maps::Linear(comm, n));
      }
    else if (mapType == "Interlaced")
      {
        return(Maps::Interlaced(comm, n));
      }
    else if (mapType == "Random")
      {
        return(Maps::Random(comm, n));
      }
    else 
    */ 
    if (mapType == "Cartesian2D")
      {

        // Get matrix dimension
        GlobalOrdinal nx = -1; if (list.isParameter("nx")) nx = list.get<GlobalOrdinal>("nx");
        GlobalOrdinal ny = -1; if (list.isParameter("ny")) ny = list.get<GlobalOrdinal>("ny");
        
        if (nx == -1 || ny == -1) 
          {
            // global parameters
            GlobalOrdinal n = -1; if (list.isParameter("n")) list.get<GlobalOrdinal>("n");

            if (n <= 0)
              throw(Exception(__FILE__, __LINE__,
                              "If nx or ny are not set, then n must be set"));

            nx = (GlobalOrdinal)sqrt((double)n);
            ny = nx;
            if (nx * ny != n) 
              throw(Exception(__FILE__, __LINE__,
                              "The number of global elements (n) must be",
                              "a perfect square, otherwise set nx and ny"));
          
            // Add computed values of 'nx' and 'ny' to the list. Users can retrieve the values after the creation of the map.
            list.get("nx", nx);
            list.get("ny", ny);
          }

        list.get("n", nx * ny); // Add computed values of 'n' to the list if not already there

        // Get the number of domains
        GlobalOrdinal mx = -1; if (list.isParameter("mx")) mx = list.get<GlobalOrdinal>("mx");  //TODO: GlobalOrdinal or LocalOrdinal?
        GlobalOrdinal my = -1; if (list.isParameter("my")) my = list.get<GlobalOrdinal>("my");

        if (mx == -1 || my == -1) 
          {
            mx = (GlobalOrdinal)(sqrt((double)(comm->getSize()+.001)));
            my =  comm->getSize()/mx;

            // simple attempt at trying to find an mx and my such 
            // mx*my = NumProc

            while ( (mx*my) != comm->getSize() ) {
              mx--;
              my = comm->getSize()/mx;
            }

            // Add computed values of 'mx' and 'my' to the list. Users can retrieve the values after the creation of the map.
            list.get("mx", mx);
            list.get("my", my);
          } 

        return(Maps::Cartesian2D<LocalOrdinal, GlobalOrdinal, Map>(comm, nx, ny, mx, my));
      }
    /* TODO
       else if (mapType == "NodeCartesian2D")
      {
        // Get matrix dimension
        int nx = list.get("nx", -1);
        int ny = list.get("ny", -1);

        if (nx == -1 || ny == -1) 
          {
            if (n <= 0)
              throw(Exception(__FILE__, __LINE__,
                              "If nx or ny are not set, then n must be set"));

            nx = (int)sqrt((double)n);
            ny = nx;
            if (nx * ny != n) 
              throw(Exception(__FILE__, __LINE__,
                              "The number of global elements (n) must be",
                              "a perfect square, otherwise set nx and ny"));
          }

        // Get the number of nodes
        int NumNodes = list.get("number of nodes",-1);
        int MyNodeID = list.get("node id",-1);
        int ndx = list.get("ndx", -1);
        int ndy = list.get("ndy", -1);
        if (ndx == -1 || ndy == -1) 
          {
            ndx = (int)(sqrt((double)(NumNodes+.001)));
            ndy = NumNodes/ndx;
            while ( (ndx*ndy) != NumNodes ) {
              ndx--;
              ndy = NumNodes/ndx;
            }
          } 

        // Get the number of processors per node
        Epetra_comm * Nodecomm= list.get("node communicator",(Epetra_comm*)0);
        int px = list.get("px", -1);
        int py = list.get("py", -1);
        if (px == -1 || py == -1) 
          {
            px = (int)(sqrt((double)(Nodecomm->NumProc()+.001)));
            py = Nodecomm->NumProc()/px;

            // simple attempt at trying to find an ndx and ndy such 
            // ndx*ndy = NumNodes

            while ( (px*py) != Nodecomm->NumProc() ) {
              px--;
              py = Nodecomm->NumProc()/px;
            }
          } 
  
        return(Maps::NodeCartesian2D(comm, *Nodecomm, MyNodeID, nx, ny, ndx, ndy, px,py));
      }  
    else if (mapType == "Cartesian3D")
      {
        // Get matrix dimension
        int nx = list.get("nx", -1);
        int ny = list.get("ny", -1);
        int nz = list.get("nz", -1);

        if (nx == -1 || ny == -1 || nz == -1) 
          {
            if (n <= 0)
              throw(Exception(__FILE__, __LINE__,
                              "If nx or ny or nz are not set, then n must be set"));

            nx = (int)pow((double)n, 0.333334);
            ny = nx;
            nz = nx;
            if (nx * ny * nz != n) 
              throw(Exception(__FILE__, __LINE__,
                              "The number of global elements n (" +
                              toString(n) + ") must be",
                              "a perfect cube, otherwise set nx, ny and nz"));
          }

        // Get the number of domains
        int mx = list.get("mx", -1);
        int my = list.get("my", -1);
        int mz = list.get("mz", -1);

        if (mx == -1 || my == -1 || mz == -1) 
          {
            mx = (int)pow((double)comm->getSize(), 0.333334);
            my = mx;
            mz = mx;

            if (mx * my * mz != comm->getSize())  {
              // simple attempt to find a set of processor assignments
              mx = 1; my = 1; mz = 1;
              int ProcTemp = comm->getSize();
              int factors[50];
              for (int jj = 0; jj < 50; jj++) factors[jj] = 0;
              for (int jj = 2; jj < 50; jj++) {
                int flag = 1;
                while (flag == 1) {
                  int temp = ProcTemp/jj;
                  if (temp*jj == ProcTemp) {
                    factors[jj]++; ProcTemp = temp;
                  }
                  else flag = 0;
                }
              }
              mx = ProcTemp;
              for (int jj = 50-1; jj > 0; jj--) {
                while (factors[jj] != 0) {
                  if (  (mx <= my) && (mx <= mz) ) mx = mx*jj;
                  else if (  (my <= mx) && (my <= mz) ) my = my*jj;
                  else mz = mz*jj;
                  factors[jj]--;
                }
              }
        
            }
          } 
        else 
          {
            if (mx * my * mz != comm->getSize()) 
              throw(Exception(__FILE__, __LINE__, 
                              "mx * my * mz != number of processes!",
                              "mx = " + toString(mx) + ", my = " + toString(my)
                              + ", mz = " + toString(mz)));
          }

        return(Maps::Cartesian3D(comm, nx, ny, nz, mx, my, mz));
      }
    */
    else 
      {
        throw(Exception(__FILE__, __LINE__,
                        "`mapType' has incorrect value (" + mapType + ")",
                        "in input to function CreateMap()",
                        "Check the documentation for a list of valid choices"));
      }
  } // CreateMap()

  } // namespace Xpetra
} // namespace Galeri

#endif
