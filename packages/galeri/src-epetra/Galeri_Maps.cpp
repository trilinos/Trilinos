// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_Maps.h"
#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_ConfigDefs.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Linear.h"
#include "Galeri_NodeCartesian2D.h"
#include "Galeri_Cartesian2D.h"
#include "Galeri_Cartesian3D.h"
#include "Galeri_Random.h"
#include "Galeri_Interlaced.h"

namespace Galeri {

template<typename int_type>
Epetra_Map* 
TCreateMap(std::string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List)
{
  // global parameters
  int_type n = List.get("n", -1);

  // Cycle on map type //
  if (MapType == "Linear")
  {
    // get the parameters
    return(Maps::TLinear<int_type>(Comm, n));
  }
  else if (MapType == "Interlaced")
  {
    return(Maps::TInterlaced<int_type>(Comm, n));
  }
  else if (MapType == "Random")
  {
    return(Maps::TRandom<int_type>(Comm, n));
  }
  else if (MapType == "Cartesian2D")
  {
    // Get matrix dimension
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    if (nx == -1 || ny == -1) 
    {
      if (n <= 0)
          throw(Exception(__FILE__, __LINE__,
                          "If nx or ny are not set, then n must be set"));

      nx = (int)sqrt((double)n);
      ny = nx;
      if ((int_type)nx * (int_type)ny != n) 
        throw(Exception(__FILE__, __LINE__,
                        "The number of global elements (n) must be",
                        "a perfect square, otherwise set nx and ny"));
    }

    // Get the number of domains
    int mx = List.get("mx", -1);
    int my = List.get("my", -1);

    if (mx == -1 || my == -1) 
    {
      mx = (int)(sqrt((double)(Comm.NumProc()+.001)));
      my =  Comm.NumProc()/mx;

      // simple attempt at trying to find an mx and my such 
      // mx*my = NumProc

      while ( (mx*my) != Comm.NumProc() ) {
        mx--;
        my = Comm.NumProc()/mx;
      }
      List.set("mx", mx);
      List.set("my", my);
    } 


    return(Maps::TCartesian2D<int_type>(Comm, nx, ny, mx, my));
  }
  else if (MapType == "NodeCartesian2D")
  {
    // Get matrix dimension
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    if (nx == -1 || ny == -1) 
    {
      if (n <= 0)
          throw(Exception(__FILE__, __LINE__,
                          "If nx or ny are not set, then n must be set"));

      nx = (int)sqrt((double)n);
      ny = nx;
      if ((int_type)nx * (int_type)ny != n) 
        throw(Exception(__FILE__, __LINE__,
                        "The number of global elements (n) must be",
                        "a perfect square, otherwise set nx and ny"));
    }

    // Get the number of nodes
    int NumNodes = List.get("number of nodes",-1);
    int MyNodeID = List.get("node id",-1);
    int ndx = List.get("ndx", -1);
    int ndy = List.get("ndy", -1);
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
    Epetra_Comm * NodeComm= List.get("node communicator",(Epetra_Comm*)0);
    int px = List.get("px", -1);
    int py = List.get("py", -1);
    if (px == -1 || py == -1) 
    {
      px = (int)(sqrt((double)(NodeComm->NumProc()+.001)));
      py = NodeComm->NumProc()/px;

      // simple attempt at trying to find an ndx and ndy such 
      // ndx*ndy = NumNodes

      while ( (px*py) != NodeComm->NumProc() ) {
        px--;
        py = NodeComm->NumProc()/px;
      }
    } 
  
    return(Maps::TNodeCartesian2D<int_type>(Comm, *NodeComm, MyNodeID, nx, ny, ndx, ndy, px,py));
  }  
  else if (MapType == "Cartesian3D")
  {
    // Get matrix dimension
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    int nz = List.get("nz", -1);

    if (nx == -1 || ny == -1 || nz == -1) 
    {
      if (n <= 0)
          throw(Exception(__FILE__, __LINE__,
                          "If nx or ny or nz are not set, then n must be set"));

      nx = (int)pow((double)n, 0.333334);
      ny = nx;
      nz = nx;
      if ((int_type)nx * (int_type)ny * (int_type)nz != n) 
        throw(Exception(__FILE__, __LINE__,
                        "The number of global elements n (" +
                        toString(n) + ") must be",
                        "a perfect cube, otherwise set nx, ny and nz"));
    }

    // Get the number of domains
    int mx = List.get("mx", -1);
    int my = List.get("my", -1);
    int mz = List.get("mz", -1);

    if (mx == -1 || my == -1 || mz == -1) 
    {
      mx = (int)pow((double)Comm.NumProc(), 0.333334);
      my = mx;
      mz = mx;

      if (mx * my * mz != Comm.NumProc())  {
	  // simple attempt to find a set of processor assignments
         mx = 1; my = 1; mz = 1;
         int ProcTemp = Comm.NumProc();
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
      List.set("mx", mx);
      List.set("my", my);
      List.set("mz", mz);
    } 
    else 
    {
      if (mx * my * mz != Comm.NumProc()) 
        throw(Exception(__FILE__, __LINE__, 
                        "mx * my * mz != number of processes!",
                        "mx = " + toString(mx) + ", my = " + toString(my)
                        + ", mz = " + toString(mz)));
    }

    return(Maps::TCartesian3D<int_type>(Comm, nx, ny, nz, mx, my, mz));
  }
  else 
  {
    throw(Exception(__FILE__, __LINE__,
                    "`MapType' has incorrect value (" + MapType + ")",
                    "in input to function CreateMap()",
                    "Check the documentation for a list of valid choices"));
  }
} // TCreateMap()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
CreateMap(std::string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List) {
  return TCreateMap<int>(MapType, Comm, List);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
CreateMap64(std::string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List) {
  return TCreateMap<long long>(MapType, Comm, List);
}
#endif


} // namespace Galeri
