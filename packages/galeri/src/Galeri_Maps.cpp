#include "Galeri_Maps.h"
#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Teuchos_ParameterList.hpp"
#include "Maps/Galeri_Linear.h"
#include "Maps/Galeri_Cartesian2D.h"
#include "Maps/Galeri_Cartesian3D.h"
#include "Maps/Galeri_Random.h"
#include "Maps/Galeri_Interlaced.h"

namespace Galeri {

Epetra_Map* 
CreateMap(string MapType, Epetra_Comm& Comm, Teuchos::ParameterList& List)
{
  // global parameters
  int n = List.get("n", -1);

  // Cycle on map type //
  if (MapType == "Linear")
  {
    // get the parameters
    return(Maps::Linear(Comm, n));
  }
  else if (MapType == "Interlaced")
  {
    return(Maps::Interlaced(Comm, n));
  }
  else if (MapType == "Random")
  {
    return(Maps::Random(Comm, n));
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
      if (nx * ny != n) 
        throw(Exception(__FILE__, __LINE__,
                        "The number of global elements (n) must be",
                        "a perfect square, otherwise set nx and ny"));
    }

    // Get the number of domains
    int mx = List.get("mx", -1);
    int my = List.get("my", -1);

    if (mx == -1 || my == -1) 
    {
      mx = (int)sqrt((double)(Comm.NumProc()));
      my = mx;

      if (mx * my != Comm.NumProc()) 
        throw(Exception(__FILE__, __LINE__,
                        "Number of processes must be a perfect square", 
                        "otherwise set mx and my"));
    } 

    return(Maps::Cartesian2D(Comm, nx, ny, mx, my));
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
      if (nx * ny * nz != n) 
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

      if (mx * my * mz != Comm.NumProc()) 
        throw(Exception(__FILE__, __LINE__,
                        "Number of processes must be a perfect cube", 
                        "otherwise set mx, my and mz"));
    } 
    else 
    {
      if (mx * my * mz != Comm.NumProc()) 
        throw(Exception(__FILE__, __LINE__, 
                        "mx * my * mz != number of processes!",
                        "mx = " + toString(mx) + ", my = " + toString(my)
                        + ", mz = " + toString(mz)));
    }

    return(Maps::Cartesian3D(Comm, nx, ny, nz, mx, my, mz));
  }
  else 
  {
    throw(Exception(__FILE__, __LINE__,
                    "`MapType' has incorrect value (" + MapType + ")",
                    "in input to function CreateMap()",
                    "Check the documentation for a list of valid choices"));
  }
} // CreateMap()

} // namespace Galeri
