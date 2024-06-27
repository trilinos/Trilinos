// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_CrsMatrices.h"
#include "Galeri_Exception.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "CrsMatrices/Galeri_Diag.h"
#include "CrsMatrices/Galeri_Tridiag.h"
#include "CrsMatrices/Galeri_Lehmer.h"
#include "CrsMatrices/Galeri_Minij.h"
#include "CrsMatrices/Galeri_Ris.h"
#include "CrsMatrices/Galeri_Hilbert.h"
#include "CrsMatrices/Galeri_JordanBlock.h"
#include "CrsMatrices/Galeri_Cauchy.h"
#include "CrsMatrices/Galeri_Fielder.h"
#include "CrsMatrices/Galeri_Hanowa.h"
#include "CrsMatrices/Galeri_KMS.h"
#include "CrsMatrices/Galeri_Parter.h"
#include "CrsMatrices/Galeri_Pei.h"
#include "CrsMatrices/Galeri_Ones.h"
#include "CrsMatrices/Galeri_Vander.h"
#include "CrsMatrices/Galeri_Cross2D.h"
#include "CrsMatrices/Galeri_Cross3D.h"
#include "CrsMatrices/Galeri_Star2D.h"
#include "CrsMatrices/Galeri_UniFlow2D.h"
#include "CrsMatrices/Galeri_Recirc2D.h"
#include "CrsMatrices/Galeri_BentPipe2D.h"
#include "CrsMatrices/Galeri_Stretched2D.h"
#include "CrsMatrices/Galeri_Laplace1DNeumann.h"
#include "CrsMatrices/Galeri_BigCross2D.h"
#include "CrsMatrices/Galeri_BigStar2D.h"

namespace Galeri {

Epetra_CrsMatrix*
CreateCrsMatrix(const std::string MatrixType, const Epetra_Map* Map,
                Teuchos::ParameterList& List)
{
  // =============== //
  // MATLAB MATRICES //
  // =============== //
  //
  if (MatrixType == "Diag")
  {
    double a = List.get("a", 1.0);
    return(Matrices::Diag(Map, a));
  }
  else if (MatrixType == "Tridiag")
  {
    double a = List.get("a", 2.0);
    double b = List.get("b", -1.0);
    double c = List.get("c", -1.0);
    return(Matrices::Tridiag(Map, a, b, c));
  }
  else if (MatrixType == "Lehmer")
  {
    return(Matrices::Lehmer(Map));
  }
  else if (MatrixType == "Minij")
  {
    return(Matrices::Minij(Map));
  }
  else if (MatrixType == "Ris")
  {
    return(Matrices::Ris(Map));
  }
  else if (MatrixType == "Hilbert")
  {
    return(Matrices::Hilbert(Map));
  }
  else if (MatrixType == "JordanBlock")
  {
    double lambda = List.get("lambda", 0.1);
    return(Matrices::JordanBlock(Map, lambda));
  }
  else if (MatrixType == "Cauchy")
  {
    return(Matrices::Cauchy(Map));
  }
  else if (MatrixType == "Fielder")
  {
    return(Matrices::Fielder(Map));
  }
  else if (MatrixType == "Hanowa")
  {
    double lambda = List.get("a", -1.0);
    return(Matrices::Hanowa(Map, lambda));
  }
  else if (MatrixType == "KMS")
  {
    double lambda = List.get("rho", -0.5);
    return(Matrices::KMS(Map, lambda));
  }
  else if (MatrixType == "Parter")
  {
    return(Matrices::Parter(Map));
  }
  else if (MatrixType == "Pei")
  {
    double lambda = List.get("alpha", 1.0);
    return(Matrices::Pei(Map, lambda));
  }
  else if (MatrixType == "Ones")
  {
    double lambda = List.get("a", 1.0);
    return(Matrices::Ones(Map, lambda));
  }
  else if (MatrixType == "Vander")
  {
    double lambda = List.get("lambda", 1.0);
    return(Matrices::Vander(Map, lambda));
  }

  // ========================== //
  // FINITE DIFFERENCE MATRICES //
  // ========================== //
  else if (MatrixType == "Cross2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    double a = List.get("a", 4.0);
    double b = List.get("b", -1.0);
    double c = List.get("c", -1.0);
    double d = List.get("d", -1.0);
    double e = List.get("e", -1.0);

    return(Matrices::Cross2D(Map, nx, ny, a, b, c, d, e));
  }
  else if (MatrixType == "Star2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    double a = List.get("a", 8.0);
    double b = List.get("b", -1.0);
    double c = List.get("c", -1.0);
    double d = List.get("d", -1.0);
    double e = List.get("e", -1.0);
    double z1 = List.get("z1", -1.0);
    double z2 = List.get("z2", -1.0);
    double z3 = List.get("z3", -1.0);
    double z4 = List.get("z4", -1.0);

    return(Matrices::Star2D(Map, nx, ny, a, b, c, d, e, z1, z2, z3, z4));
  }
  else if (MatrixType == "BigStar2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    double a = List.get("a", 20.0);
    double b = List.get("b", -8.0);
    double c = List.get("c", -8.0);
    double d = List.get("d", -8.0);
    double e = List.get("e", -8.0);
    double z1 = List.get("z1", 2.0);
    double z2 = List.get("z2", 2.0);
    double z3 = List.get("z3", 2.0);
    double z4 = List.get("z4", 2.0);
    double bb = List.get("bb", 1.0);
    double cc = List.get("cc", 1.0);
    double dd = List.get("dd", 1.0);
    double ee = List.get("ee", 1.0);

    return(Matrices::BigStar2D(Map, nx, ny, a, b, c, d, e,
                               z1, z2, z3, z4, bb, cc, dd, ee));
  }
  else if (MatrixType == "Biharmonic2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    // double a = List.get("a", 20.0); // unused
    // double b = List.get("b", -8.0);
    // double c = List.get("c", -8.0);
    // double d = List.get("d", -8.0);
    // double e = List.get("e", -8.0);
    // double z1 = List.get("z1", 2.0);
    // double z2 = List.get("z2", 2.0);
    // double z3 = List.get("z3", 2.0);
    // double z4 = List.get("z4", 2.0);
    // double bb = List.get("bb", 1.0);
    // double cc = List.get("cc", 1.0);
    // double dd = List.get("dd", 1.0);
    // double ee = List.get("ee", 1.0);

    return(Matrices::BigStar2D(Map, nx, ny, 20, -8, -8, -8, -8,
                               2, 2, 2, 2, 1, 1, 1, 1));
  }
  else if (MatrixType == "BigCross2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);

    double a = List.get("a", 60.0);
    double b = List.get("b", -16.0);
    double c = List.get("c", -16.0);
    double d = List.get("d", -16.0);
    double e = List.get("e", -16.0);
    double bb = List.get("bb", 1.0);
    double cc = List.get("cc", 1.0);
    double dd = List.get("dd", 1.0);
    double ee = List.get("ee", 1.0);

    return(Matrices::BigCross2D(Map, nx, ny, a, b, c, d, e,
                                bb, cc, dd, ee));
  }
  else if (MatrixType == "Cross3D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    int nz = List.get("nz", -1);

    double a = List.get("a", 6.0);
    double b = List.get("b", -1.0);
    double c = List.get("c", -1.0);
    double d = List.get("d", -1.0);
    double e = List.get("e", -1.0);
    double f = List.get("f", -1.0);
    double g = List.get("g", -1.0);

    return(Matrices::Cross3D(Map, nx, ny, nz, a, b, c, d, e, f, g));
  }
  else if (MatrixType == "Laplace1D")
  {
    return(Matrices::Tridiag(Map, 2.0, -1.0, -1.0));
  }
  else if (MatrixType == "Laplace1DNeumann")
  {
    return(Matrices::Laplace1DNeumann(Map));
  }
  else if (MatrixType == "Laplace2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    return(Matrices::Cross2D(Map, nx, ny, 4.0, -1.0, -1.0, -1.0, -1.0));
  }
  else if (MatrixType == "Laplace2DFourthOrder")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    return(Matrices::BigCross2D(Map, nx, ny, 60.0, -16.0, -16.0, -16.0, -16.0,
                                1.0, 1.0, 1.0, 1.0));
  }
  else if (MatrixType == "Stretched2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    double epsilon = List.get("epsilon", 0.1);
    // double lx = List.get("lx", 1.0); // unused
    // double ly = List.get("ly", 1.0);

    return(Matrices::Stretched2D(Map, nx, ny, epsilon));
  }
  else if (MatrixType == "UniFlow2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    double conv = List.get("conv", 1.0);
    double diff = List.get("diff", 1.0e-5);
    double alpha = List.get("alpha", 0.0);
    double lx = List.get("lx", 1.0);
    double ly = List.get("ly", 1.0);

    return(Matrices::UniFlow2D(Map, nx, ny, lx, ly, conv, diff, alpha));
  }
  else if (MatrixType == "Recirc2D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    double conv = List.get("conv", 1.0);
    double diff = List.get("diff", 1.0e-5);
    double lx = List.get("lx", 1.0);
    double ly = List.get("ly", 1.0);

    return(Matrices::Recirc2D(Map, nx, ny, lx, ly, conv, diff));
  }
  else if (MatrixType == "BentPipe2D")
  {
    /*
     * to visualize the convective field in MATLAB:
     >> [x,y] = meshgrid(0:0.02:1, 0:0.02:1);
     >> bx = 4 * x .* (x - 1) .* (1 - 2 * y);
     >> by = -4 * y .* (y - 1) .* (1 - 2 * x);
     >> % select the starting point of the streamlines
     >> streamline(stream2(x,y,bx,by,[0.5 0.5 0.5], [0.2 0.4 0.6]))
    */
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    if (nx == -1 || ny == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)sqrt((double)n);
      ny = nx;
      if (((long long) nx) * ny != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }

    double conv = List.get("conv", 1.0);
    double diff = List.get("diff", 1.0e-5);
    double lx = List.get("lx", 1.0);
    double ly = List.get("ly", 1.0);

    return(Matrices::BentPipe2D(Map, nx, ny, lx, ly, conv, diff));
  }
  else if (MatrixType == "Laplace3D")
  {
    int nx = List.get("nx", -1);
    int ny = List.get("ny", -1);
    int nz = List.get("nz", -1);
    if (nx == -1 || ny == -1 || nz == -1)
    {
      long long n = Map->NumGlobalElements64();
      nx = (int)pow((double)n, 0.33334);
      ny = nx; nz = nx;
      if (((long long) nx) * ny * nz != n)
        throw(Exception(__FILE__, __LINE__,
                        "You need to specify nx and ny"));
    }


    return(Matrices::Cross3D(Map, nx, ny, nz, 6.0, -1.0, -1.0,
                             -1.0, -1.0, -1.0, -1.0));
  }
  else
  {
    throw(Exception(__FILE__, __LINE__,
                    "`MatrixType' has incorrect value (" + MatrixType + ")",
                    "in input to function CreateMatrix()",
                    "Check the documentation for a list of valid choices"));
  }
} // CreateMatrix()

} // namespace Galeri
