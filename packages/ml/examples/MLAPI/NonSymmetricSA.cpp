
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "ml_config.h"
#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI.h"

using namespace Teuchos;
using namespace MLAPI;
//
// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  // Initialize the workspace and set the output level
  Init();

  try {

    int NX = 100;       // number of nodes on the X-axis
    int NY = 100;       // number of nodes on the Y-axis
    double LX = 1.0;    // length of the X-axis
    double LY = 1.0;    // lenght of the Y-axis
    double conv = 1e-5; // convection coefficient
    double diff = 1.0;  // diffusion coefficient

    Space FineSpace(NX * NY);

    DistributedMatrix MatA(FineSpace, FineSpace);

    if (GetMyPID() == 0) {

      double HX = LX / (NX + 1);  
      double HY = LY / (NY + 1);

      for (int ix = 0 ; ix < NX ; ++ix) {
        for( int iy = 0 ; iy < NY ; ++iy) {

          double X = HX * (ix + 1);
          double Y = HY * (iy + 1);

          double ConvX =  conv * 4 * X * (X - 1.) *(1. - 2 * Y) / HX;
          double ConvY = -conv * 4 * Y * (Y - 1.) *(1. - 2 * X) / HY;

          double lower = 0.0, upper = 0.0;
          double left = 0.0,  right = 0.0, center = 0.0;

          if( ConvX<0 ) {
            right  += ConvX;
            center -= ConvX;
          } else {
            left   -= ConvX;
            center += ConvX;
          }

          if( ConvY<0 ) {
            upper  += ConvY;
            center -= ConvY;
          } else {
            lower  -= ConvY;
            center += ConvY;
          }

          center += diff * 2. / (HX * HX) + diff * 2. / (HY * HY);
          left   -= diff / (HX * HX);
          right  -= diff / (HX * HX);
          lower  -= diff / (HY * HY);
          upper  -= diff / (HY * HY);

          int row = iy * NX + ix;

          if (ix != 0)
            MatA.SetElement(row, row - 1, left);

          if (ix != NX - 1)
            MatA.SetElement(row, row + 1, right);

          if (iy != 0)
            MatA.SetElement(row, row - NX, lower);

          if (iy != NY - 1)
            MatA.SetElement(row, row + NX, upper);

          MatA.SetElement(row, row, center);
        }
      }
    }

    MatA.FillComplete();

    // wrap MatA as an Operator
    Operator A(FineSpace, FineSpace, &MatA, false);

    Teuchos::ParameterList List;
    List.set("smoother: type", "symmetric Gauss-Seidel");
    List.set("smoother: sweeps", 1);
    List.set("smoother: damping factor", 1.0);
    List.set("coarse: max size", 32);

    MultiLevelAdaptiveSA Prec(A, List, true);

    // test the solver
    MultiVector LHS(FineSpace);
    MultiVector RHS(FineSpace);

    LHS.Random();
    RHS = 0.0;

    List.set("krylov: type", "GMRES");
    Krylov(A, LHS, RHS, Prec, List);

    Finalize(); 

  }
  catch (const int e) {
    cerr << "Caught integer exception, code = " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(0);

}

#else

#include "ml_include.h"

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("The ML API requires the following configuration options:");
  puts("\t--enable-epetra");
  puts("\t--enable-teuchos");
  puts("\t--enable-ifpack");
  puts("\t--enable-amesos");
  puts("Please check your configure line.");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
}
#endif // #if defined(HAVE_ML_MLAPI)
