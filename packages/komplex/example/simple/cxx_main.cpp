// @HEADER
// ************************************************************************
//
//           Komplex: An Equivalent Real Solver for Complex Linear Systems
//                 Copyright (2008) Sandia Corporation
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
//
// Questions about Komplex? Contact Michael Heroux (maherou at sandia.gov)
//
// ************************************************************************
// @HEADER
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "Komplex_LinearProblem.h"
#include "AztecOO.h"

using namespace Galeri;

void computeRHS(double c0r, double c0i, const Epetra_CrsMatrix & A0,
		double c1r, double c1i, const Epetra_CrsMatrix & A1,
		const Epetra_Vector & XrExact, const Epetra_Vector & XiExact, 
		Epetra_Vector & Br, Epetra_Vector & Bi) {

  // Computing 
  // Br = Ar*xr - Ai*xi
  // Bi = Ai*xr + Ar*xi
  // where 
  // Ar = c0r*A0 + c1r*A1
  // Ai = c0i*A0 + c1i*A1
  // Thus the full computation is
  // Br = (c0r*A0 + c1r*A1)*xr - (c0i*A0 + c1i*A1)*xi
  // Bi = (c0i*A0 + c1i*A1)*xr + (c0r*A0 + c1r*A1)*xi

  // Compute A0/1 * xr/i
  Epetra_Vector A0xr(Br), A0xi(Br), A1xr(Br), A1xi(Br);
  A0.Multiply(false,XrExact,A0xr);
  A0.Multiply(false,XiExact,A0xi);
  A1.Multiply(false,XrExact,A1xr);
  A1.Multiply(false,XiExact,A1xi);
  Br.Update(c0r, A0xr, c1r, A1xr, 0.0);
  Br.Update(-c0i, A0xi, -c1i, A1xi, 1.0);  
  Bi.Update(c0i, A0xr, c1i, A1xr, 0.0);
  Bi.Update(c0r, A0xi, c1r, A1xi, 1.0);

  return;
}
double computeResidualError(const Epetra_Vector & XrExact, const Epetra_Vector & XiExact, 
			    const Epetra_Vector & XrComp, const Epetra_Vector & XiComp) {
  Epetra_Vector vecr(XrExact), veci(XrExact);
  double resr, resi;
  vecr.Update(1.0, XrExact, -1.0, XrComp, 0.0);
  veci.Update(1.0, XiExact, -1.0, XiComp, 0.0);
  vecr.Dot(vecr, &resr);
  veci.Dot(veci, &resi);
  double result = resr+resi;
  return(result);

}

// =========== //
// main driver //
// =========== //

int main(int argv, char* argc[])
{
#ifdef HAVE_MPI
  MPI_Init(&argv, &argc);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // Here we create the linear problem
  //
  //   Matrix * LHS = RHS
  //
  // with Matrix arising from a 5-point formula discretization.
  
  Epetra_Map*         Map = 0;
  Epetra_CrsMatrix   *A0 = 0, *A1 = 0;

  Teuchos::ParameterList GaleriList;
  // dimension of the problem is nx x ny
  GaleriList.set("nx", 10 * Comm.NumProc());
  GaleriList.set("ny", 10);
  // total number of processors is mx x my
  GaleriList.set("mx", Comm.NumProc());
  GaleriList.set("my", 1);

  Map = CreateMap("Cartesian2D", Comm, GaleriList);
  A0 = CreateCrsMatrix("Biharmonic2D", Map, GaleriList);
  A1 = CreateCrsMatrix("Laplace2D", Map, GaleriList);
  Epetra_Vector XrExact(*Map); XrExact.PutScalar(1.0);
  Epetra_Vector XiExact(*Map); XiExact.PutScalar(1.0);
  Epetra_Vector Xr(*Map); Xr.PutScalar(0.0);
  Epetra_Vector Xi(*Map); Xi.PutScalar(0.0);
  Epetra_Vector Br(*Map);
  Epetra_Vector Bi(*Map);
  double c0r = 1.0, c0i = 0.0001;
  double c1r = 0.000001, c1i = 0.1;
  computeRHS(c0r, c0i, *A0, c1r, c1i, *A1, XrExact, XiExact, Br, Bi);
  Komplex_LinearProblem KP(c0r, c0i, *A0, c1r, c1i, *A1, Xr, Xi, Br, Bi);

  for (int i = 0; i<4; ++i) {

    switch(i) {
    case 0:
      break;
    case 1:
      c0r = 10.0; c0i = 1.0;
      c1r = 0.0; c1i = 1.0;
      break;
    case 2:
      c0r = 1.0; c0i = 0.0;
      c1r = 0.1; c1i = 0.0;
      break;
    case 3:
      c0r = 5.0; c0i = 0.1;
      c1r = 0.1; c1i = 1.0;
      break;
    default:
      exit(1);
    }
    
    if (i>0) {
      Xr.PutScalar(0.0);
      Xi.PutScalar(0.0);
      computeRHS(c0r, c0i, *A0, c1r, c1i, *A1, XrExact, XiExact, Br, Bi);
      KP.UpdateValues(c0r, c0i, *A0, c1r, c1i, *A1, Xr, Xi, Br, Bi);
    }

    bool debug = false;
    if (debug) {
      cout << "c0r = " << c0r << "   c0i = " << c0i << "   c1r = " << c1r << "   c1i = " << c1i << endl << endl;
      cout << "A0 = " << *A0 << endl;
      cout << "A1 = " << *A1 << endl;
      
      cout << "XrExact = " << XrExact << endl;
      cout << "XiExact = " << XiExact << endl;
      cout << "Br = " << Br << endl;
      cout << "Bi = " << Bi << endl;
    }
    // Get pointer to linear problem
    Epetra_LinearProblem * LP = KP.KomplexProblem();
    AztecOO solver(*LP);
    
    
    solver.Iterate(1000,1.0e-6);
    
    KP.ExtractSolution(Xr, Xi);
    
    // We compute sum of squares between exact and computed. 
    double residualNorm = computeResidualError(Xr, Xi, XrExact, XiExact);
    
    if (Comm.MyPID() == 0)
      cout << residualNorm << endl;
  }
  delete Map;
  delete A0;
  delete A1;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}
