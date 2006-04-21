
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
/*#############################################################################
# CVS File Information
#    Current revision: $Revision$
#    Branch:           $Branch$
#    Last modified:    $Date$
#    Modified by:      $Author$
#############################################################################*/

// Sample driver for Maxwell equation AMG solver in the ML package. The
// software is tested by setting up a 2-dimensional uniform grid example on 
// a square. For details about the problem at hand, please refer to file
// ml_simple_max.c, of which this file is the C++ counterpart.
//
// This file shows how to use the class
// ML_Epetra::MultiLevelPreconditioner to solve this formulation of the
// Maxwell equations. The class takes care of building the node and edge
// hierarchy, definining the Hiptmair smoother, and setting the coarse
// solver. More information about MultiLevelPreconditioner can be found in
// the ML User's Guide.

#include "ml_include.h"

// ML_Epetra::MultiLevelPreconditioner() requires Trilinos to be
// configured with --enable-epetra --enable-teuchos. This example
// requires --enable-triutils for the definition of the linear systems.

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util.h"

int main(int argc, char *argv[])
{

#ifdef ML_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  char *datafile;

  if (argc != 5) {
    if (Comm.MyPID() == 0) {
      cout << "usage: ml_maxwell.exe <S> <M> <T> <Kn>" << endl;
      cout << "        S = edge stiffness matrix file" << endl;
      cout << "        M = edge mass matrix file" << endl;
      cout << "        T = discrete gradient file" << endl;
      cout << "       Kn = auxiliary nodal FE matrix file" << endl;
    }
#ifdef ML_MPI
    MPI_Finalize();
#endif
    exit(1);
  }

  Epetra_Map *map;
  Epetra_Map *edgemap, *nodemap;
  Epetra_CrsMatrix *CurlCurl, *Mass, *T, *Kn;
  Epetra_Vector *epx;
  Epetra_Vector *epb;
  Epetra_Vector *xexact;

  // ==================================================== //
  // READ IN PROBLEM FROM FILE                            //
  // ==================================================== //

  //curlcurl matrix
  datafile = argv[1];
  cout << "Reading in " << *datafile << " ...";
  Trilinos_Util_ReadMatrixMarket2Epetra(datafile, Comm, map,
                      CurlCurl,
                      epx, epb, xexact );
  cout << "   "
       << CurlCurl->NumMyRows() << " rows, "
       << CurlCurl->NumMyCols() << " columns" << endl;
  delete map; delete epx; delete epb; delete xexact;

  //mass matrix
  datafile = argv[2];
  cout << "Reading in " << *datafile << " ...";
  Trilinos_Util_ReadMatrixMarket2Epetra(datafile, Comm, map,
                      Mass,
                      epx, epb, xexact );
  cout << "   "
       << Mass->NumMyRows() << " rows, "
       << Mass->NumMyCols() << " columns" << endl;
  delete map; delete epx; delete epb; delete xexact;

  //T matrix
  datafile = argv[3];
  cout << "Reading in " << *datafile << " ...";
  Trilinos_Util_ReadMatrixMarket2Epetra(datafile, Comm, nodemap,
                      T,
                      epx, epb, xexact );
  cout << "   "
       << T->NumMyRows() << " rows, "
       << T->NumMyCols() << " columns" << endl;
  delete epx; delete epb; delete xexact;

  //node matrix
  datafile = argv[4];
  cout << "Reading in " << datafile << " ...";
  Trilinos_Util_ReadMatrixMarket2Epetra(datafile, Comm, map,
                      Kn,
                      epx, epb, xexact );
  cout << "   "
       << Kn->NumMyRows() << " rows, "
       << Kn->NumMyCols() << " columns" << endl;
  delete map; delete epx; delete epb; delete xexact;

  // ==================================================== //
  // S E T U P   O F    M L   P R E C O N D I T I O N E R //
  // ==================================================== //

  Teuchos::ParameterList MLList;
  int *options    = new int[AZ_OPTIONS_SIZE];
  double *params  = new double[AZ_PARAMS_SIZE];
  ML_Epetra::SetDefaults("maxwell", MLList, options, params);

  // how verbose ML is
  MLList.set("output", 10);

  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("coarse: max size", 30);
  MLList.set("aggregation: threshold", 0.0);
  MLList.set("max levels",2);
  //MLList.set("R and P smoothing: damping", "non-smoothed");

  // coarse level solve
  //MLList.set("coarse: type", "Amesos-KLU");
  MLList.set("coarse: type", "Hiptmair");

  // this creates the multilevel hierarchy, sets the smoothers,
  // prepares the coarse solver, etc...

  ML_Epetra::MultiLevelPreconditioner * MLPrec =
    new ML_Epetra::MultiLevelPreconditioner(*CurlCurl, *Mass, *T, *Kn, MLList);
    //new ML_Epetra::MultiLevelPreconditioner(*CurlCurl, *T, *Kn, MLList);

  MLPrec->PrintUnused(0);

  // ========================================================= //
  // D E F I N I T I O N   O F   A Z T E C O O   P R O B L E M //
  // ========================================================= //

  // create left-hand side and right-hand side, and populate them with
  // data from file. Both vectors are defined on the domain map of the
  // edge matrix.
  // Epetra_Vectors can be created in View mode, to accept pointers to
  // double vectors.

  Epetra_Vector rhs(CurlCurl->DomainMap());
  rhs.Random();
  Epetra_Vector x(CurlCurl->DomainMap());
  x.PutScalar(0.0);

  // for AztecOO, we need an Epetra_LinearProblem
  Epetra_CrsMatrix *Combined = Epetra_MatrixAdd(CurlCurl,Mass,1.0);
  Epetra_LinearProblem Problem(Combined,&x,&rhs);
  // AztecOO Linear problem
  AztecOO solver(Problem);
  // set MLPrec as precondititoning operator for AztecOO linear problem
  solver.SetPrecOperator(MLPrec);

  // a few options for AztecOO
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 1);

  solver.Iterate(150, 1e-12);

  // =============== //
  // C L E A N   U P //
  // =============== //

  delete edgemap;
  delete nodemap;
  delete MLPrec;    // destroy phase prints out some information

#ifdef ML_MPI
  MPI_Finalize();
#endif
        
  return 0;
        
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
#if !defined(HAVE_ML_EPETRA)
  puts("--enable-epetra");
#endif
#if !defined(HAVE_ML_TEUCHOS)
  puts("--enable-teuchos");
#endif
#if !defined(HAVE_ML_TRIUTILS)
  puts("--enable-triutils");
#endif
#if !defined(HAVE_ML_AZTECOO)
  puts("--enable-aztecoo");
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}

#endif
