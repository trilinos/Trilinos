// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_TEUCHOS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack_AdditiveSchwarz.h"
#include "AztecOO.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_gIct.h"
#include "Ifpack_vIct.h"
#include "Ifpack_gRiluk.h"
#include "Ifpack_vRiluk.h"

using namespace Trilinos_Util;

int TestWithIFPACK(Epetra_LinearProblem* Problem, const string what)
{

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  Teuchos::ParameterList List;
  Epetra_RowMatrix* A = Problem->GetMatrix();

  Ifpack_Preconditioner* Prec = 0;
  
  if (what == "Jacobi(1)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(A);
    List.set("point: type", "Jacobi");
  }
  if (what == "Jacobi(5)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation>(A);
    List.set("point: sweeps", 5);
  }
  // Cholesky
  else if (what == "vICT(0)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_vIct>(A);
    List.set("fact: level-of-fill", 0);
    List.set("schwarz: use reordering", true);
  }
  else if (what == "vICT(4)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_vIct>(A);
    List.set("fact: level-of-fill", 4);
    List.set("schwarz: use reordering", true);
  }
  else if (what == "gICT(0)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_gIct>(A);
    List.set("fact: level-of-fill", 0);
    List.set("schwarz: use reordering", true);
  }
  else if (what == "gICT(4)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_gIct>(A);
    List.set("fact: level-of-fill", 4);
    List.set("schwarz: use reordering", true);
  }
  // classical ILU
  else if (what == "vILU(0)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_vRiluk>(A);
    List.set("fact: level-of-fill", 0);
    List.set("schwarz: use reordering", true);
  }
  else if (what == "vILU(4)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_vRiluk>(A);
    List.set("fact: level-of-fill", 4);
    List.set("schwarz: use reordering", true);
//    List.set("schwarz: reordering type", "metis");
  }
  else if (what == "gILU(0)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_gRiluk>(A);
    List.set("fact: level-of-fill", 0);
    List.set("schwarz: use reordering", true);
  }
  else if (what == "gILU(4)") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_gRiluk>(A);
    List.set("fact: level-of-fill", 4);
    List.set("schwarz: use reordering", true);
  }
  // without reordering
  else if (what == "vILU(0) no reorder") {
    Prec = new Ifpack_AdditiveSchwarz<Ifpack_vRiluk>(A);
    List.set("fact: level-of-fill", 0);
    List.set("schwarz: use reordering", false);
  }

  assert(Prec != 0);
 

  IFPACK_CHK_ERR(Prec->SetParameters(List));
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());

  // create the AztecOO solver
  AztecOO AztecOOSolver(*Problem);

  // specify solver
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_cg);
  AztecOOSolver.SetAztecOption(AZ_output,32);

  AztecOOSolver.SetPrecOperator(Prec);

  AztecOOSolver.Iterate(150,1e-8);

  Prec->Print(cout);
  delete Prec;

  return(AztecOOSolver.NumIters());

}

// ====================================================================== 
int TestWithAztecOO(Epetra_LinearProblem* Problem, const string what)
{

  Epetra_MultiVector& RHS = *(Problem->GetRHS());
  Epetra_MultiVector& LHS = *(Problem->GetLHS());
  LHS.PutScalar(0.0);

  Teuchos::ParameterList List;
  Epetra_RowMatrix* A = Problem->GetMatrix();

  // create the AztecOO solver
  AztecOO AztecOOSolver(*Problem);

  // specify solver
  AztecOOSolver.SetAztecOption(AZ_solver,AZ_gmres);
  AztecOOSolver.SetAztecOption(AZ_overlap,0);
  AztecOOSolver.SetAztecOption(AZ_output,32);

  if (what == "Jacobi(1)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_Jacobi);
  }
  else if (what == "Jacobi(5)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_Jacobi);
    AztecOOSolver.SetAztecOption(AZ_poly_ord,5);
  }
  else if (what == "vICT(0)" || what == "gICT(0)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_icc);
    AztecOOSolver.SetAztecOption(AZ_graph_fill,0);
  }
  else if (what == "vICT(4)" || what == "gICT(4)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_icc);
    AztecOOSolver.SetAztecOption(AZ_graph_fill,4);
  }
  else if (what == "vILU(0)" || what == "gILU(0)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    AztecOOSolver.SetAztecOption(AZ_graph_fill,0);
  }
  else if (what == "vILU(4)" || what == "gILU(4)") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    AztecOOSolver.SetAztecOption(AZ_graph_fill,4);
  }
  else if (what == "vILU(0) no reorder") {
    AztecOOSolver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    AztecOOSolver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    AztecOOSolver.SetAztecOption(AZ_graph_fill,0);
    AztecOOSolver.SetAztecOption(AZ_reorder,0);
  }
  else
    exit(EXIT_FAILURE);

  AztecOOSolver.Iterate(1550,1e-8);

  return(AztecOOSolver.NumIters());

}
// ====================================================================== 
int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_Time Time(Comm);

  // size of the global matrix. 
  const int NumPoints = 27000;

  CrsMatrixGallery Gallery("laplace_3d", Comm);
  Gallery.Set("problem_size", NumPoints);
  Gallery.Set("map_type", "linear");

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();

  int TestPassed = true;

  vector<string> Tests;
  //Tests.push_back("Jacobi(5)");
  //Tests.push_back("Jacobi(5)");
  //Tests.push_back("gILU(4)");
  //Tests.push_back("vILU(0) no reorder");
  Tests.push_back("vILU(0)");
  //Tests.push_back("vILU(4)");
  //Tests.push_back("gICT(0)");

  for (int i = 0 ; i < Tests.size() ; ++i) {
    int AztecOOIters, IFPACKIters;
    double AztecOOTime, IFPACKTime;

    Time.ResetStartTime();

//    AztecOOIters =  TestWithAztecOO(Problem,Tests[i]);
    AztecOOTime = Time.ElapsedTime();

    Time.ResetStartTime();
    IFPACKIters = TestWithIFPACK(Problem,Tests[i]);
    IFPACKTime = Time.ElapsedTime();

    cout << endl;
    cout << "Testing `" << Tests[i] << "'" << endl;
    cout << "Total time for AztecOO = " << AztecOOTime << " (s)" << endl;
    cout << "Total time for IFPACK  = " << IFPACKTime << " (s)" << endl;
    cout << "Iterations required by AztecOO = " << AztecOOIters << endl;
    cout << "Iterations required by IFPACK  = " << IFPACKIters << endl;

    if (AztecOOIters != IFPACKIters)
      TestPassed = false;
  }

  if (!TestPassed) {
    cerr << "TEST FAILED!!!!!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif
  cout << "TEST PASSED" << endl;
  exit(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
  puts("to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
