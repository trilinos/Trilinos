
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

// Shows the usage of ml_ParameterList.xml file reader.
// The user must create an XML file called ``ml_ParameterList.xml'', and
// locate it in the working directory. Every time an ML preconditioner is
// built, the code checks for the presence of this file; if found, the XML
// content is parsed into a Teuchos::ParameterList that by default is added to
// the already existing list.
// The following two additional parameters are parsed:
//
// - ResetList [bool]: if true, the new list read from file overwrites
//                     the already existing list; is false, the new list
//                     is added to the existing list.
//
// - SetDefaults [string]: if present, the parsed list is first set to
//                         default values.
//
// Please check file ml/examples/BasicExamples/ml_ParameterList.xml
// for an example of input format. Note that the XML file must contain
// on object, called ParameterList, whose name is "MultiLevelPreconditioner".
// Also consult Teuchos' XML I/O capabilities.
//
// This example requires Teuchos to be compiled with --enable-teuchos-expat.
//
// WARNING: File ml_ParameterList.xml is now empty! You should edit this file;
// for example, if you want to use IFPACK smoothers for all levels except the
// finest, where you want symmetric Gauss-Seidel, the content of the file
// might be:
/*---(beginning of ml_ParameterList.xml)---
<ParameterList name="MultiLevelPreconditioner">
  <Parameter name="aggregation: damping factor" type="double" value="1.333"/>
  <Parameter name="smoother: type (level 0)" type="string" value="symmetric Gauss-Seidel"/>
  <Parameter name="smoother: type (level 1)" type="string" value="IFPACK"/>
  <Parameter name="smoother: ifpack type" type="string" value="ILU"/>
  <Parameter name="smoother: ifpack overlap" type="int" value="1"/>
  <ParameterList name="smoother: ifpack list">
    <Parameter name="fact: level-of-fill" type="int" value="1"/>
  </ParameterList>
</ParameterList>
---(end of ml_ParameterList.xml)---*/
//
// \author Marzio Sala, ETHZ/COLAB
//
// \data Last modified on 19-Mar-06.

#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) && defined(HAVE_ML_AZTECOO)

// epetra objects
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
// required to build the example matrix
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
// required by the linear system solver
#include "AztecOO.h"
// required by ML
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Galeri;

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx = 128;
  int ny = nx * Comm.NumProc(); // each subdomain is a square

  ParameterList GaleriList;
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);
  GaleriList.set("mx", 1);
  GaleriList.set("my", Comm.NumProc());

  Epetra_Map* Map = CreateMap("Cartesian2D", Comm, GaleriList);
  Epetra_CrsMatrix* A = CreateCrsMatrix("Laplace2D", Map, GaleriList);
    
  // Build a linear system with trivial solution, using a random vector
  // as starting solution.
  Epetra_Vector LHS(*Map); LHS.Random();
  Epetra_Vector RHS(*Map); RHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, &RHS);

  // As we wish to use AztecOO, we need to construct a solver object 
  // for this problem
  AztecOO solver(Problem);

  // =========================== begin of ML part ===========================
  
  // Create an almost empty list. Note: you *must* add
  //      MLList.set("read XML", true);
  // to your MLList in order to read the XML file.
  Teuchos::ParameterList MLList;
  MLList.set("ML output", 10);
  MLList.set("read XML", true);

  // Since a file called ml_ParameterList.xml is present in the working
  // directory, ML will parse this XML file and use it to specify the
  // parameters. 
  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  // =========================== end of ML part =============================
  
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  // destroy the preconditioner
  delete MLPrec;
  
  delete A;
  delete Map;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
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
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-galeri");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return(EXIT_SUCCESS);
}

#endif
