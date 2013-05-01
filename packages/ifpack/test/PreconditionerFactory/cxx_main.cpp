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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "Ifpack_ConfigDefs.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack.h"
#include "AztecOO.h"
#ifdef HAVE_IFPACK_AMESOS
#include "Amesos_TestRowMatrix.h"
#endif

// =======================================================================
// GOAL: test that the names in the factory do not change. This test
//       will not solve any linear system.
//
int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;
  const int n = 9; 
  GaleriList.set("n", n);
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap("Linear", Comm, GaleriList) );
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A = Teuchos::rcp( Galeri::CreateCrsMatrix("Minij", &*Map, GaleriList) );
  
  Ifpack Factory;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Prec;

  Prec = Teuchos::rcp( Factory.Create("point relaxation", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("point relaxation stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("block relaxation", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("block relaxation stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("IC", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ICT", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ILU", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ILUT", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("IC stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ICT stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ILU stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("ILUT stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

#ifdef HAVE_IFPACK_AMESOS
  Prec = Teuchos::rcp( Factory.Create("Amesos", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("Amesos stand-alone", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;
#endif
  
  Prec = Teuchos::rcp( Factory.Create("Chebyshev", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("Krylov", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  if (Comm.MyPID() == 0)
    cout << "Test `PrecondititonerFactory.exe' passed!" << endl;

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
