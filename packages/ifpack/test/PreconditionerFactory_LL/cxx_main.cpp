/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

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
  Teuchos::RefCountPtr<Epetra_Map> Map = Teuchos::rcp( Galeri::CreateMap64("Linear", Comm, GaleriList) );
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

  Prec = Teuchos::rcp( Factory.Create("Polynomial", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  Prec = Teuchos::rcp( Factory.Create("Krylov", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;

  
  


#if defined (HAVE_IFPACK_SUPPORTGRAPH) && defined (HAVE_IFPACK_AMESOS)
  Prec = Teuchos::rcp( Factory.Create("MSF Amesos", &*A) );
  assert (Prec !=Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;


#endif



#ifdef HAVE_IFPACK_SUPPORTGRAPH
  Prec = Teuchos::rcp( Factory.Create("MSF IC", &*A) );
  assert (Prec != Teuchos::null);
  IFPACK_CHK_ERR(Prec->Initialize());
  IFPACK_CHK_ERR(Prec->Compute());
  cout << *Prec;


#endif

  if (Comm.MyPID() == 0)
    cout << "Test `PrecondititonerFactory_LL.exe' passed!" << endl;



#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

  return(EXIT_SUCCESS);
}
