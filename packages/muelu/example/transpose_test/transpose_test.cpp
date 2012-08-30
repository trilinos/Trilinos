// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * transpose_test.cpp
 *
 *  Created on: 22.09.2011
 *      Author: tobias
 */

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Version.h"

#include "Epetra_CrsMatrix.h"

#ifdef HAVE_MUELU_EPETRAEXT
#include "EpetraExt_CrsMatrixIn.h"
#include <EpetraExt_MatrixMatrix.h>
#endif

#include "Epetra_RowMatrixTransposer.h"

#include <Teuchos_RCP.hpp>

#include "Xpetra_EpetraMap.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_EpetraVector.hpp"

#include "MueLu_Utilities.hpp"

Epetra_CrsMatrix* TridiagMatrix(const Epetra_Map* Map, const double a, const double b, const double c)
{
  Epetra_CrsMatrix* Matrix = new Epetra_CrsMatrix(Copy, *Map,  3);

  int NumGlobalElements = Map->NumGlobalElements();
  int NumMyElements = Map->NumMyElements();
  int* MyGlobalElements = Map->MyGlobalElements();

  std::vector<double> Values(2);
  std::vector<int> Indices(2);
  int NumEntries;

  for (int i = 0; i < NumMyElements; ++i)
  {
    if (MyGlobalElements[i] == 0)
    {
      // off-diagonal for first row
      Indices[0] = 1;
      NumEntries = 1;
      Values[0] = c;
    }
    else if (MyGlobalElements[i] == NumGlobalElements - 1)
    {
      // off-diagonal for last row
      Indices[0] = NumGlobalElements - 2;
      NumEntries = 1;
      Values[0] = b;
    }
    else
    {
      // off-diagonal for internal row
      Indices[0] = MyGlobalElements[i] - 1;
      Values[1] = b;
      Indices[1] = MyGlobalElements[i] + 1;
      Values[0] = c;
      NumEntries = 2;
    }
    Matrix->InsertGlobalValues(MyGlobalElements[i], NumEntries, &Values[0],
                               &Indices[0]);
    // Put in the diagonal entry
    Values[0] = a;
    Matrix->InsertGlobalValues(MyGlobalElements[i], 1, &Values[0],
                               MyGlobalElements + i);
  }

  // Finish up, trasforming the matrix entries into local numbering,
  // to optimize data transfert during matrix-vector products
  Matrix->FillComplete();
  Matrix->OptimizeStorage();

  return(Matrix);
}

int main(int argc, char *argv[])
{

  std::cout << Epetra_Version() << std::endl << std::endl;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif


  Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));

  // Construct a Map with NumElements and index base of 0
  Teuchos::RCP<Epetra_Map> rgMap = Teuchos::rcp(new Epetra_Map(63, 0, Comm));
  Teuchos::RCP<Epetra_Map> doMap = Teuchos::rcp(new Epetra_Map(20, 0, Comm));

  Epetra_CrsMatrix* Pin = NULL;
  EpetraExt::MatrixMarketFileToCrsMatrix("Test.mat",
      *rgMap, *rgMap, *doMap,
      Pin, false,
      true);
  Teuchos::RCP<Epetra_CrsMatrix> P = Teuchos::rcp(Pin);

  Epetra_CrsMatrix* A = NULL;
  A = TridiagMatrix(rgMap.get(), 1.0, -2.0, 1.0);
  Teuchos::RCP<Epetra_CrsMatrix> rcpA = Teuchos::rcp(A);

  ////////////////////////////////////////////
  // plain Epetra
  ////////////////////////////////////////////

  // Create x and b vectors
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*doMap));
  Teuchos::RCP<Epetra_Vector> x2 = Teuchos::rcp(new Epetra_Vector(*rgMap));
  Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(*rgMap));
  Teuchos::RCP<Epetra_Vector> b2 = Teuchos::rcp(new Epetra_Vector(*rgMap));

  x->PutScalar(1.0);
  x2->PutScalar(1.0);
  double normx = 0.0;
  x->Norm1(&normx);
  if (Comm.MyPID() == 0) std::cout << "||x|| = " << normx << std::endl;
  x2->Norm1(&normx);
  if (Comm.MyPID() == 0) std::cout << "||x2|| = " << normx << std::endl;

  /*P->Apply(*x,*b);
  normx = 0.0;
  b->Norm1(&normx);*/
  //if (Comm.MyPID() == 0) std::cout << "||Px||_1 = " << normx << std::endl;

  /*Epetra_RowMatrixTransposer et(&*P);
        Epetra_CrsMatrix* PT;
        int rv = et.CreateTranspose(true,PT);
        if (rv != 0) {
          std::ostringstream buf;
          buf << rv;
          std::string msg = "Utils::Transpose: Epetra::RowMatrixTransposer returned value of " + buf.str();
          std::cout << msg << std::endl;
        }

  Teuchos::RCP<Epetra_CrsMatrix> rcpPT(PT);

  rcpPT->Apply(*x2,*b2);
  normx = 0.0;
  b2->Norm1(&normx);*/
  //if (Comm.MyPID() == 0) std::cout << "||P^T x||_1 = " << normx << std::endl;

  // matrix-matrix multiply
  Teuchos::RCP<Epetra_CrsMatrix> AP = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rcpA->RangeMap(),1));
  EpetraExt::MatrixMatrix::Multiply(*rcpA,false,*P,false,*AP,false);
  AP->FillComplete(P->DomainMap(),rcpA->RangeMap());
  //std::cout << *AP << std::endl;
  AP->Apply(*x,*b2);
  normx = 0.0;
  b2->Norm1(&normx);
  if (Comm.MyPID() == 0) std::cout << "Epetra: ||AP x||_1 = " << normx << std::endl;

  // build A^T explicitely
  Epetra_RowMatrixTransposer etA(&*rcpA);
        Epetra_CrsMatrix* AT;
        int rv3 = etA.CreateTranspose(true,AT);
        if (rv3 != 0) {
          std::ostringstream buf;
          buf << rv3;
          std::string msg = "Utils::Transpose: Epetra::RowMatrixTransposer returned value of " + buf.str();
          std::cout << msg << std::endl;
        }
  Teuchos::RCP<Epetra_CrsMatrix> rcpAT(AT);

  // calculate A^T Px
  Teuchos::RCP<Epetra_CrsMatrix> APexpl = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rcpA->DomainMap(),20));
  EpetraExt::MatrixMatrix::Multiply(*rcpAT,false,*P,false,*APexpl,false);
  APexpl->FillComplete(P->DomainMap(),rcpA->DomainMap()); // check me
  APexpl->Apply(*x,*b2);
  normx = 0.0;
  b2->Norm1(&normx);
  if (Comm.MyPID() == 0) std::cout << "Epetra: ||A^T_expl P x||_1 = " << normx << std::endl;

  // calculate A^T Px
  Teuchos::RCP<Epetra_CrsMatrix> APimpl = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rcpA->RangeMap(),20));
  EpetraExt::MatrixMatrix::Multiply(*rcpA,true,*P,false,*APimpl,false);
  APimpl->FillComplete(P->DomainMap(),APimpl->RangeMap()); // check me
  APimpl->Apply(*x,*b2);
  normx = 0.0;
  b2->Norm1(&normx);
  if (Comm.MyPID() == 0) std::cout << "Epetra: ||A^T_impl P x||_1 = " << normx << std::endl;


  ///////////////////////////////////////
  // Xpetra
  ///////////////////////////////////////


   // wrap original matrix to Matrix
   Teuchos::RCP<Xpetra::EpetraMap> xrgMap = Teuchos::rcp(new Xpetra::EpetraMap(rgMap));
   Teuchos::RCP<Xpetra::EpetraMap> xdoMap = Teuchos::rcp(new Xpetra::EpetraMap(doMap));

   Teuchos::RCP<Xpetra::EpetraCrsMatrix> Pop = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(P) );
   Teuchos::RCP<Xpetra::CrsMatrix<double> > crsP = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<double> >(Pop);
   Teuchos::RCP<Xpetra::CrsMatrixWrap<double> > crsPOp = Teuchos::rcp( new Xpetra::CrsMatrixWrap<double>(crsP) );
   crsPOp->fillComplete(xdoMap,xrgMap);

   Teuchos::RCP<Xpetra::EpetraCrsMatrix> Aop = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(rcpA) );
   Teuchos::RCP<Xpetra::CrsMatrix<double> > crsA = Teuchos::rcp_implicit_cast<Xpetra::CrsMatrix<double> >(Aop);
   Teuchos::RCP<Xpetra::CrsMatrixWrap<double> > crsAOp = Teuchos::rcp( new Xpetra::CrsMatrixWrap<double>(crsA) );
   crsAOp->fillComplete(xrgMap,xrgMap);

   // wrap vectors
   Teuchos::RCP<Xpetra::EpetraVector> xx = Teuchos::rcp(new Xpetra::EpetraVector(x));
   Teuchos::RCP<Xpetra::EpetraVector> xx2 = Teuchos::rcp(new Xpetra::EpetraVector(x2));
   Teuchos::RCP<Xpetra::EpetraVector> bb1 = Teuchos::rcp(new Xpetra::EpetraVector(b));
   Teuchos::RCP<Xpetra::EpetraVector> bb2 = Teuchos::rcp(new Xpetra::EpetraVector(b2));

   bb1->putScalar(0.0);
   bb2->putScalar(0.0);

   //crsPOp->apply(*xx,*bb1);

   // if (Comm.MyPID() == 0) std::cout << "||Pop x||_1 = " << bb1->norm1() << std::endl;

   //Teuchos::RCP<Xpetra::Matrix<double> > crsPOpT = MueLu::Utils2<double>::Transpose(crsPOp,false);

   //crsPOpT->apply(*xx2,*bb2);
   //if (Comm.MyPID() == 0) std::cout << "||Pop^T x||_1 = " << bb2->norm1() << std::endl;

   // calculate APx
   Teuchos::RCP<Xpetra::Matrix<double> > crsAP = MueLu::Utils<double>::TwoMatrixMultiply(crsAOp,false,crsPOp,false);
   //crsAP->describe(*fos,Teuchos::VERB_EXTREME);
   bb1->putScalar(0.0);
   crsAP->apply(*xx,*bb1);
   normx = bb1->norm1();
   if (Comm.MyPID() == 0) std::cout << "Xpetra: ||A Pop x||_1 = " << normx << std::endl;

   // calculate A^T P x explicitely
   Teuchos::RCP<Xpetra::Matrix<double> > crsAOpT = MueLu::Utils2<double>::Transpose(crsAOp,false);
   Teuchos::RCP<Xpetra::Matrix<double> > AtPexpl = MueLu::Utils<double>::TwoMatrixMultiply(crsAOpT,false,crsPOp,false);
   bb1->putScalar(0.0);
   AtPexpl->apply(*xx,*bb1);
   normx = bb1->norm1();
   if (Comm.MyPID() == 0)
       std::cout << "Xpetra: ||A^T_expl Pop x||_1 = " << normx << std::endl;

   // calculate A^T P x
   Teuchos::RCP<Xpetra::Matrix<double> > AtPimpl = MueLu::Utils<double>::TwoMatrixMultiply(crsAOp,true,crsPOp,false);
   bb1->putScalar(0.0);
   AtPimpl->apply(*xx,*bb1);
   normx = bb1->norm1();
   if (Comm.MyPID() == 0)
       std::cout << "Xpetra: ||A^T_impl Pop x||_1 = " << normx << std::endl;


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
