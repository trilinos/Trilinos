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
#ifndef __MATRIX_VECTOR_CHECKER_HPP__
#define __MATRIX_VECTOR_CHECKER_HPP__

#include <cstdlib>
#include <iostream>
#include <iomanip>

#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_MatrixMatrix.h>

int MatrixVectorChecker(const Epetra_RowMatrix & Mat) {
// Check the matrix-vector product, getrow, and matrix-matrix
// multiply associated with Mat by doing the following:
//
//      1) create random vector v where .5 < v_i < 1   or  -1 < v_i < -.5
//
//      2) y <--  Mat*v    via a standard matrix-vector product
//
//      3) w <--  Mat*v    via a matrix-vector product implemented with getrow
//
//      4) z <-- |Mat|*|v| via a matrix-vector product implemented with getrow
//
//      5) print if abs(y[i] - w[i])/z[i] > 1.e-8
//
// If Mat is a Epetra_CrsMatrix ...
//
//      6) tCRS <-- Mat*v via matrix-matrix multiply
//
//      7) print if  abs(y[i] - tCRS[i])/z[i] > 1.e-8

  const Epetra_Map & ColMap = Mat.RowMatrixColMap();
  const Epetra_Map & RowMap = Mat.RowMatrixRowMap();
  const Epetra_Map & DomMap = Mat.MatrixDomainMap();
  const Epetra_Map & RanMap = Mat.MatrixRangeMap();
  const int           MyPID = RowMap.Comm().MyPID();
  Epetra_Vector v(DomMap),y(RanMap), w(RowMap,true), z(RowMap,true), t(RowMap);
  v.Random();
  for (int i = 0; i < DomMap.NumMyElements(); i++) {
     if ((v[i] > -.5 ) && (v[i] < 0. )) v[i] -= .5;
     if ((v[i] <  .5 ) && (v[i] > 0. )) v[i] += .5;
  }

  // y <-- Mat*v via standard multiply
  Mat.Multiply(false,v,y);

  // Copy v and add ghost stuff imported from other processors

  double *vdata;
  Epetra_Vector vhat(ColMap);

  if (Mat.RowMatrixImporter()) {
     vhat.Import(v, *(Mat.RowMatrixImporter()), Insert);
     vdata = vhat.Values();
  }
  else vdata = v.Values();

  // w <-- A*v and z <-- |A||v| via a getrow multiply

  double *Values  = new double  [ColMap.NumMyElements()+1];
  int    *Indices = new int     [ColMap.NumMyElements()+1];
  int    RowLength;
  for (int i = 0; i < RowMap.NumMyElements(); i++) {
    Mat.ExtractMyRowCopy(i, ColMap.NumMyElements(), RowLength, Values, Indices);
    for (int j = 0; j < RowLength; j++) w[i] += (Values[j]*vdata[Indices[j]]);
    for (int j = 0; j < RowLength; j++) z[i] += fabs(Values[j]*vdata[Indices[j]]);
    if (z[i] == 0.) z[i] = 1.e-20;
  }

  // We need import/export stuff on a what/zhat somewhat similar to vhat above?
  if (!RanMap.PointSameAs(RowMap))
    std::cout << "Range and Rowmaps do not match some import/export code needs to be added" << std::endl;

  // print if abs(y[i] - w[i])/z[i] > 1.e-8

  int Count = 0;
  for (int i = 0; i < RowMap.NumMyElements(); i++) {
    if (( fabs(y[i] - w[i])/z[i] > 1.e-8) && (Count++ < 20) ) {
       std::cout << std::setw(5) << MyPID << ":matvec/getrow mismatch,row (GID=" <<
       std::setw(6) << RowMap.GID(i) << ",LID=" << std::setw(6) << i <<
       "), y w & z =" << std::setw(9) << std::scientific <<
       std::setprecision(2) << y[i] << " " << std::setw(9) << w[i]
       << " " << std::setw(9) << z[i] << std::endl;
    }
  }
  const Epetra_CrsMatrix *MatCRS =
              dynamic_cast<const Epetra_CrsMatrix *>(&Mat);

  // tCRS <-- Mat*v via matrix-matrix multiply if Mat is CRS
  if (MatCRS) {
     if (MyPID == 0) std::cout << "crs" << std::endl;
     Epetra_CrsMatrix vCRS( Copy, DomMap, 1);
     Epetra_CrsMatrix tCRS( Copy, RanMap, 1);

     // create vector vCRS which is matrix version of v.
     int *mygids = new int [DomMap.NumMyElements()+1];
     DomMap.MyGlobalElements(mygids);
     double *ptr;
     v.ExtractView(&ptr);
     int zero = 0;
     for (int i = 0; i < DomMap.NumMyElements(); i++)
        vCRS.InsertGlobalValues(mygids[i], 1, &(ptr[DomMap.LID(mygids[i])]), &zero);

     // create map with only element assigned to one processor.
     Epetra_Map       VecAsMatMap(1,0,RowMap.Comm());
     vCRS.FillComplete(VecAsMatMap,DomMap);
     vCRS.OptimizeStorage(); //JJH this is done by default now
     EpetraExt::MatrixMatrix::Multiply(*MatCRS,false,vCRS,false,tCRS);

     // now getrow tCRS and compare with y
     Count = 0;
     for (int i = 0; i < RowMap.NumMyElements(); i++) {
        tCRS.ExtractMyRowCopy(i, ColMap.NumMyElements(), RowLength, Values, Indices);
        if (RowLength == 0) {
           RowLength = 1; Values[0] = 0.; Indices[0] = 0;
         }
        if (RowLength != 1) {
          std::cout << std::setw(5) << MyPID << ":matvec/matmult row (GID=" <<
          std::setw(6) << RowMap.GID(i) << ",LID=" << std::setw(6) << i <<
          ") rowlength problem" << std::endl;
        }
        if (Indices[0] != 0) {
          std::cout << std::setw(5) << MyPID << ":matvec/matmult row (GID=" <<
          std::setw(6) << RowMap.GID(i) << ",LID=" << std::setw(6) << i <<
          ") has Col Id = " << Indices[0] << std::endl;
        }
        t[i] = Values[0];
     }
     // might need to do some import/export stuff if RangeMap and RowMap
     // do not match. This all depends on what matmat mult does?

     // print if abs(y[i] - t[i])/z[i] > 1.e-8

     for (int i = 0; i < RowMap.NumMyElements(); i++) {
        if (( fabs(y[i] - t[i])/z[i] > 1.e-8) && (Count++ < 20) ) {
          std::cout << std::setw(5) << MyPID << ":matvec/matmult mismatch,row (GID=" <<
          std::setw(6) << RowMap.GID(i) << ",LID=" << std::setw(6) << i <<
          "), y t & z =" << std::setw(9) << std::scientific <<
          std::setprecision(2) << y[i] << " " << std::setw(9) << Values[0]
          << " " << std::setw(9) << z[i] << std::endl;
        }
     } // for (int i = 0; ...
  } // if (MatCRS)
  else
    if (MyPID == 0) std::cout << "not crs" << std::endl;
  return 0;

} // MatrixVectorChecker()


#endif
