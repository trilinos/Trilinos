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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayViewDecl.hpp>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

namespace MueLu {

// Check the matrix-vector product, getrow, and matrix-matrix
// multiply associated with mat by doing the following:
//
//      1) create random vector v where .5 < v_i < 1   or  -1 < v_i < -.5
//
//      2) y <--  mat*v    via a standard matrix-vector product
//
//      3) w <--  mat*v    via a matrix-vector product implemented with getrow
//
//      4) z <-- |mat|*|v| via a matrix-vector product implemented with getrow
//
//      5) print if abs(y[i] - w[i])/z[i] > 1.e-8
//
// If mat is a Epetra_CrsMatrix ...
//
//      6) tCRS <-- mat*v via matrix-matrix multiply
//
//      7) print if  abs(y[i] - tCRS[i])/z[i] > 1.e-8
template <typename Scalar,typename LocalOrdinal,typename GlobalOrdinal,typename Node, typename LocalMatOps>
int MatrixVectorChecker(const RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > & mat) {
#include "MueLu_UseShortNames.hpp"

  // const RCP<const Map> & colMap = mat->getColMap(); // unused
  const Teuchos::RCP<const Map> & rowMap = mat->getRowMap();
  const Teuchos::RCP<const Map> & domMap = mat->getDomainMap();
  const Teuchos::RCP<const Map> & ranMap = mat->getRangeMap();
  const int myPID = rowMap->getComm()->getRank();

  // MyVector vv(domMap), yy(ranMap), ww(rowMap,true), zz(rowMap,true); 
  Teuchos::RCP<Vector> vv = VectorFactory::Build(domMap); //TODO: RCP ?
  Teuchos::RCP<Vector> yy = VectorFactory::Build(ranMap);
  Teuchos::RCP<Vector> ww = VectorFactory::Build(rowMap,true);
  Teuchos::RCP<Vector> zz = VectorFactory::Build(rowMap,true);

  Teuchos::ArrayRCP<SC> v = vv->getDataNonConst(0);
  Teuchos::ArrayRCP<SC> y = yy->getDataNonConst(0);
  Teuchos::ArrayRCP<SC> w = ww->getDataNonConst(0);
  Teuchos::ArrayRCP<SC> z = zz->getDataNonConst(0);
  
  vv->randomize();

  for (unsigned int i = 0; i < domMap->getNodeNumElements(); i++) {
     if ((v[i] > -.5 ) && (v[i] < 0. )) v[i] -= .5;
     if ((v[i] <  .5 ) && (v[i] > 0. )) v[i] += .5;
  }

  // y <-- mat*v via standard multiply
  mat->apply(*vv,*yy, Teuchos::NO_TRANS, 1, 0);

  // TODO ??
  // Copy v and add ghost stuff imported from other processors
  // MyVector vhat(colMap);  //TODO: replace MyVector
  //   if (mat->RowMatrixImporter()) {
  //      vhat.Import(v, *(mat->RowMatrixImporter()), Insert);
  //      v = vhat.Values();
  //   }
  //   else v = v.Values();

  // w <-- A*v and z <-- |A||v| via a getrow multiply

  Teuchos::ArrayView<const LO> indices;
  Teuchos::ArrayView<const SC> values;

  LO RowLength;
  for (size_t i = 0; i < rowMap->getNodeNumElements(); i++) {

    mat->getLocalRowView(i, indices, values);
    RowLength = indices.size();

    for (LO j = 0; j < RowLength; j++) {
      w[i] += (values[j]*v[indices[j]]);
      // ww->sumIntoLocalValue(i, 0, values[j]*v[indices[j]]); // Not implemented
    }

    for (LO j = 0; j < RowLength; j++) {
      z[i] += fabs(values[j]*v[indices[j]]);
      // zz->sumIntoLocalValue(i, 0, fabs(values[j]*v[indices[j]])); // Not implemented
    }
    if (z[i] == 0.) z[i] = 1.e-20;
  }

  // We need import/export stuff on a what/zhat somewhat similar to vhat above?
  if (!ranMap->isCompatible(*rowMap)) 
    std::cout << "Range and Rowmaps do not match. Some import/export code needs to be added" << std::endl;

  // print if abs(y[i] - w[i])/z[i] > 1.e-8

  int Count = 0;
  for (size_t i = 0; i < rowMap->getNodeNumElements(); i++) {
    if (( fabs(y[i] - w[i])/z[i] > 1.e-8) && (Count++ < 20) ) {
       std::cout << std::setw(5) << myPID << ":matvec/getrow mismatch,row (GID=" << 
       std::setw(6) << rowMap->getGlobalElement(i) << ",LID=" << std::setw(6) << i << 
       "), y w & z =" << std::setw(9) << std::scientific << 
       std::setprecision(2) << y[i] << " " << std::setw(9) << w[i] 
       << " " << std::setw(9) << z[i] << std::endl;
    }
  }

#ifdef TODO_TODO // TODO
  const Epetra_CrsMatrix *matCRS =
              dynamic_cast<const Epetra_CrsMatrix *>(&mat);

  // tCRS <-- mat*v via matrix-matrix multiply if mat is CRS
  if (matCRS) {
     if (myPID == 0) std::cout << "crs" << std::endl;
     Epetra_CrsMatrix vCRS( Copy, domMap, 1);
     Epetra_CrsMatrix tCRS( Copy, ranMap, 1);

     // create vector vCRS which is matrix version of v.
     int *mygids = new int [domMap.getNodeNumElements()+1];
     domMap.MyGlobalElements(mygids);
     double *ptr;   
     v.ExtractView(&ptr);
     int zero = 0;
     for (int i = 0; i < domMap.getNodeNumElements(); i++) 
        vCRS.InsertGlobalValues(mygids[i], 1, &(ptr[domMap.LID(mygids[i])]), &zero);

     // create map with only element assigned to one processor.
     Epetra_Map       VecAsMatMap(1,0,rowMap.Comm());
     vCRS.FillComplete(VecAsMatMap,domMap);
     vCRS.OptimizeStorage(); //JJH this is done by default now
     EpetraExt::MatrixMatrix::Multiply(*matCRS,false,vCRS,false,tCRS);

     // now getrow tCRS and compare with y 
     Count = 0; 
     for (int i = 0; i < rowMap.getNodeNumElements(); i++) {
        tCRS.ExtractMyRowCopy(i, colMap.getNodeNumElements(), RowLength, values, indices);
        if (RowLength == 0) {
           RowLength = 1; values[0] = 0.; indices[0] = 0;
         }
        if (RowLength != 1) {
          std::cout << std::setw(5) << myPID << ":matvec/matmult row (GID=" << 
          std::setw(6) << rowMap.GID(i) << ",LID=" << std::setw(6) << i << 
          ") rowlength problem" << std::endl;
        }
        if (indices[0] != 0) {
          std::cout << std::setw(5) << myPID << ":matvec/matmult row (GID=" << 
          std::setw(6) << rowMap.GID(i) << ",LID=" << std::setw(6) << i << 
          ") has Col Id = " << indices[0] << std::endl;
        }
        t[i] = values[0];
     }
     // might need to do some import/export stuff if RangeMap and rowMap 
     // do not match. This all depends on what matmat mult does?

     // print if abs(y[i] - t[i])/z[i] > 1.e-8

     for (int i = 0; i < rowMap.getNodeNumElements(); i++) {
        if (( fabs(y[i] - t[i])/z[i] > 1.e-8) && (Count++ < 20) ) {
          std::cout << std::setw(5) << myPID << ":matvec/matmult mismatch,row (GID=" << 
          std::setw(6) << rowMap.GID(i) << ",LID=" << std::setw(6) << i << 
          "), y t & z =" << std::setw(9) << std::scientific << 
          std::setprecision(2) << y[i] << " " << std::setw(9) << values[0] 
          << " " << std::setw(9) << z[i] << std::endl;
        }
     } // for (int i = 0; ...
  } // if (matCRS)
  else
    if (myPID == 0) std::cout << "not crs" << std::endl;

#endif

  return 1;

} // MatrixVectorChecker()

} // namespace MueLu

#endif

// JG TODO:
// - Add the test with matrix-matrix multiply.
// - Do we need Import/Export operation ?
// - Can we always use CrsMatrix->getLocalRowView() (on GPU etc.) ?
// - Can we always use Vector->getDataNonConst() (on GPU etc.) ?
